#include "cart.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <ostream>
#include "conf.hh"
#include "input.hh"
#include "output.hh"
#include "ptcl.hh"
#include "random.hh"

CartNode::CartNode(Conf &conf): conf_(conf) {
  using namespace std;
  // setup timer
  t_cfrc.Label("cart force").Comm(conf_.cart_comm);
  t_comm.Label("cart comm").Comm(conf_.cart_comm);
  t_step.Label("cart step").Comm(conf_.cart_comm);
  t_oput.Label("cart output").Comm(conf_.cart_comm);
  t_init.Label("cart init").Comm(conf_.cart_comm).Start();

  // make local copy
  cart_comm = conf_.cart_comm;
  sys_size = conf_.sys_size;
  sys_size_2 = 0.5 * sys_size;
  sys_min = conf_.sys_min;
  sys_max = conf_.sys_max;
  dt = conf_.delta_t;
  cutoff2 = conf_.cutoff * conf_.cutoff;

  MPI_Comm_rank(cart_comm, &cart_rank);
  MPI_Comm_size(cart_comm, &cart_size);
  MPI_Cart_coords(cart_comm, cart_rank, 3, cart_pos);

  periodic = conf_.periodic;

  if (cart_rank == 0) {
    if (conf_.verbose > 0)
      cout << setw(24) << left << "# cart_size" << cart_size << "\n";

    // prepare input stream
    if (!conf_.ifname.empty()) {
      is.open(conf_.ifname.c_str());
      if (!is) {
        cout << "pmt: cannot open input file '" << conf_.ifname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }

    // prepare output stream
    if (!conf_.ofname.empty()) {
      os.open(conf_.ofname.c_str());
      if (!os) {
        cout << "pmt: cannot open output file '" << conf_.ofname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }

    // prepare restart save stream
    if (!conf_.rfname.empty()) {
      rs.open(conf_.rfname.c_str());
      if (!rs) {
        cout << "pmt: cannot open restart save file '" << conf_.rfname
            << "'. abort\n" << flush;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }
  }

  // setup node
  div_min = div_max = sys_min;
  div_min += v3r(cart_pos    ) * sys_size / v3r(conf_.cart_num);
  div_max += v3r(cart_pos + 1) * sys_size / v3r(conf_.cart_num);

  // topology
  InitConnect();

  // generate particles
  GenerateParticles();

  // write particle coordinates at step 0
  if (!conf_.ofname.empty() || !conf_.rfname.empty()) {
    t_oput.Start();
    if (conf_.write_step0 || ((conf_.max_step == 0) && rs)) {
      OutputXYZ(os, rs, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
                0, conf_.max_step, cart_rank, cart_size, cart_comm);
    }
    steps_to_write = conf_.write_interval;
    t_oput.Stop();
  }

  // calculate initial force
  if (conf_.max_step > 0) {
    t_cfrc.Start();
    CalculateForce();
    t_cfrc.Stop();
  }

  t_init.Stop();
}

CartNode::~CartNode() {
  // print timer
  if (conf_.verbose > 1) t_cfrc.PrintAll("# ");
  if (conf_.verbose > 0) t_cfrc.PrintMax("# max ");
  if (conf_.verbose > 1) t_comm.PrintAll("# ");
  if (conf_.verbose > 0) t_comm.PrintMax("# max ");
  if (conf_.verbose > 1) t_step.PrintAll("# ");
  if (conf_.verbose > 0) t_step.PrintMax("# max ");
  if (conf_.verbose > 1) t_oput.PrintAll("# ");
  if (conf_.verbose > 0) t_oput.PrintMax("# max ");
  if (conf_.verbose > 1) t_init.PrintAll("# ");
  if (conf_.verbose > 0) t_init.PrintMax("# max ");
}

void CartNode::StepForward(int t) {
  using namespace std;
  t_step.Start();
  Ptcls::size_type p_size = ptcls.size();

  // update positions and velocities
  // v(t+0.5dt) = v(t) + 0.5dt*Fi(t)/mi
  // r(t+dt) = r(t) + dt*v(t+0.5dt)
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    Ptcl &ptcl = ptcls[p];
    ptcl.attr = 0;
    ptcl.crd += dt * (ptcl.vel += (force[p] * ptcl.inv_2mass));

    // embarkation check & periodic shift
    for (int i = 0; i < 3; ++i) {
      ptcl.attr <<= 2;
      if (ptcl.crd[i] < div_min[i]) {
        switch (periodic[i]) {
          case 2: {
            if (ptcl.crd[i] < sys_min[i]) {
              ptcl.crd[i] = 2.0 * sys_min[i] - ptcl.crd[i];
              ptcl.vel[i] *= -1.0;
            }
            else ptcl.attr |= LOWER_DIR;
            break;
          }
          case 1: {
            if (ptcl.crd[i] < sys_min[i]) ptcl.crd[i] += sys_size[i];
            ptcl.attr |= LOWER_DIR;
            break;
          }
          case 0:
          default: {
            if (ptcl.crd[i] >= sys_min[i]) ptcl.attr |= LOWER_DIR;
          }
        }
      }
      else if (ptcl.crd[i] >= div_max[i]) {
        switch (periodic[i]) {
          case 2: {
            if (ptcl.crd[i] >= sys_max[i]) {
              ptcl.crd[i] = 2.0 * sys_max[i] - ptcl.crd[i];
              ptcl.vel[i] *= -1.0;
            }
            else ptcl.attr |= UPPER_DIR;
            break;
          }
          case 1: {
            if (ptcl.crd[i] >= sys_max[i]) ptcl.crd[i] -= sys_size[i];
            ptcl.attr |= UPPER_DIR;
            break;
          }
          case 0:
          default: {
            if (ptcl.crd[i] < sys_max[i]) ptcl.attr |= UPPER_DIR;
          }
        }
      }
    }
  }

  // calculate force F(t+0.5dt) using r(t+dt)
  // F(t+0.5dt)
  t_cfrc.Start();
  CalculateForce();
  t_cfrc.Stop();

  // update velocities
  // v(t+dt) = v(t+0.5dt) + 0.5dt*Fi(t+dt)/mi
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    Ptcl &ptcl = ptcls[p];
    ptcl.vel += force[p] * ptcl.inv_2mass;
  }

  t_comm.Start();
  ExchangeParticles();
  t_comm.Stop();
  t_step.Stop();

  // write particle coordinates at step t
  if (!conf_.ofname.empty() || !conf_.rfname.empty()) {
    t_oput.Start();
    if ((steps_to_write <= 1) || ((conf_.max_step == t) && rs)) {
      OutputXYZ(os, rs, conf_.cmd_line.c_str(), ptcls, conf_.total_ptcl,
                t, conf_.max_step, cart_rank, cart_size, cart_comm);
      steps_to_write = conf_.write_interval;
    } else {
      --steps_to_write;
    }
    t_oput.Stop();
  }
}

void CartNode::InitConnect() {
  const unsigned tag[3] = { MIDDLE_DIR, LOWER_DIR, UPPER_DIR };
  const int offset[3] = { 0, -1, 1 };
  v3i coords;
  conns.reserve(26);
  for (int x = 0; x < 3; ++x) {
    for (int y = 0; y < 3; ++y) {
      for (int z = 0; z < 3; ++z) {
        if (x == 0 && y == 0 && z == 0) continue;
        Connect conn;
        conn.dir_tag = ((tag[x] << 4) | (tag[y] << 2) | tag[z]);
        //
        coords = cart_pos;
        coords[0] += offset[x];
        coords[1] += offset[y];
        coords[2] += offset[z];
        for (int i = 0; i < 3; ++i) {
          if (periodic[i] != 1) {  // non-periodic
            if (coords[i] < 0) {
              coords[i] = 0;
            } else if (coords[i] >= conf_.cart_num[i]) {
              coords[i] = conf_.cart_num[i] - 1;
            }
          }
        }
        MPI_Cart_rank(cart_comm, coords, &conn.send_to);
        //
        coords = cart_pos;
        coords[0] -= offset[x];
        coords[1] -= offset[y];
        coords[2] -= offset[z];
        for (int i = 0; i < 3; ++i) {
          if (periodic[i] != 1) {  // non-periodic
            if (coords[i] < 0) {
              coords[i] = 0;
            } else if (coords[i] >= conf_.cart_num[i]) {
              coords[i] = conf_.cart_num[i] - 1;
            }
          }
        }
        MPI_Cart_rank(cart_comm, coords, &conn.recv_from);
        //
        conns.push_back(conn);
      }
    }
  }
}

void CartNode::GenerateParticles() {
  using namespace std;
  if (!conf_.ifname.empty()) {
    // read from XYZ input file
    InputXYZ(is, &ptcls, &conf_.total_ptcl, div_min, div_max,
             cart_rank, cart_size, cart_comm);
  } else {
    // generate at random
    real_t chg = 0.0;
    const int last_id = conf_.total_ptcl - 1;
    for (int n = 0; n < conf_.total_ptcl; ++n) {
      Ptcl p;
      for (int i = 0; i < 3; ++i) {
        p.crd[i] = Rand() * sys_size[i] + sys_min[i];
        p.vel[i] = Gaussian(0.0, 0.5);
      }

      if (conf_.uniformize) {
        if (conf_.neutralize) {
          chg = ((n & 1) == 0)
              ? ((n < last_id) ? 20.0 : 0.0)
              : -chg;
        } else {
          chg = 20.0;
        }
        p.chg = chg;
        p.inv_2mass = 0.05;
      } else {
        if (conf_.neutralize) {
          chg = ((n & 1) == 0)
              ? ((n < last_id) ? (Gaussian(0.0, 0.2) * 100) : 0.0)
              : -chg;
        } else {
          chg = Gaussian(0.0, 0.2) * 100;
        }
        p.chg = chg;
        p.inv_2mass = 0.5 / (Rand() * 15.0 + 1.0);
      }

      if (div_min <= p.crd && p.crd < div_max) {
        p.id = n;
        ptcls.push_back(p);
      }
    }
  }
  force.resize(ptcls.size());

  if (cart_rank == 0) {
    if (conf_.verbose > 0)
      cout << setw(24) << left << "# total_ptcl" << conf_.total_ptcl << "\n";
  }
}

void CartNode::ExchangeParticles() {
  using namespace std;
  sort(ptcls.begin(), ptcls.end(), Less<Ptcl::DIR>());

  // find lower bound of outgoing ptcl
  Ptcl stayin;
  stayin.attr = 0x00;
  Ptcls::iterator lb_outgoing =
      upper_bound(ptcls.begin(), ptcls.end(), stayin, Less<Ptcl::DIR>());

  // copy outgoing ptcls
  Ptcls send_buff(lb_outgoing, ptcls.end());
  // trim ptcls
  ptcls.erase(lb_outgoing, ptcls.end());

  for (int i = 0; i < 26; ++i) {
    // get range of ptcl whose attr equals a dir_tag
    Ptcl p;
    p.attr = conns[i].dir_tag;
    pair<Ptcls::iterator, Ptcls::iterator> range =
        equal_range(send_buff.begin(), send_buff.end(), p, Less<Ptcl::DIR>());

    // get offsets from ptcls[0] and number of Ptcl to send
    int first = static_cast<int>(distance(send_buff.begin(), range.first));
    int count = static_cast<int>(distance(range.first, range.second));
    int bytes = count * sizeof(Ptcl);

    // send immediately
    MPI_Isend(&send_buff[first], bytes, MPI_BYTE, conns[i].send_to,
              static_cast<int>(p.attr), cart_comm, &conns[i].req);
  }

  for (int i = 0; i < 26; ++i) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, cart_comm, &status);

    // get size of message
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    int count = bytes / sizeof(Ptcl);

    // receive using i-th buffer
    recv_buff[i].resize(count);
    MPI_Recv(&recv_buff[i][0], bytes, MPI_BYTE, status.MPI_SOURCE,
             status.MPI_TAG, cart_comm, MPI_STATUS_IGNORE);

    // append incoming ptcls to ptcls
    ptcls.insert(ptcls.end(), recv_buff[i].begin(), recv_buff[i].end());
  }

  // wait all send
  for (int i = 0; i < 26; ++i) {
    MPI_Wait(&conns[i].req, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(cart_comm);

  force.resize(ptcls.size());
}

void CartNode::CalculateForce() {
#if 1
  CalculateForceCutoffPeriodic();
#else
  force.assign(ptcls.size(), 0.0);
#endif
}

void CartNode::CalculateForceCutoffPeriodic() {
  const double kInv4Pi = 0.25 / acos(-1.0);
  Ptcls::size_type p_size = ptcls.size();

  // clear force
  force.assign(p_size, 0.0);

  // prepare recieve buffer
  std::vector<v4r> j_ptcls(conf_.total_ptcl);

  // pack my particles in j_ptcls
#pragma omp parallel for schedule(static)
  for (Ptcls::size_type p = 0; p < p_size; ++p) {
    j_ptcls[p][0] = ptcls[p].crd[0];
    j_ptcls[p][1] = ptcls[p].crd[1];
    j_ptcls[p][2] = ptcls[p].crd[2];
    j_ptcls[p][3] = ptcls[p].chg;
  }

  Ptcls::size_type j_size = p_size;

  // calculate force from particles of local node
#pragma omp parallel for schedule(dynamic, 8)
  for (Ptcls::size_type i = 0; i < p_size; ++i) {
    for (Ptcls::size_type j = 0; j < i; ++j) {
      v3r r;
      for (int d = 0; d < 3; ++d) {
        r[d] = ptcls[i].crd[d] - j_ptcls[j][d];
        if (periodic[d] == 1) {
          if      (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
          else if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
        }
      }
      double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      if (r2 < cutoff2) {
        double r6 = r2 * r2 * r2;
        double _r3 = 1.0 / sqrt(r6);  // divzero ?
        double chgs = ptcls[i].chg * j_ptcls[j][3];
        r *= (kInv4Pi * chgs) * _r3;
        force[i] += r;
        force[j] -= r;
      }
    }
  }

  // prepare send buffer
  std::vector<v4r> send_buffer;

  // calculate force from particles of external node
  for (int i = 1; i < cart_size; ++i) {
    // copy j_ptcls to send_buffer
    send_buffer.resize(j_size);
    for (Ptcls::size_type j = 0; j < j_size; ++j) {
      send_buffer[j] = j_ptcls[j];
    }

    // send
    int dest = (cart_rank + 1) % cart_size;
    MPI_Request req;
    MPI_Isend(&send_buffer[0], j_size * sizeof(v4r), MPI_BYTE,
              dest, cart_rank, cart_comm, &req);

    // recv
    int source = (cart_rank + cart_size - 1) % cart_size;
    MPI_Status status;
    MPI_Recv(&j_ptcls[0], conf_.total_ptcl * sizeof(v4r), MPI_BYTE,
             source, source, cart_comm, &status);

    // wait
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    // get size of receive message
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    j_size = bytes / sizeof(v4r);

    // calculate force
#pragma omp parallel for collapse(2), schedule(dynamic, 4)
    for (Ptcls::size_type i = 0; i < p_size; ++i) {
      for (Ptcls::size_type j = 0; j < j_size; ++j) {
        v3r r;
        for (int d = 0; d < 3; ++d) {
          r[d] = ptcls[i].crd[d] - j_ptcls[j][d];
          if (periodic[d] == 1) {
            if      (r[d] >  sys_size_2[d]) r[d] -= sys_size[d];
            else if (r[d] < -sys_size_2[d]) r[d] += sys_size[d];
          }
        }
        double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        if (r2 < cutoff2) {
          double r6 = r2 * r2 * r2;
          double _r3 = 1.0 / sqrt(r6);  // divzero ?
          double chgs = ptcls[i].chg * j_ptcls[j][3];
          r *= (kInv4Pi * chgs) * _r3;
          force[i] += r;
        }
      }
    }
  }
}
