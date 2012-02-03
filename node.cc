#include "node.hh"
#include "conf.hh"

#include "cart.hh"

Node::Node(Conf &conf): conf_(conf) {
  switch (conf_.node_type) {
    case Conf::CART_NODE: work_node_ = new CartNode(conf_); break;
    case Conf::IDLE_NODE: /* fall down */
    default: work_node_ = new WorkNode;
  }
}

Node::~Node() {
  if (work_node_) delete work_node_;
}
