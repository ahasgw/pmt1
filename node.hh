#ifndef NODE_HH_
#define NODE_HH_ 1

class Conf;

class WorkNode {
 public:
  virtual ~WorkNode() {}
  virtual void StepForward(int t) {}
};

class Node {
  Conf &conf_;
  WorkNode *work_node_;

 public:
  Node(Conf &conf);
  ~Node();

  void StepForward(int t) { work_node_->StepForward(t); }
};

#endif  // NODE_HH_
