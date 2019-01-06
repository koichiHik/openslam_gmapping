
#include <gmapping/gridfastslam/pose_node.h>

namespace GMapping {

  TNode::TNode(const OrientedPoint& p, double w, TNode* n, unsigned int c) :
    pose(p), weight(w), childs(c), parent(n), reading(0), flag(0), accWeight(0)
  {
    if (n){
      n->childs++;
    }
  }

  TNode::~TNode(){
    if (parent && (--parent->childs)<=0)
      delete parent;
    assert(!childs);
  }

}