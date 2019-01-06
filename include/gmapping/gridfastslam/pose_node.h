
#ifndef _POSE_NODE_H_
#define _POSE_NODE_H_

#include <vector>
#include <deque>
#include <gmapping/utils/point.h>
#include <gmapping/sensor/sensor_range/rangereading.h>

namespace GMapping {

  /**This class defines the the node of reversed tree in which the trajectories are stored.
     Each node of a tree has a pointer to its parent and a counter indicating the number of childs of a node.
      The tree is updated in a way consistent with the operation performed on the particles.
  */
  struct TNode{
    /**Constructs a node of the trajectory tree.
     @param pose:      the pose of the robot in the trajectory
      @param weight:    the weight of the particle at that point in the trajectory
      @param accWeight: the cumulative weight of the particle
      @param parent:    the parent node in the tree
      @param childs:    the number of childs
    */
    TNode(const OrientedPoint& pose, double weight, TNode* parent=0, unsigned int childs=0);   

    /**Destroys a tree node, and consistently updates the tree. If a node whose parent has only one child is deleted,
     also the parent node is deleted. This because the parent will not be reacheable anymore in the trajectory tree.*/
    ~TNode();
    /**
    {
      if (parent && (--parent->childs)<=0)
        delete parent;
       assert(!childs);
    }
    **/

    /**The pose of the robot*/
    OrientedPoint pose; 
    
    /**The weight of the particle*/
    double weight;

    /**The sum of all the particle weights in the previous part of the trajectory*/
    double accWeight;

    /**The parent*/
    TNode* parent;

    /**The range reading to which this node is associated*/
    const RangeReading* reading;

    /**The number of childs*/
    unsigned int childs;

    /**counter in visiting the node (internally used)*/
    mutable unsigned int visitCounter;

    /**visit flag (internally used)*/
    mutable bool flag;
  };

  typedef std::vector<TNode*> TNodeVector;
  typedef std::deque<TNode*> TNodeDeque;

};

#endif