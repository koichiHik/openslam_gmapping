
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <vector>
#include <deque>
#include <gmapping/scanmatcher/smmap.h>
#include <gmapping/gridfastslam/pose_node.h>

namespace GMapping {

  struct Particle{
    /**constructs a particle, given a map */
    Particle(const ScanMatcherMap& map);

    /** @returns the weight of a particle */
    inline operator double() const {return weight;}

    /** @returns the pose of a particle */    
    inline operator OrientedPoint() const {return pose;}

    /** Sets the weight of the particle. */
    inline void setWeight(double w) {weight=w;}
    /** The map */
    ScanMatcherMap map;
    /** The pose of the robot */
    OrientedPoint pose;
    /** The pose of the robot at the previous time frame (used for computing thr odometry displacements) */
    OrientedPoint previousPose;

    /** The weight of the particle */
    double weight;

    /** The cumulative weight of the particle */
    double weightSum;

    /** Entry to the trajectory tree */
    TNode* node; 
  };

  typedef std::vector<Particle> ParticleVector;

};

#endif