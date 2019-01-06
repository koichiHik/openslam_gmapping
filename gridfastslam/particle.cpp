
#include <gmapping/gridfastslam/particle.h>

namespace GMapping {

  //HERE STARTS THE BEEF
  Particle::Particle(const ScanMatcherMap& m)
  : map(m), pose(0,0,0), weight(0), weightSum(0) {
  node = 0;
  }

};