
#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif

/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/
inline void GridSlamProcessor::scanMatch(const double* plainReading){

  // Sample a new pose from each scan in the reference
  double sumScore=0;

  // Iterate for all particles existing.
  for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
    OrientedPoint corrected;
    double score, l, s;

    // Improve pose of the particle based on scan matching.
    score = m_matcher.optimize(corrected, it->map, it->pose, plainReading);

    // If the score gets better, replace pose.
    if (score > m_params.gridSlamProcParams.minimumScore){
      it->pose = corrected;

    } else if (m_infoStream) {
      m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l <<std::endl;
      m_infoStream << "lp:" << m_lastPartPose.x << " "  << m_lastPartPose.y << " "<< m_lastPartPose.theta <<std::endl;
      m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " "<< m_odoPose.theta <<std::endl;
    }

    // Calculate score for weighting calculation.
    // This weighting is for particle.
    m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);
    sumScore += score;
    it->weight += l;
    it->weightSum += l;

    // Set up the selective copy of the active area
    // by detaching the areas that will be updated
    m_matcher.invalidateActiveArea();
    m_matcher.computeActiveArea(it->map, it->pose, plainReading);
  }

  if (m_infoStream) {
    m_infoStream << "Average Scan Matching Score=" << sumScore/m_particles.size() << std::endl;
  }
}

inline void GridSlamProcessor::normalize(){

  // Search maximum weight of the particle.
  double lmax = -std::numeric_limits<double>::max();
  for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
    lmax = std::max(it->weight, lmax);
  }
  
  // Create weights vector that is normalized.
  m_weights.clear();
  double wcum = 0;
  double gain = 1.0 / (m_params.gridSlamProcParams.obsSigmaGain * m_particles.size());
  for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++){
    m_weights.push_back(exp(gain*(it->weight - lmax)));
    wcum += m_weights.back();
  }
  
  m_neff = 0;
  for (std::vector<double>::iterator it = m_weights.begin(); it != m_weights.end(); it++){
    *it = *it / wcum;
    double w = *it;
    m_neff += w * w;
  }
  m_neff = 1.0 / m_neff;
}

inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading){
  
  bool hasResampled = false;

  // If "m_neff" gets smaller than the threshold, start resampling procedure.
  if (m_neff < m_params.gridSlamProcParams.resampleThreshold * m_particles.size()){
    
    if (m_infoStream) {
      m_infoStream  << "*************RESAMPLE***************" << std::endl;
    }
    
    // Use low-variance sampler for resampling.
    // m_weights contains normalized weight of each particle.
    std::vector<u_int32_t> indexes;
    uniform_resampler<double, double> resampler;
    indexes = resampler.resampleIndexes(m_weights, adaptSize);
    
    ParticleVector newlyGeneratedParticles;
    uint32_t j = 0;

    // This is for deleteing the particles which have been resampled away.
    std::vector<uint32_t> deletedParticles;
    
    // Stroing index of deleted particles.
    for (uint32_t i = 0; i < indexes.size(); i++){

      while(j < indexes[i]){ 
        deletedParticles.push_back(j);
        j++;
			}
      if (j == indexes[i]) {
	      j++;
      }

      // Selected particle.
      Particle& p = m_particles[indexes[i]];
      newlyGeneratedParticles.push_back(p);
      newlyGeneratedParticles.back().setWeight(0);
    }

    while(j < indexes.size()){
      deletedParticles.push_back(j);
      j++;
    }

    std::cerr << "Deleting Nodes:";
    for (unsigned int i = 0; i < deletedParticles.size(); i++){
      std::cerr << " " << deletedParticles[i];
      // Delete trajectory nodes because its on heap.
      delete m_particles[deletedParticles[i]].node;
      m_particles[deletedParticles[i]].node = 0;
    }

    // Re-generation of m_particles.
    m_particles.clear();
    for (ParticleVector::iterator itr=newlyGeneratedParticles.begin(); 
          itr!= newlyGeneratedParticles.end(); itr++) {
      m_particles.push_back(*itr);
    }
    //std::copy(newlyGeneratedParticles.begin(), newlyGeneratedParticles.end(), m_particles.begin());
    hasResampled = true;
  }

  for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++) {
    // Register scan for the map owned by THIS particle.
    TNode* node = new TNode(it->pose, 0.0, it->node, 0);
    node->reading = reading;
    it->node = node;
    m_matcher.invalidateActiveArea();
    m_matcher.registerScan(it->map, it->pose, plainReading);
  } 

  return hasResampled;
}
