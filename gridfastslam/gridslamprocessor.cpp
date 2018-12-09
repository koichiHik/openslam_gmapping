#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include <gmapping/gridfastslam/gridslamprocessor.h>

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping {

const double m_distanceThresholdCheck = 20;
 
using namespace std;

  GridSlamProcessor::GridSlamProcessor()
    : m_infoStream(cout) {}

  GridSlamProcessor::GridSlamProcessor(std::ostream& infoS)
    : m_infoStream(infoS){}
  
  GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor& gsp) 
    :last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout){
    
    m_params = gsp.m_params;

    
    m_beams=gsp.m_beams;
    m_indexes=gsp.m_indexes;
    m_motionModel=gsp.m_motionModel;
    m_matcher=gsp.m_matcher;
    
    m_count=gsp.m_count;
    m_readingCount=gsp.m_readingCount;
    m_lastPartPose=gsp.m_lastPartPose;
    m_pose=gsp.m_pose;
    m_odoPose=gsp.m_odoPose;
    m_linearDistance=gsp.m_linearDistance;
    m_angularDistance=gsp.m_angularDistance;
    m_neff=gsp.m_neff;
		
    m_xmin=gsp.m_xmin;
    m_ymin=gsp.m_ymin;
    m_xmax=gsp.m_xmax;
    m_ymax=gsp.m_ymax;
    m_delta=gsp.m_delta;  
    
    TNodeVector v=gsp.getTrajectories();
    for (unsigned int i=0; i<v.size(); i++){
		m_particles[i].node=v[i];
    }

    cerr  << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ..." ;
    updateTreeWeights(false);
    cerr  << ".done!" <<endl;
  }
  

  GridSlamProcessor* GridSlamProcessor::clone() const {
  	GridSlamProcessor* cloned=new GridSlamProcessor(*this);
    return cloned;
  }
  
  GridSlamProcessor::~GridSlamProcessor(){
    cerr << __PRETTY_FUNCTION__ << ": Start" << endl;
    cerr << __PRETTY_FUNCTION__ << ": Deleting tree" << endl;
    for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++){
      if (it->node) {
	      delete it->node;
      }
    }   
  }

  //HERE STARTS THE BEEF
  GridSlamProcessor::Particle::Particle(const ScanMatcherMap& m)
    : map(m), pose(0,0,0), weight(0), weightSum(0), gweight(0), previousIndex(0){
    node = 0;
  }
   
  void GridSlamProcessor::setSensorMap(const SensorMap& smap){
    
    /*
      Construct the angle table for the sensor
      FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot 
    */
    
    SensorMap::const_iterator laser_it=smap.find(std::string("FLASER"));
    if (laser_it==smap.end()){
      cerr << "Attempting to load the new carmen log format" << endl;
      laser_it=smap.find(std::string("ROBOTLASER1"));
      assert(laser_it!=smap.end());
    }
    const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>((laser_it->second));
    assert(rangeSensor && rangeSensor->beams().size());
    
    m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
    double* angles=new double[rangeSensor->beams().size()];
    for (unsigned int i=0; i<m_beams; i++){
      angles[i]=rangeSensor->beams()[i].pose.theta;
    }
    m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
    delete [] angles;
  }
  
  void GridSlamProcessor::init(const GridSlamProcessor::Params& params, OrientedPoint initialPose) {
    
    m_params = params;
    m_matcher.setMatchingParameters(params.scanMatcherParams);
    m_motionModel.setMotionModelParams(params.motionModelParams);

    m_xmin = m_params.gridSlamProcParams.x_min;
    m_xmax = m_params.gridSlamProcParams.x_max;
    m_ymin = m_params.gridSlamProcParams.y_min;
    m_ymax = m_params.gridSlamProcParams.y_max;
    m_delta = m_params.gridSlamProcParams.delta;

    if (m_infoStream) {
      m_infoStream 
      << " -xmin " << m_xmin 
      << " -xmax "<< m_xmax 
      << " -ymin "<< m_ymin 
      << " -ymax "<< m_ymax	
      << " -delta "<< m_delta
      << " -particles "<< m_params.gridSlamProcParams.particle_num << endl;
    }   

    // Initialization of Particles.
    m_particles.clear();
    TNode* node = new TNode(initialPose, 0, 0, 0);
    ScanMatcherMap lmap(Point(m_xmin + m_xmax, m_ymin + m_ymax)*.5, m_xmax - m_xmin, m_ymax - m_ymin, m_delta);
    for (unsigned int i = 0; i < m_params.gridSlamProcParams.particle_num; i++) {
      m_particles.push_back(Particle(lmap));
      m_particles.back().pose = initialPose;
      m_particles.back().previousPose = initialPose;
      m_particles.back().setWeight(0);
      m_particles.back().previousIndex = 0;

  		// we use the root directly
	  	m_particles.back().node= node;
    }

    m_neff = (double)m_params.gridSlamProcParams.particle_num;
    m_count = 0;
    m_readingCount = 0;
    m_linearDistance = m_angularDistance=0;

    m_infoStream << "m_particles.size() at init" << m_particles.size() << " : " << m_neff << endl;

  }

  void GridSlamProcessor::processTruePos(const OdometryReading& o){

    const OdometrySensor* os=dynamic_cast<const OdometrySensor*>(o.getSensor());
    if (os && os->isIdeal() && m_outputStream){
      m_outputStream << setiosflags(ios::fixed) << setprecision(3);
      m_outputStream <<  "SIMULATOR_POS " <<  o.getPose().x << " " << o.getPose().y << " " ;
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " <<  o.getTime() << endl;
    }
  }

  bool GridSlamProcessor::processScan(const RangeReading & reading, int adaptParticles){

    // All reading is with pose. Retireve the position from the reading, and compute the odometry.
    OrientedPoint relPose = reading.getPose();
    if (!m_count){
      m_lastPartPose = m_odoPose = relPose;
    }
    
    // Generate the next generation of the particles.
    for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++){
      OrientedPoint& pose(it->pose);
      pose = m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);
    }

    // Output the updated particle pose.
    printUpdatedParticlePose(reading);
    
    // Invoke the callback.
    onOdometryUpdate();   

    // Accumulate the robot translation and rotation
    OrientedPoint move = relPose - m_odoPose;
    move.theta = atan2(sin(move.theta), cos(move.theta));
    m_linearDistance += sqrt(move*move);
    m_angularDistance += fabs(move.theta);
    
    // If the robot moves more than threshold, raise warning message.
    if (m_linearDistance > m_distanceThresholdCheck){
      printOdomJumpWarnings(relPose);
    }

    // Update the last odometry.
    m_odoPose=relPose;
    bool processed=false;

    // Process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
    if (!m_count 
      || m_linearDistance >= m_params.gridSlamProcParams.linearThresholdDistance
      || m_angularDistance >= m_params.gridSlamProcParams.angularThresholdDistance
      || (0.0 <= m_params.gridSlamProcParams.period && m_params.gridSlamProcParams.period < (reading.getTime() - last_update_time_))) {
      
      last_update_time_ = reading.getTime();      

      printBasicInfoForScanUpdate();
      
      printScanFrameInfo(reading);

      // This is for converting the reading in a scan-matcher feedable form
      assert(reading.size() == m_beams);
      double* plainReading = new double[m_beams];
      for(unsigned int i = 0; i < m_beams; i++) {
      	plainReading[i] = reading[i];
      }
      
      RangeReading* reading_copy = 
          new RangeReading( reading.size(), 
                            &(reading[0]),
                            static_cast<const RangeSensor*>(reading.getSensor()),
                            reading.getTime());

      if (m_count > 0) {

        // Update likelihood of each particles by matching the current scan.
        // Update generative model
        m_infoStream << "Before Scan Match"<< endl;
        scanMatch(plainReading);
        m_infoStream << "After Scan Match"<< endl;

        if (m_outputStream.is_open()){
          printScanMatchUpdateInfo(reading);
        }

        // Invoke callback.
        onScanmatchUpdate();
        
        updateTreeWeights(false);

        printUpdatedNeff();              

        // Resample Particle
        // TODO Scan Registration is currently called inside. Take outside.
        resample(plainReading, adaptParticles, reading_copy);
        m_infoStream << "After resample"<< endl;
	
      } else {

        m_infoStream << "Registering First Scan"<< endl;
        for (ParticleVector::iterator it = m_particles.begin(); it != m_particles.end(); it++){

          // Create a new node for each particle. This reading will be treated as "landmark".
          // This node has its own scan.
          TNode* node = new TNode(it->pose, 0.0, it->node, 0);
          node->reading = reading_copy;
          it->node = node;

          m_infoStream << "invalidateActiveArea"<< endl;
          m_matcher.invalidateActiveArea();
          m_infoStream << "computeActiveArea"<< endl;
          m_matcher.computeActiveArea(it->map, it->pose, plainReading);
          m_infoStream << "registerScan"<< endl;
          m_matcher.registerScan(it->map, it->pose, plainReading);          
        }

        m_infoStream << "Out side of loop."<< endl;

      }

      updateTreeWeights(false);
      
      m_infoStream << "After Update Tree Weight"<< endl;

      delete [] plainReading;
      // Update the past pose for the next iteration.
      m_lastPartPose = m_odoPose;
      m_linearDistance = 0;
      m_angularDistance = 0;
      m_count++;
      processed = true;
      
      // Preparation for next step. Update the past pose for all existing particles.
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	      it->previousPose = it->pose;
      }
      
    }

    printFlush();

    m_infoStream << "End of processScan"<< endl;

    m_readingCount++;
    return processed;
  }
  
  std::ofstream& GridSlamProcessor::outputStream(){
    return m_outputStream;
  }
  
  std::ostream& GridSlamProcessor::infoStream(){
    return m_infoStream;
  }

  void GridSlamProcessor::printOdomJumpWarnings(const GMapping::OrientedPoint & new_pose) {

    cerr << "***********************************************************************" << endl;
    cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
    cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
    cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y 
    << " " <<m_odoPose.theta << endl;
    cerr << "New Odometry Pose (reported from observation)= " << new_pose.x << " " << new_pose.y 
    << " " <<new_pose.theta << endl;
    cerr << "***********************************************************************" << endl;
    cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
    cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
    cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
    cerr << "***********************************************************************" << endl;
    
  }

  void GridSlamProcessor::printBasicInfoForScanUpdate() {
    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "FRAME " <<  m_readingCount;
      m_outputStream << " " << m_linearDistance;
      m_outputStream << " " << m_angularDistance << endl;
    }
  }


  void GridSlamProcessor::printUpdatedParticlePose(const RangeReading & reading) {

    if (m_outputStream.is_open()) {
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODOM ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(3) << m_odoPose.x << " " << m_odoPose.y << " ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << m_odoPose.theta << " ";
      m_outputStream << reading.getTime();
      m_outputStream << endl;
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODO_UPDATE "<< m_particles.size() << " ";
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
        OrientedPoint& pose(it->pose);
        m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
        m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it-> weight << " ";
      }
      m_outputStream << reading.getTime();
      m_outputStream << endl;
    }

  }
  
  void GridSlamProcessor::printUpdatedNeff() {
    if (m_infoStream){
      m_infoStream << "neff= " << m_neff  << endl;
    }

    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "NEFF " << m_neff << endl;
    }
  }

  void GridSlamProcessor::printScanMatchUpdateInfo(const RangeReading & reading) {
    if (m_outputStream.is_open()){
      m_outputStream << "LASER_READING "<< reading.size() << " ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(2);
      for (RangeReading::const_iterator b=reading.begin(); b!=reading.end(); b++){
        m_outputStream << *b << " ";
      }
      OrientedPoint p=reading.getPose();
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << p.x << " " << p.y << " " << p.theta << " " << reading.getTime()<< endl;
      m_outputStream << "SM_UPDATE "<< m_particles.size() << " ";
      for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
        const OrientedPoint& pose=it->pose;
        m_outputStream << setiosflags(ios::fixed) << setprecision(3) <<  pose.x << " " << pose.y << " ";
        m_outputStream << setiosflags(ios::fixed) << setprecision(6) <<  pose.theta << " " << it-> weight << " ";
      }
      m_outputStream << endl;
    }
  }

  void GridSlamProcessor::printScanFrameInfo(const RangeReading & reading) {

    if (m_infoStream) {
      m_infoStream << "update frame " <<  m_readingCount << endl
      << "update ld=" << m_linearDistance 
      << " ad=" << m_angularDistance << endl
      << "Laser Pose= " << reading.getPose().x << " " << reading.getPose().y
      << " " << reading.getPose().theta << endl
      << "m_count " << m_count << endl;
    }

  }

  void GridSlamProcessor::printFlush() {
    if (m_outputStream.is_open()) {
      m_outputStream << flush;
    }
  }
  
  int GridSlamProcessor::getBestParticleIndex() const{
    unsigned int bi=0;
    double bw=-std::numeric_limits<double>::max();
    for (unsigned int i=0; i<m_particles.size(); i++)
      if (bw<m_particles[i].weightSum){
	bw=m_particles[i].weightSum;
	bi=i;
      }
    return (int) bi;
  }

  void GridSlamProcessor::onScanmatchUpdate(){}
  void GridSlamProcessor::onResampleUpdate(){}
  void GridSlamProcessor::onOdometryUpdate(){}

  
};// end namespace




