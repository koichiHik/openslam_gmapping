

#include <string>
#include <fstream>
#include <iomanip>
#include <gmapping/gridfastslam/gridslamprocessor.h>


namespace GMapping {

  using namespace std;

  void GridSlamProcessor::printOdomJumpWarnings(const GMapping::OrientedPoint & new_pose, double dist_thresh) {

    cerr << "***********************************************************************" << endl;
    cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
    cerr << "m_distanceThresholdCheck=" << dist_thresh << endl;
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

}