#include <iostream>
#include "scanmatcherprocessor.h"
#include "eig3.h"

//#define SCANMATHCERPROCESSOR_DEBUG
namespace GMapping {

using namespace std;

ScanMatcherProcessor::ScanMatcherProcessor(const ScanMatcherMap& m)
  : m_map(m.getCenter(), m.getWorldSizeX(), m.getWorldSizeY(), m.getResolution()), 
    m_pose(0,0,0) {}


ScanMatcherProcessor::ScanMatcherProcessor 
  (double xmin, double ymin, double xmax, double ymax, double delta, double patchdelta):
  m_map(Point((xmax+xmin)*.5, (ymax+ymin)*.5), xmax-xmin, ymax-ymin, patchdelta), 
	m_pose(0,0,0) {}

ScanMatcherProcessor::~ScanMatcherProcessor () {}


void ScanMatcherProcessor::setSensorMap(const SensorMap& smap, std::string sensorName) {

	m_sensorMap = smap;
	
	/*
	 Construct the angle table for the sensor
	 FIXME has to be extended to more than one laser... 
	*/
	
	SensorMap::const_iterator laser_it = m_sensorMap.find(sensorName);
	assert(laser_it! = m_sensorMap.end());
	const RangeSensor* rangeSensor = dynamic_cast<const RangeSensor*>((laser_it->second));
	assert(rangeSensor && rangeSensor->beams().size());
	
	m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
	double* angles = new double[rangeSensor->beams().size()];
	for (unsigned int i=0; i < m_beams; i++){
		angles[i]=rangeSensor->beams()[i].pose.theta;
	}
	m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
	delete [] angles;
	
}

void ScanMatcherProcessor::init(){
	m_first = true;
	m_pose = OrientedPoint(0,0,0);
	m_count = 0;
}

void ScanMatcherProcessor::processScan(const RangeReading & reading){

	// Retireve the position from the reading, and compute the odometry.
	OrientedPoint relPose = reading.getPose();
	if (!m_count){
		m_odoPose = relPose;
	}

	// Compute the move in the scan matcher. 
	OrientedPoint move = relPose - m_odoPose;
	double dth = m_odoPose.theta - m_pose.theta;

	double lin_move = move*move;
	if (lin_move > m_params.scan_matcher_proc_param.maxMove){
		cerr << "Too big jump in the log file: " << lin_move << endl;
		cerr << "relPose=" << relPose.x << " " <<relPose.y << endl;
		cerr << "ignoring" << endl;
		return;
	}
	
	double s = sin(dth), c = cos(dth);
	OrientedPoint dPose;
	dPose.x = c * move.x - s * move.y;
	dPose.y = s * move.x + c * move.y;
	dPose.theta = move.theta;

	// Update pose by adding difference.
	m_pose = m_pose + dPose;
	m_pose.theta = atan2(sin(m_pose.theta), cos(m_pose.theta));
	m_odoPose = relPose;
	
	assert(reading.size()==m_beams);

	double * plainReading = new double[m_beams];
	reading.rawView(plainReading, m_map.getDelta());
	
	//the final stuff: scan match the pose
	double score=0;
	OrientedPoint newPose=m_pose;
	if (m_count){
		if(m_params.scan_matcher_proc_param.computeCovariance){
			ScanMatcher::CovarianceMatrix cov;
			score=m_matcher.optimize(newPose, cov, m_map, m_pose, plainReading);

			double m[3][3];
			double evec[3][3];
			double eval[3];
			m[0][0] = cov.xx; 
			m[0][1] = cov.xy; 
			m[0][2] = cov.xt;
			m[1][0] = cov.xy;
			m[1][1] = cov.yy;
			m[1][2] = cov.yt;
			m[2][0] = cov.xt;
			m[2][1] = cov.yt;
			m[2][2] = cov.tt;
			eigen_decomposition(m,evec,eval);

		} else {

			if (m_params.scan_matcher_proc_param.useICP){
				cerr << "USING ICP" << endl;
				score = m_matcher.icpOptimize(newPose, m_map, m_pose, plainReading);
			}else {
				score = m_matcher.optimize(newPose, m_map, m_pose, plainReading);
			}
		}
	}

	// Register Scan.
	if (!m_count || score < m_params.scan_matcher_proc_param.regScore){
		m_matcher.invalidateActiveArea();

		// If the score is below the threshold, use odom base pose.
		if (score < m_params.scan_matcher_proc_param.critScore){
			m_matcher.registerScan(m_map, m_pose, plainReading);
		} else {
			m_matcher.registerScan(m_map, newPose, plainReading);
		}	
	}

	m_pose = newPose;
	delete [] plainReading;
	m_count++;
}

void ScanMatcherProcessor::setScanMatcherProcParams(const ScanMatcherProcessor::Params& params) {
	m_params = params;
	m_matcher.setMatchingParameters(params.scan_matcher_param);
}

OrientedPoint ScanMatcherProcessor::getPose() const{
	return m_pose;
}

};

