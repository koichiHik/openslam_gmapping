#include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include <gmapping/scanmatcher/scanmatcher.h>
#include "gridlinetraversal.h"
//#define GENERATE_MAPS

namespace GMapping {

using namespace std;

ScanMatcher::ScanMatcher() 
	: m_params(), m_laserBeams(0), 
		m_activeAreaComputed(false), 
		m_linePoints(new IntPoint[20000])
{}

ScanMatcher::~ScanMatcher(){
	delete [] m_linePoints;
}

void ScanMatcher::invalidateActiveArea(){
	m_activeAreaComputed=false;
}

void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (m_activeAreaComputed) {
		return;
	}
	OrientedPoint lp = p;
	lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
	lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
	lp.theta += m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	
	Point min(map.map2world(0,0));
	Point max(map.map2world(map.getMapSizeX() - 1, map.getMapSizeY() - 1));
	       
	min.x = std::min(lp.x, min.x);
	min.y = std::min(lp.y, min.y);
	max.x = std::max(lp.x, max.x);
	max.y = std::max(lp.y, max.y);
	
	/*determine the size of the area*/
	const double* angle = m_laserAngles + m_params.initialBeamSkip;
	for (const double* r = readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++){
		if (*r > m_params.laserMaxRange || *r == 0.0 || isnan(*r)) {
			continue;
		}
		double d = *r > m_params.usableRange ? m_params.usableRange : *r;
		Point phit = lp;
		phit.x += d * cos(lp.theta + *angle);
		phit.y += d * sin(lp.theta + *angle);
		if (phit.x<min.x) min.x=phit.x;

		min.x = std::min(min.x, phit.x);
		min.y = std::min(min.y, phit.y);
		max.x = std::max(max.x, phit.x);
		max.y = std::max(max.y, phit.y);
	}
	
	// Point is outside of map.
	if ( !map.isInside(min)	|| !map.isInside(max)){
		Point lmin(map.map2world(0,0));
		Point lmax(map.map2world(map.getMapSizeX() - 1, map.getMapSizeY() - 1));
		min.x = ( min.x >= lmin.x ) ? lmin.x : min.x - m_params.enlargeStep;
		max.x = ( max.x <= lmax.x ) ? lmax.x : max.x + m_params.enlargeStep;
		min.y = ( min.y >= lmin.y ) ? lmin.y : min.y - m_params.enlargeStep;
		max.y = ( max.y <= lmax.y ) ? lmax.y : max.y + m_params.enlargeStep;
		map.resize(min.x, min.y, max.x, max.y);
	}
	
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	// Allocate the active area.
	angle = m_laserAngles + m_params.initialBeamSkip;
	// Iterate for all sensor readings.
	for (const double* r = readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++) {

		// When map generation is requested.
		if (m_generateMap){
			double d = *r;
			if (d > m_params.laserMaxRange || d == 0.0 || isnan(d)) {
				continue;
			}
			d = std::min(d, m_params.usableRange);
			Point phit = lp + Point(d * cos(lp.theta + *angle), d * sin(lp.theta + *angle));
			IntPoint p0=map.world2map(lp);
			IntPoint p1=map.world2map(phit);
			
			GridLineTraversalLine line;
			line.points = m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i = 0; i < line.num_points - 1; i++){
				assert(map.isInside(m_linePoints[i]));
				activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
				assert(m_linePoints[i].x >= 0 && m_linePoints[i].y >= 0);
			}
			if (d < m_params.usableRange){
				IntPoint cp = map.storage().patchIndexes(p1);
				assert(cp.x >= 0 && cp.y >= 0);
				activeArea.insert(cp);
			}

		} else {
			if (*r > m_params.laserMaxRange || *r > m_params.usableRange || *r == 0.0 || isnan(*r)) {
				continue;
			}

			Point phit = lp;
			phit.x += *r * cos(lp.theta + *angle);
			phit.y += *r * sin(lp.theta + *angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x >= 0 && p1.y >= 0);
			IntPoint cp = map.storage().patchIndexes(p1);
			assert(cp.x >= 0 && cp.y >= 0);
			activeArea.insert(cp);
		}
	}
	
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;
}

double ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	
	if (!m_activeAreaComputed) {
		computeActiveArea(map, p, readings);
	}
		
	// This operation replicates the cells that will be changed in the registration operation.
	map.storage().allocActiveArea();
	
	// Transform from sensor to base_link.
	OrientedPoint lp = p;
	lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
	lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
	lp.theta += m_laserPose.theta;

	// Point of radiation.
	IntPoint p0 = map.world2map(lp);
	
	const double * angle = m_laserAngles + m_params.initialBeamSkip;
	double esum = 0;

	// Iterate through all beams.
	for (const double* r = readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++) {

		if (m_generateMap){
			double d = *r;
			// Reject all invalid readings.
			if (d > m_params.laserMaxRange || d == 0.0 || isnan(d)) {
				continue;
			}
			d = std::min(d, m_params.usableRange);
			Point phit = lp + Point(d * cos(lp.theta + *angle), d * sin(lp.theta + *angle));
			IntPoint p1 = map.world2map(phit);
			GridLineTraversalLine line;
			line.points = m_linePoints;

			// Collect all grid traversed by this line.
			// Kind of ray casting??
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i = 0; i < line.num_points - 1; i++){
				PointAccumulator& cell = map.cell(line.points[i]);
				double e =- cell.entropy();
				cell.update(false, Point(0,0));
				e += cell.entropy();
				esum += e;
			}
			if (d < m_params.usableRange){
				double e =- map.cell(p1).entropy();
				map.cell(p1).update(true, phit);
				e += map.cell(p1).entropy();
				esum += e;
			}

		} else {

			if (*r > m_params.laserMaxRange || *r > m_params.usableRange || *r == 0.0 || isnan(*r)) {
				continue;
			}

			Point phit = lp;
			phit.x += *r * cos(lp.theta + *angle);
			phit.y += *r * sin(lp.theta + *angle);
			IntPoint p1 = map.world2map(phit);
			assert(p1.x >= 0 && p1.y >= 0);
			map.cell(p1).update(true,phit);
		}
	}

	return esum;
}


double ScanMatcher::icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	double currentScore;
	double sc = score(map, init, readings);
	OrientedPoint start = init;
	pnew = init;
	int iterations = 0;
	do{
		currentScore = sc;
		sc = icpStep(pnew, map, start, readings);
		start = pnew;
		iterations++;
	} while (sc > currentScore);
	cerr << "i=" << iterations << endl;
	return currentScore;
}

double ScanMatcher::optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const {
	double bestScore = -1;
	OrientedPoint currentPose = init;
	double currentScore = score(map, currentPose, readings);
	double adelta = m_params.optAngularDelta, ldelta = m_params.optLinearDelta;
	unsigned int refinement=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	int c_iterations=0;

	do{
		if (bestScore >= currentScore){
			refinement++;
			adelta *= .5;
			ldelta *= .5;
		}
		bestScore = currentScore;

		OrientedPoint bestLocalPose = currentPose;
		OrientedPoint localPose = currentPose;

		Move move = Front;
		do {
			localPose = currentPose;
			switch(move){
				case Front:
					localPose.x += ldelta;
					move = Back;
					break;
				case Back:
					localPose.x -= ldelta;
					move=Left;
					break;
				case Left:
					localPose.y -= ldelta;
					move=Right;
					break;
				case Right:
					localPose.y += ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta += adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta -= adelta;
					move=Done;
					break;
				default:;
			}
			
			// Score weighting based on odometry reliability.
			double odo_gain = 1;
			if (m_params.angularOdometryReliability > 0.0){
				double dth = init.theta - localPose.theta;
				dth = atan2(sin(dth), cos(dth));
				dth *= dth;
				odo_gain *= exp(-m_params.angularOdometryReliability * dth);
			}
			if (m_params.linearOdometryReliability > 0.0){
				double dx = init.x - localPose.x;
				double dy = init.y - localPose.y;
				double drho = dx * dx + dy * dy;
				odo_gain *= exp(-m_params.linearOdometryReliability * drho);
			}
			double localScore = odo_gain * score(map, localPose, readings);
			
			// If local score is the highest.
			if (localScore>currentScore){
				currentScore = localScore;
				bestLocalPose = localPose;
			}
			c_iterations++;
		} while(move != Done);

		currentPose = bestLocalPose;

	}while (currentScore > bestScore || refinement < m_params.optRecursiveIterations);

	pnew = currentPose;
	return bestScore;
}

struct ScoredMove{
	OrientedPoint pose;
	double score;
	double likelihood;
};

typedef std::list<ScoredMove> ScoredMoveList;

double ScanMatcher::optimize(OrientedPoint& _mean, ScanMatcher::CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	ScoredMoveList moveList;
	double bestScore = -1;
	OrientedPoint currentPose = init;
	ScoredMove sm = {currentPose, 0, 0};
	unsigned int matched = likelihoodAndScore(sm.score, sm.likelihood, map, currentPose, readings);
	double currentScore = sm.score;
	moveList.push_back(sm);
	double adelta = m_params.optAngularDelta, ldelta = m_params.optLinearDelta;
	unsigned int refinement=0;
	int count=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	do{
		if (bestScore >= currentScore){
			refinement++;
			adelta *= 0.5;
			ldelta *= 0.5;
		}
		bestScore = currentScore;
		OrientedPoint bestLocalPose = currentPose;
		OrientedPoint localPose = currentPose;

		Move move = Front;
		do {
			localPose = currentPose;
			switch(move){
				case Front:
					localPose.x += ldelta;
					move = Back;
					break;
				case Back:
					localPose.x -= ldelta;
					move = Left;
					break;
				case Left:
					localPose.y -= ldelta;
					move = Right;
					break;
				case Right:
					localPose.y += ldelta;
					move = TurnLeft;
					break;
				case TurnLeft:
					localPose.theta += adelta;
					move = TurnRight;
					break;
				case TurnRight:
					localPose.theta -= adelta;
					move = Done;
					break;
				default:;
			}
			double localScore, localLikelihood;
			
			double odo_gain = 1;
			if (m_params.angularOdometryReliability > 0.0){
				double dth = init.theta - localPose.theta;
				dth = atan2(sin(dth), cos(dth));
				dth *= dth;
				odo_gain *= exp(-m_params.angularOdometryReliability * dth);
			}
			if (m_params.linearOdometryReliability > 0.0){
				double dx = init.x - localPose.x;
				double dy = init.y - localPose.y;
				double drho = dx * dx + dy * dy;
				odo_gain *= exp(-m_params.linearOdometryReliability * drho);
			}
			localScore = odo_gain * score(map, localPose, readings);
			//update the score
			count++;
			matched = likelihoodAndScore(localScore, localLikelihood, map, localPose, readings);
			if (localScore > currentScore){
				currentScore = localScore;
				bestLocalPose = localPose;
			}
			sm.score = localScore;
			sm.likelihood = localLikelihood;
			sm.pose = localPose;
			moveList.push_back(sm);
			//update the move list
		} while(move != Done);
		currentPose = bestLocalPose;
		//here we look for the best move;
	}while (currentScore > bestScore || refinement < m_params.optRecursiveIterations);
	
	//normalize the likelihood
	double lmin = 1e9;
	double lmax = -1e9;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		lmin = std::min(it->likelihood, lmin);
		lmax = std::max(it->likelihood, lmax);
	}

	for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++){
		it->likelihood = exp(it->likelihood - lmax);
	}
	//compute the mean
	OrientedPoint mean(0,0,0);
	double lacc = 0;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		mean = mean + it->pose * it->likelihood;
		lacc += it->likelihood;
	}
	mean = mean * (1./lacc);
	CovarianceMatrix cov = {0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		OrientedPoint delta = it->pose-mean;
		delta.theta = atan2(sin(delta.theta), cos(delta.theta));
		cov.xx += delta.x * delta.x * it->likelihood;
		cov.yy += delta.y * delta.y * it->likelihood;
		cov.tt += delta.theta * delta.theta *it->likelihood;
		cov.xy += delta.x * delta.y * it->likelihood;
		cov.xt += delta.x * delta.theta * it->likelihood;
		cov.yt += delta.y * delta.theta * it->likelihood;
	}
	cov.xx /= lacc, cov.xy /= lacc, cov.xt /= lacc, cov.yy /= lacc, cov.yt /= lacc, cov.tt /= lacc;
	
	_mean = currentPose;
	_cov = cov;
	return bestScore;
}

void ScanMatcher::setLaserParameters
	(unsigned int beams, double* angles, const OrientedPoint& lpose){

	assert(beams<LASER_MAXBEAMS);
	m_laserPose = lpose;
	m_laserBeams = beams;
	memcpy(m_laserAngles, angles, sizeof(double)*m_laserBeams);

}

double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings){

	ScoredMoveList moveList;
	
	// Create list of movement with scores. 
	for (double xx = -m_params.llSamplerange; xx <= m_params.llSamplerange; xx += m_params.llSampleStep) {
		for (double yy = -m_params.llSamplerange; yy <= m_params.llSamplerange; yy += m_params.llSampleStep) {
			for (double tt = -m_params.laSampleRange; tt <= m_params.laSampleRange; tt += m_params.laSampleStep) {
				
				OrientedPoint rp=p;
				rp.x += xx;
				rp.y += yy;
				rp.theta += tt;
				
				ScoredMove sm;
				sm.pose = rp;
				
				likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
				moveList.push_back(sm);
			}
		}
	}
	
	//normalize the likelihood
	double lmax = -1e9;
	double lcum = 0;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		lmax = std::max(it->likelihood, lmax);
	}
	for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++){
		lcum += exp(it->likelihood - lmax);
		it->likelihood = exp(it->likelihood - lmax);
	}
	
	OrientedPoint mean(0,0,0);
	double s = 0,c = 0;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		mean = mean + it->pose * it->likelihood;
		s += it->likelihood * sin(it->pose.theta);
		c += it->likelihood * cos(it->pose.theta);
	}

	mean = mean * (1./lcum);
	s /= lcum;
	c /= lcum;
	mean.theta = atan2(s,c);
	
	CovarianceMatrix cov = {0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++) {
		OrientedPoint delta = it->pose-mean;
		delta.theta = atan2(sin(delta.theta), cos(delta.theta));
		cov.xx += delta.x * delta.x * it->likelihood;
		cov.yy += delta.y * delta.y * it->likelihood;
		cov.tt += delta.theta * delta.theta * it->likelihood;
		cov.xy += delta.x * delta.y * it->likelihood;
		cov.xt += delta.x * delta.theta * it->likelihood;
		cov.yt += delta.y * delta.theta * it->likelihood;
	}
	cov.xx /= lcum, cov.xy /= lcum, cov.xt /= lcum, cov.yy /= lcum, cov.yt /= lcum, cov.tt /= lcum;
	
	_mean = mean;
	_cov = cov;
	_lmax = lmax;
	return log(lcum) + lmax;
}

double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p,
	Gaussian3& odometry, const double* readings, double gain){

	ScoredMoveList moveList;
	for (double xx = -m_params.llSamplerange; xx <= m_params.llSamplerange; xx += m_params.llSampleStep) {
		for (double yy = -m_params.llSamplerange; yy <= m_params.llSamplerange; yy += m_params.llSampleStep) {
			for (double tt = -m_params.laSampleRange; tt <= m_params.laSampleRange; tt += m_params.laSampleStep) {	
				OrientedPoint rp = p;
				rp.x += xx;
				rp.y += yy;
				rp.theta += tt;
				
				ScoredMove sm;
				sm.pose = rp;
				
				likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
				sm.likelihood += odometry.eval(rp)/gain;
				assert(!isnan(sm.likelihood));
				moveList.push_back(sm);
			}
		}
	}
	
	//OrientedPoint delta=mean-currentPose;
	//normalize the likelihood
  double lmax = -std::numeric_limits<double>::max();
	double lcum = 0;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		lmax = std::max(it->likelihood, lmax);
	}
	for (ScoredMoveList::iterator it = moveList.begin(); it != moveList.end(); it++){
		lcum += exp(it->likelihood - lmax);
		it->likelihood = exp(it->likelihood - lmax);
	}
	
	OrientedPoint mean(0,0,0);
	double s = 0, c = 0;
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		mean = mean + it->pose * it->likelihood;
		s += it->likelihood * sin(it->pose.theta);
		c += it->likelihood * cos(it->pose.theta);
	}
	mean = mean * (1./lcum);
	s /= lcum;
	c /= lcum;
	mean.theta = atan2(s,c);
	
	// Create covariance matrix.
	CovarianceMatrix cov = {0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it = moveList.begin(); it != moveList.end(); it++){
		OrientedPoint delta = it->pose - mean;
		delta.theta = atan2(sin(delta.theta), cos(delta.theta));
		cov.xx += delta.x * delta.x * it->likelihood;
		cov.yy += delta.y * delta.y * it->likelihood;
		cov.tt += delta.theta * delta.theta * it->likelihood;
		cov.xy += delta.x * delta.y * it->likelihood;
		cov.xt += delta.x * delta.theta * it->likelihood;
		cov.yt += delta.y * delta.theta * it->likelihood;
	}
	cov.xx /= lcum, cov.xy /= lcum, cov.xt /= lcum, cov.yy /= lcum, cov.yt /= lcum, cov.tt /= lcum;
	
	_mean = mean;
	_cov = cov;
	_lmax = lmax;
	double v = log(lcum) + lmax;
	assert(!isnan(v));
	return v;
}

void ScanMatcher::setMatchingParameters(const ScanMatcher::Params &params) {
	m_params = params;
}

};

