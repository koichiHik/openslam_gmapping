#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#include "icp.h"
#include "smmap.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>
#define LASER_MAXBEAMS 2048

namespace GMapping {

class ScanMatcher{
	public:

		struct Params {
			Params() 
				: //laserPose(0, 0, 0), 
				optRecursiveIterations(3),
				llSamplerange(0.01),
				llSampleStep(0.01),
				laSampleRange(0.005),
				laSampleStep(0.005),
				enlargeStep(10),
				fullnessThreshold(0.1),
				angularOdometryReliability(0.0),
				linearOdometryReliability(0.0),
				freeCellRatio(sqrt(2)),
				initialBeamSkip(0),
				nullLikelihood(-0.5)
			{}

			//OrientedPoint laserPose;
			double laserMaxRange, usableRange;
			int kernelSize;
			double optLinearDelta;
			double optAngularDelta;
			double gaussianSigma, likelihoodSigma;
			int optRecursiveIterations;	
			unsigned int likelihoodSkip;

			double llSamplerange, llSampleStep, laSampleRange, laSampleStep;
			double enlargeStep;
			double fullnessThreshold;
			double angularOdometryReliability, linearOdometryReliability;
			double freeCellRatio;
			unsigned int initialBeamSkip;
			double nullLikelihood;
		};

		typedef Covariance3 CovarianceMatrix;
		
		ScanMatcher();
		~ScanMatcher();
		
		double icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double optimize(OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		// This function register scan. Meaning associate reading to the particle of location p.
		double registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);

		void setLaserParameters(unsigned int beams, double* angles, const OrientedPoint& lpose);
		
		void setMatchingParameters (const ScanMatcher::Params& params);

		inline void setGenerateMap(bool val) { m_generateMap = val; }

		inline bool getGenerateMap() const {return m_generateMap; }

		void invalidateActiveArea();
	
		void computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);

		inline double icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		inline double score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		inline unsigned int likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double likelihood(double& lmax, OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		
		double likelihood(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, Gaussian3& odometry, const double* readings, double gain=180.);


		inline const double* laserAngles() const { return m_laserAngles; }
		
		inline unsigned int laserBeams() const { return m_laserBeams; }
			
	protected:
		//state of the matcher
		bool m_activeAreaComputed;
		
		/**laser parameters*/
		OrientedPoint m_laserPose;
		unsigned int m_laserBeams;
		double       m_laserAngles[LASER_MAXBEAMS];

		Params m_params;

	public:
	
		// allocate this large array only once
		IntPoint* m_linePoints;

	private:
		bool m_generateMap;

};

inline double ScanMatcher::icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	const double* angle = m_laserAngles + m_params.initialBeamSkip;
	
	// Convert from sensor to baselink.
	OrientedPoint lp = p;
	lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
	lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
	lp.theta += m_laserPose.theta;

	unsigned int skip = 0;
	double freeDelta = map.getDelta() * m_params.freeCellRatio;
	std::list<PointPair> pairs;
	
	// Iterate all beams.
	for (const double* r=readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++){
		skip++;
		skip = skip > m_params.likelihoodSkip ? 0 : skip;
		if (*r > m_params.usableRange || *r == 0.0) {
			continue;
		}
		if (skip) {
			continue;
		}

		// Occupied Labeling
		Point phit = lp;
		phit.x += *r * cos(lp.theta + *angle);
		phit.y += *r * sin(lp.theta + *angle);
		IntPoint iphit = map.world2map(phit);

		// Free Labeling
		Point pfree = lp;
		pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle);
		pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);
 		pfree = pfree - phit;
		IntPoint ipfree = map.world2map(pfree);

		bool found = false;
		Point bestMu(0.,0.);
		Point bestCell(0.,0.);

		// Matching within kernel.
		for (int xx = -m_params.kernelSize; xx <= m_params.kernelSize; xx++) {
			for (int yy = -m_params.kernelSize; yy <= m_params.kernelSize; yy++) {

				IntPoint pr = iphit + IntPoint(xx,yy);
				IntPoint pf = pr + ipfree;
				const PointAccumulator& cell = map.cell(pr);
				const PointAccumulator& fcell = map.cell(pf);

				if (((double)cell )> m_params.fullnessThreshold && ((double)fcell ) < m_params.fullnessThreshold){

					Point mu = phit - cell.mean();
					if (!found){
						bestMu = mu;
						bestCell = cell.mean();
						found = true;
					}else {
						if((mu * mu) < (bestMu * bestMu)){
							bestMu = mu;
							bestCell = cell.mean();
						} 
					}

				}
			}
		}

		if (found){
			pairs.push_back(std::make_pair(phit, bestCell));
		}
	}
	
	OrientedPoint result(0,0,0);
	//double icpError=icpNonlinearStep(result,pairs);
	std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta << std::endl;
	pret.x = p.x + result.x;
	pret.y = p.y + result.y;
	pret.theta = p.theta + result.theta;
	pret.theta = atan2(sin(pret.theta), cos(pret.theta));

	return score(map, p, readings);
}

inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	double s = 0;
	const double* angle = m_laserAngles + m_params.initialBeamSkip;
	OrientedPoint lp = p;
	lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
	lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
	lp.theta += m_laserPose.theta;
	unsigned int skip = 0;
	double freeDelta = map.getDelta() * m_params.freeCellRatio;

	// Iterate for all scans.
	for (const double* r = readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++) {
		skip++;
		skip = skip > m_params.likelihoodSkip ? 0 : skip;
		if (skip || *r > m_params.usableRange || *r == 0.0) {
			continue;
		}
		
		// Occupied Labling
		Point phit = lp;
		phit.x += *r * cos(lp.theta + *angle);
		phit.y += *r * sin(lp.theta + *angle);
		IntPoint iphit = map.world2map(phit);

		// Free Labeling.
		Point pfree = lp;
		pfree.x += (*r - map.getDelta() * freeDelta) * cos(lp.theta + *angle);
		pfree.y += (*r - map.getDelta() * freeDelta) * sin(lp.theta + *angle);
 		pfree = pfree - phit;
		IntPoint ipfree = map.world2map(pfree);

		bool found = false;
		Point bestMu(0.,0.);
		for (int xx = -m_params.kernelSize; xx <= m_params.kernelSize; xx++) {
			for (int yy = -m_params.kernelSize; yy <= m_params.kernelSize; yy++) {
				IntPoint pr=iphit+IntPoint(xx,yy);
				IntPoint pf=pr+ipfree;
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )> m_params.fullnessThreshold && ((double)fcell ) < m_params.fullnessThreshold){
					Point mu = phit - cell.mean();
					if (!found){
						bestMu = mu;
						found = true;
					} else {
						bestMu = (mu*mu) < (bestMu*bestMu) ? mu:bestMu;
					}
				}
			}
		}

		if (found) {
			s += exp(-1.0 / m_params.gaussianSigma * bestMu * bestMu);
		}
	}
	return s;
}

inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	using namespace std;
	l = 0;
	s = 0;
	const double * angle = m_laserAngles + m_params.initialBeamSkip;

	// Reflection Point.
	OrientedPoint lp=p;
	lp.x += cos(p.theta) * m_laserPose.x - sin(p.theta) * m_laserPose.y;
	lp.y += sin(p.theta) * m_laserPose.x + cos(p.theta) * m_laserPose.y;
	lp.theta += m_laserPose.theta;

	unsigned int skip = 0;
	unsigned int c = 0;
	double noHit = m_params.nullLikelihood / (m_params.likelihoodSigma);
	double freeDelta = map.getDelta() * m_params.freeCellRatio;

	for (const double* r = readings + m_params.initialBeamSkip; r < readings + m_laserBeams; r++, angle++){
		skip++;
		skip = skip > m_params.likelihoodSkip ? 0 : skip;
		if (*r > m_params.usableRange)  {
			continue;
		}
		if (skip) {
			continue;
		}

		// Occupied Labeling
		Point phit = lp;
		phit.x += *r * cos(lp.theta + *angle);
		phit.y += *r * sin(lp.theta + *angle);
		IntPoint iphit = map.world2map(phit);

		// Free Labeling
		Point pfree = lp;
		pfree.x += (*r - freeDelta) * cos(lp.theta+*angle);
		pfree.y += (*r - freeDelta) * sin(lp.theta+*angle);
		pfree = pfree - phit;
		IntPoint ipfree = map.world2map(pfree);

		bool found = false;
		Point bestMu(0.,0.);
		for (int xx = -m_params.kernelSize; xx <= m_params.kernelSize; xx++) {
			for (int yy = -m_params.kernelSize; yy <= m_params.kernelSize; yy++){
				IntPoint pr = iphit + IntPoint(xx,yy);
				IntPoint pf = pr+ipfree;
				const PointAccumulator& cell = map.cell(pr);
				const PointAccumulator& fcell = map.cell(pf);
				if (((double)cell ) > m_params.fullnessThreshold && ((double)fcell ) < m_params.fullnessThreshold){
					Point mu = phit - cell.mean();
					if (!found){
						bestMu = mu;
						found = true;
					}else
						bestMu = (mu*mu) < (bestMu*bestMu) ? mu : bestMu;
				}
			}
		}

		if (found){
			s += exp(-1.0 / m_params.gaussianSigma * bestMu * bestMu);
			c++;
		}

		if (!skip){
			double f = (-1.0 / m_params.likelihoodSigma) * bestMu * bestMu;
			l += (found) ? f : noHit;
		}
	}
	return c;
}

};

#endif
