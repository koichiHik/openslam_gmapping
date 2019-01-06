#ifndef SCANMATCHERPROCESSOR_H
#define SCANMATCHERPROCESSOR_H

//#include <gmapping/log/sensorlog.h>
#include <gmapping/sensor/sensor_range/rangesensor.h>
#include <gmapping/sensor/sensor_range/rangereading.h>
//#include <gsl/gsl_eigen.h>
#include <gmapping/scanmatcher/scanmatcher.h>

namespace GMapping {

class ScanMatcherProcessor{

	struct ScanMatcherProcParams {
		ScanMatcherProcParams() 
			: useICP(false), computeCovariance(false), 
			numBeams(0), regScore(300), critScore(regScore * 0.5),
			maxMove(1.0)
		{}

		bool useICP;
		bool computeCovariance;
		int numBeams;
		double regScore;
		double critScore;
		double maxMove;
	};

	struct Params {
		ScanMatcherProcParams scan_matcher_proc_param;
		ScanMatcher::Params scan_matcher_param;
	};

	public:
  ScanMatcherProcessor(const ScanMatcherMap& m);
  ScanMatcherProcessor (double xmin, double ymin, double xmax, double ymax, double delta, double patchdelta);
		virtual ~ScanMatcherProcessor ();
		virtual void processScan(const RangeReading & reading);
		void setSensorMap(const SensorMap& smap, std::string sensorName="FLASER");
		void init();
		void setScanMatcherProcParams(const ScanMatcherProcessor::Params& params);
		OrientedPoint getPose() const;
		inline const ScanMatcherMap& getMap() const {return m_map;}
		inline ScanMatcher& matcher() {return m_matcher;}

	protected:
		ScanMatcher m_matcher;
		bool m_first;
		SensorMap m_sensorMap;
		unsigned int m_beams;
		//state
		ScanMatcherMap m_map;
		OrientedPoint m_pose;
		OrientedPoint m_odoPose;
		int  m_count;
		//gsl_eigen_symmv_workspace * m_eigenspace;
		Params m_params;
};

};

#endif


