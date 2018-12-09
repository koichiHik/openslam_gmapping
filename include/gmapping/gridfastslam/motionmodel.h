#ifndef MOTIONMODEL_H
#define MOTIONMODEL_H

#include <gmapping/utils/point.h>
#include <gmapping/utils/stat.h>
#include <gmapping/utils/macro_params.h>

namespace  GMapping { 

struct MotionModel{

	struct Params {
		double err_ratio_rho_per_rho;
		double err_ratio_theta_per_rho;
		double err_ratio_rho_per_theta;
		double err_ratio_theta_per_theta;
	};

	/** Update pose of the particle. **/
	OrientedPoint drawFromMotion(const OrientedPoint& p, double linearMove, double angularMove) const;
	OrientedPoint drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const;
	Covariance3 gaussianApproximation(const OrientedPoint& pnew, const OrientedPoint& pold) const;
	void setMotionModelParams(const MotionModel::Params& params);
	Params m_params;
};

};

#endif
