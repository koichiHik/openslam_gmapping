#include <gmapping/gridfastslam/motionmodel.h>
#include <gmapping/utils/stat.h>
#include <iostream>

#define MotionModelConditioningLinearCovariance 0.01
#define MotionModelConditioningAngularCovariance 0.001

namespace GMapping {

OrientedPoint 
MotionModel::drawFromMotion (const OrientedPoint& p, double linearMove, double angularMove) const{
	OrientedPoint n(p);
	
	double lm=linearMove  + fabs( linearMove ) * sampleGaussian( m_params.err_ratio_rho_per_rho ) +
		fabs( angularMove ) * sampleGaussian( m_params.err_ratio_theta_per_rho );
	
	double am=angularMove + fabs( linearMove ) * sampleGaussian( m_params.err_ratio_rho_per_theta ) +
		fabs( angularMove ) * sampleGaussian( m_params.err_ratio_theta_per_theta );
	
	n.x += lm*cos(n.theta + .5 * am);
	n.y += lm*sin(n.theta + .5 * am);
	n.theta += am;
	n.theta = atan2(sin(n.theta), cos(n.theta));

	return n;
}

OrientedPoint 
MotionModel::drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const{
	double sxy = 0.3 * m_params.err_ratio_rho_per_rho;
	OrientedPoint delta = absoluteDifference(pnew, pold);
	OrientedPoint noisypoint(delta);

	noisypoint.x += sampleGaussian(m_params.err_ratio_rho_per_rho * fabs(delta.x) + 
		m_params.err_ratio_theta_per_rho * fabs(delta.theta) + sxy * fabs(delta.y));
	
	noisypoint.y += sampleGaussian(m_params.err_ratio_rho_per_rho * fabs(delta.y) +
		m_params.err_ratio_theta_per_rho * fabs(delta.theta) + sxy * fabs(delta.x));
	
	noisypoint.theta += sampleGaussian(m_params.err_ratio_theta_per_theta * fabs(delta.theta) + 
		m_params.err_ratio_rho_per_theta * sqrt(delta.x * delta.x + delta.y * delta.y));
	
	noisypoint.theta = fmod(noisypoint.theta, 2*M_PI);

	// Normalization.
	if (noisypoint.theta > M_PI) {
		noisypoint.theta -= 2 * M_PI;
	}	

	return absoluteSum(p,noisypoint);
}

Covariance3 MotionModel::gaussianApproximation(const OrientedPoint& pnew, const OrientedPoint& pold) const{
	OrientedPoint delta = absoluteDifference(pnew,pold);
	double linearMove = sqrt(delta.x * delta.x + delta.y * delta.y);
	double angularMove = fabs(delta.x);

	double s11 = m_params.err_ratio_rho_per_rho * m_params.err_ratio_rho_per_rho * linearMove * linearMove;
	double s22 = m_params.err_ratio_theta_per_theta * m_params.err_ratio_theta_per_theta * angularMove * angularMove;
	double s12 = m_params.err_ratio_theta_per_rho * m_params.err_ratio_theta_per_rho * angularMove * linearMove;
	
	Covariance3 cov;
	double s = sin(pold.theta), c = cos(pold.theta);
	cov.xx = c * c * s11 + MotionModelConditioningLinearCovariance;
	cov.yy = s * s * s11 + MotionModelConditioningLinearCovariance;
	cov.tt = s22 + MotionModelConditioningAngularCovariance;
	cov.xy = s * c * s11;
	cov.xt = c * s12;
	cov.yt = s * s12;
	return cov;
}

void MotionModel::setMotionModelParams(const MotionModel::Params& params) {
	this->m_params = params;
}

};

