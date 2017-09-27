#ifndef __BEZIER_CURVE_H__
#define __BEZIER_CURVE_H__
#include <Eigen/Core>
#include <vector>

class BezierCurve
{
	Eigen::Vector2d p0,p1,p2;
	double			invT;
public:
	BezierCurve();
	BezierCurve(const Eigen::Vector2d& p0,const Eigen::Vector2d& p1,const Eigen::Vector2d& p2,const double& t=1.0);
	void Initialize(const Eigen::Vector2d& p0,const Eigen::Vector2d& p1,const Eigen::Vector2d& p2,const double& t=1.0);
	double GetT() {return 1.0/invT;};
	Eigen::Vector2d GetPosition(double t);
	Eigen::Vector2d GetVelocity(double t);
};


#endif