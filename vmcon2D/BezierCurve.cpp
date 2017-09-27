#include "BezierCurve.h"

BezierCurve::
BezierCurve()
	:p0(Eigen::Vector2d::Zero()),p1(Eigen::Vector2d::Zero()),p2(Eigen::Vector2d::Zero()),invT(1.0)
{

}
BezierCurve::
BezierCurve(const Eigen::Vector2d& _p0,const Eigen::Vector2d& _p1,const Eigen::Vector2d& _p2,const double& _T)
	:p0(_p0),p1(_p1),p2(_p2),invT(1.0/_T)
{
}
void
BezierCurve::
Initialize(const Eigen::Vector2d& p0,const Eigen::Vector2d& p1,const Eigen::Vector2d& p2,const double& t)
{
	this->p0 = p0;
	this->p1 = p1;
	this->p2 = p2;
	this->invT = 1.0/t;
}
Eigen::Vector2d
BezierCurve::
GetPosition(double t)
{
	double s = invT*t;
	return 	p0*(1-s)*(1-s)+
			p1*2*s*(1-s)+
			p2*s*s;
}

Eigen::Vector2d
BezierCurve::
GetVelocity(double t)
{
	double s = invT*t;

	return invT*(
			2*(1.0-s)*p0 +
			(4.0*s+2.0)*p1 +
			2*s*p2);
}