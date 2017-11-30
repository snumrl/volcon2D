#include "BallInfo.h"

using namespace dart::dynamics;
using namespace dart::simulation;

BallInfo::
BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel)
	:isReleased(true),constraint(cons),skeleton(skel)
{

}

void
BallInfo::
ComputeFallingPosition(double h,Eigen::Vector3d& fp)
{
	Eigen::Vector3d p = skeleton->getBodyNode(0)->getCOM();
	Eigen::Vector3d v = skeleton->getBodyNode(0)->getCOMLinearVelocity();

	double dx = h-p[1];
	double g = -9.8;
	double v2_plus_2gdx = v[1]*v[1] + 2.0*g*dx;
	if(v2_plus_2gdx<0)
	{
		std::cout<<"no solution"<<std::endl;
		return;
	}

	double t1 = (-v[1] - sqrt(v2_plus_2gdx))/g;
	double t2 = (-v[1] + sqrt(v2_plus_2gdx))/g;
	double t = (t1>t2?t1:t2);
	if(t<0.0)
		t=0.05;

	fp = p+t*v;
	if(fp[1]>h)
	fp[1] = h;
}
void
BallInfo::
GetPosition(Eigen::Vector3d& p)
{
	p =  skeleton->getBodyNode(0)->getCOM();
}
void
BallInfo::
GetVelocity(Eigen::Vector3d& v)
{
	v = skeleton->getBodyNode(0)->getCOMLinearVelocity();
}
void
BallInfo::
Release(const WorldPtr& world)
{
	if(!isReleased){
		world->getConstraintSolver()->removeConstraint(constraint);
		isReleased = true;
	}
}

void
BallInfo::
Attach(const WorldPtr& world,BodyNode* bn)
{
	if(isReleased)
	{
		constraint.reset();
		constraint = std::make_shared<dart::constraint::WeldJointConstraint>(skeleton->getBodyNode(0),bn);
		isReleased = false;

		world->getConstraintSolver()->addConstraint(constraint);
	}
	
	
}
