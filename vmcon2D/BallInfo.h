#ifndef __BALL_INFO_H__
#define __BALL_INFO_H__

#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
class BallInfo
{
public:

	dart::dynamics::SkeletonPtr 				skeleton;
	dart::constraint::WeldJointConstraintPtr 	constraint;

	bool 										isReleased;
public:
	BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel);

	void ComputeFallingPosition(double h,Eigen::Vector3d& fp);
	void Release(const dart::simulation::WorldPtr& world);
	void GetPosition(Eigen::Vector3d& p);
	void GetVelocity(Eigen::Vector3d& v);
	void Attach(const dart::simulation::WorldPtr& world,dart::dynamics::BodyNode* bn);
	void TimeStepping();
};

#endif