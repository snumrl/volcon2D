#include "MuscleConstraint.h"

using namespace FEM;
bool
MuscleConstraint::
ComputeDeformationGradient(const Eigen::VectorXd& x)
{
	FEMConstraint::ComputeDeformationGradient(x);
	mCacheL = (mCacheF*mFiberDirection).norm();
	return true;
}
MuscleConstraint::
MuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:FEMConstraint(stiffness,0.0,i0,i1,i2,vol,invDm),mFiberDirection(fiber_direction),mActivationLevel(0.0)
{
	mCacheddT = mFiberDirection*mFiberDirection.transpose();
}

void 
MuscleConstraint::
SetActivationLevel(const double& a) 
{
	mActivationLevel = a;
}
const double& 
MuscleConstraint::
GetActivationLevel() 
{
	return mActivationLevel;
}
const Eigen::Vector2d&
MuscleConstraint::
GetFiberDirection()
{
	return mFiberDirection;
}