#ifndef __MUSCLE_CONSTRAINT_H__
#define __MUSCLE_CONSTRAINT_H__
#include "FEMConstraint.h"
namespace FEM
{
class MuscleConstraint : public FEMConstraint
{
protected:
	Eigen::Vector2d	mFiberDirection;
	Eigen::Matrix2d mCacheddT;
	double mCacheL;
	double mActivationLevel;
	
	virtual bool ComputeDeformationGradient(const Eigen::VectorXd& x);
public:
	MuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);

    virtual void Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b)=0;

	void SetActivationLevel(const double& a);
	const double& GetActivationLevel();
	const Eigen::Vector2d& GetFiberDirection();
};
};
#endif