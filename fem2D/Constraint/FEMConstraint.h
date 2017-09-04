#ifndef __FEM_CONSTRAINT_H__
#define __FEM_CONSTRAINT_H__
#include "Constraint.h"

namespace FEM
{
class FEMConstraint : public Constraint
{
protected:

	int mi0,mi1,mi2;
	double mVol;
	double mMu,mLambda;
	Eigen::Matrix2d mInvDm,mCacheDs,mCacheF;

	virtual void ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP) = 0;
	virtual bool ComputeDeformationGradient(const Eigen::VectorXd& x);
public:
	FEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);
   
   	int GetDof() override;
   	int GetNumHessianTriplets() override;

    void AddInversionFreePosition(Eigen::VectorXd& x);
	void AddOffset(const int& offset);
	const int& GetI0();
	const int& GetI1();
	const int& GetI2();
};
};
#endif
