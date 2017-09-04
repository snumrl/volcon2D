#include "FEMConstraint.h"
#include <iostream>
#include <Eigen/LU>
using namespace FEM;
bool
FEMConstraint::
ComputeDeformationGradient(const Eigen::VectorXd& x)
{
	Eigen::Vector2d x0(x.block<2,1>(mi0*2,0));

	Eigen::Matrix2d Ds;

	Ds.block<2,1>(0,0) = x.block<2,1>(mi1*2,0)-x0;
	Ds.block<2,1>(0,1) = x.block<2,1>(mi2*2,0)-x0;

	mCacheDs = Ds;
	mCacheF = mCacheDs * mInvDm;

	return true;
}

FEMConstraint::
FEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:Constraint(stiffness),
	mMu(stiffness/((1.0+poisson_ratio))),mLambda(stiffness*poisson_ratio/((1.0+poisson_ratio)*(1-2.0*poisson_ratio))),
	mi0(i0),mi1(i1),mi2(i2),mVol(vol),mInvDm(invDm),
	mCacheDs(Eigen::Matrix2d::Zero()),
	mCacheF(Eigen::Matrix2d::Zero())
{
}

void
FEMConstraint::
AddInversionFreePosition(Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	double vol = mCacheF.determinant();
	if(vol<0)
	{
		Eigen::Vector2d proj_0;
		Eigen::Vector2d proj_1;
		Eigen::Vector2d proj_2;

		Eigen::Vector2d p0,p1,p2;
		p0 = x.block<2,1>(2*mi0,0);
		p1 = x.block<2,1>(2*mi1,0);
		p2 = x.block<2,1>(2*mi2,0);

		proj_0 = (((p2-p1).dot(p0-p1))/((p2-p1).squaredNorm()))*(p2-p1) + p1;
		proj_1 = (((p2-p0).dot(p1-p0))/((p2-p0).squaredNorm()))*(p2-p0) + p0;
		proj_2 = (((p1-p0).dot(p2-p0))/((p1-p0).squaredNorm()))*(p1-p0) + p0;
		
		double len_0 = (p0-proj_0).squaredNorm();		
		double len_1 = (p1-proj_1).squaredNorm();		
		double len_2 = (p2-proj_2).squaredNorm();
		double eps = 1E-7;

		if(len_0<len_1&&len_0<len_2)
			x.block<2,1>(mi0*2,0) = proj_0+eps*(proj_0-p0);
		else if(len_1<len_0&&len_1<len_2)
			x.block<2,1>(mi1*2,0) = proj_1+eps*(proj_1-p1);
		else
			x.block<2,1>(mi2*2,0) = proj_2+eps*(proj_2-p2);
	}
}
int
FEMConstraint::
GetDof()
{
	return 2;
}
int
FEMConstraint::
GetNumHessianTriplets()
{
	return 36;
}
void
FEMConstraint::
AddOffset(const int& offset) 
{
	mi0+=offset;
	mi1+=offset;
	mi2+=offset;
}

const int&
FEMConstraint::
GetI0()
{
	return mi0;
}
const int&
FEMConstraint::
GetI1()
{
	return mi1;
}
const int&
FEMConstraint::
GetI2()
{
	return mi2;
}