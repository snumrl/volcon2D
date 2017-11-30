#include "PWLinearMuscleConstraint.h"
#include <iostream>

using namespace FEM;
double
PWLinearMuscleConstraint::
ComputeL(const Eigen::VectorXd& x)
{
	double l = 0.0;

	for(int i =0;i<mIndex.size()-1;i++)
		l += (x.block<2,1>(mIndex[i]*2,0) - x.block<2,1>(mIndex[i+1]*2,0)).norm();
	return l;
}
double
PWLinearMuscleConstraint::
ComputeDesiredL()
{
	return (1.0 - 0.5*mActivationLevel)*ml0;
}
double
PWLinearMuscleConstraint::
ComputeU(const double& tau)
{
	double u =0.0;
	double e = 0.01;
	double k = mStiffness;
	if(tau<-e)
		u = 0.0;
	else if(tau<=e)
		u = k*(1.0/(6.0*e)*tau*tau*tau + 0.5*tau*tau + 0.5*e*tau + 1.0/6.0*e*e);
	else
		u = k*(tau*tau + e*e/3.0);
	
	return u;
}
double
PWLinearMuscleConstraint::
ComputedU(const double& tau)
{
	double du =0.0;
	double e = 0.01;
	double k = mStiffness;
	if(tau<-e)
		du = 0.0;
	else if(tau<=e)
		du = k*(0.5/e*tau*tau + tau + 0.5*e);
	else
		du = k*(2.0*tau);

	return du;	
}
double
PWLinearMuscleConstraint::
ComputeddU(const double& tau)
{
	double ddu =0.0;
	double e = 0.01;
	double k = mStiffness;
	if(tau<-e)
		ddu = 0.0;
	else if(tau<=e)
		ddu = k*(1.0/e*tau + 1.0);
	else
		ddu = k*(2.0);

	return ddu;
}

PWLinearMuscleConstraint::
PWLinearMuscleConstraint(const double& stiffness,const std::vector<int>& index,const double& l0)
	:Constraint(stiffness),mIndex(index),mActivationLevel(0.0),ml0(l0)
{
}
double
PWLinearMuscleConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	double tau = ComputeL(x) - ComputeDesiredL();
	tau -= 0.01;
	return ComputeU(tau);
}
void
PWLinearMuscleConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	double tau = ComputeL(x) - ComputeDesiredL();
	tau -= 0.01;
	double du = ComputedU(tau);
	for(int i =0;i<mIndex.size()-1;i++)
	{
		Eigen::Vector2d xi_xi1 = x.block<2,1>(mIndex[i]*2,0) - x.block<2,1>(mIndex[i+1]*2,0);
		double l = (xi_xi1).norm();
		gradient.block<2,1>(mIndex[i]*2,0) += du*(xi_xi1)/l;
		gradient.block<2,1>(mIndex[i+1]*2,0) += -du*(xi_xi1)/l;
	}
}
void
PWLinearMuscleConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
	double tau = ComputeL(x) - ComputeDesiredL();
	tau -= 0.01;
	double du = ComputedU(tau);
	double ddu = ComputeddU(tau);
	for(int i =0;i<mIndex.size()-1;i++)
	{
		Eigen::Vector2d xi_xi1 = x.block<2,1>(mIndex[i]*2,0) - x.block<2,1>(mIndex[i+1]*2,0);
		Eigen::Matrix2d xxT = xi_xi1*xi_xi1.transpose();
		double l = (xi_xi1).norm();
		double invl = 1.0/l;
		Eigen::Matrix2d H_ii = ddu*invl*xxT + du*(invl*Eigen::Matrix2d::Identity() - (invl*invl*invl)*xxT);

		Eigen::Vector2d dx_ii1 = dx.block<2,1>(mIndex[i]*2,0) - dx.block<2,1>(mIndex[i+1]*2,0);
		dg.block<2,1>(mIndex[i]*2,0) += H_ii*dx_ii1;
		dg.block<2,1>(mIndex[i+1]*2,0) += -H_ii*dx_ii1;
	}
}
void
PWLinearMuscleConstraint::
Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b)
{
	double save = mActivationLevel;
	mActivationLevel +=0.03;
	double tau = ComputeL(x) - ComputeDesiredL();
	tau -= 0.01;
	double ddu_da = ComputedU(tau);
	mActivationLevel = save;
	tau = ComputeL(x) - ComputeDesiredL();
	tau -= 0.01;
	ddu_da -= ComputedU(tau);
	ddu_da /= 0.03;
	
	for(int i =0;i<mIndex.size()-1;i++)
	{
		Eigen::Vector2d xi_xi1 = x.block<2,1>(mIndex[i]*2,0) - x.block<2,1>(mIndex[i+1]*2,0);
		double l = (xi_xi1).norm();
		b.block<2,1>(mIndex[i]*2,0) += ddu_da*(xi_xi1)/l;
		b.block<2,1>(mIndex[i+1]*2,0) += -ddu_da*(xi_xi1)/l;
	}
}
void
PWLinearMuscleConstraint::
EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	std::cout<<"PWLinearMuscleConstraint not supported"<<std::endl;
}
void
PWLinearMuscleConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	std::cout<<"PWLinearMuscleConstraint not supported"<<std::endl;
}
void
PWLinearMuscleConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	std::cout<<"PWLinearMuscleConstraint not supported"<<std::endl;
}

int
PWLinearMuscleConstraint::
GetNumHessianTriplets()
{
	return 16;
}
int
PWLinearMuscleConstraint::
GetDof()
{
	return 1;
}
ConstraintType
PWLinearMuscleConstraint::
GetType()
{
	return ConstraintType::PW_LINEAR_MUSCLE;
}
void
PWLinearMuscleConstraint::
AddOffset(const int& offset)
{
	for(auto& i : mIndex)
		i+=offset;
}
void 
PWLinearMuscleConstraint::
SetActivationLevel(const double& a)
{
	if(a<0)
		mActivationLevel=0;
	else if(a>1.0)
		mActivationLevel=1.0;
	else
		mActivationLevel = a;
}

const double&
PWLinearMuscleConstraint::
GetActivationLevel()
{
	return mActivationLevel;
}

const std::vector<int>&
PWLinearMuscleConstraint::
GetIndex()
{
	return mIndex;
}