#include "Constraint.h"
#include "AttachmentConstraint.h"
#include <Eigen/Dense>
#include <iostream>
using namespace FEM;
AttachmentConstraint::
AttachmentConstraint(const double& stiffness,int i0,const Eigen::Vector2d& p)
	:Constraint(stiffness),mi0(i0),mp(p)
{

}

double
AttachmentConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	Eigen::Vector2d x_p0 = x.block<2,1>(mi0*2,0) - mp;

    return 0.5*mStiffness*x_p0.squaredNorm();
}
void
AttachmentConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	Eigen::Vector2d x_p0 = x.block<2,1>(mi0*2,0) - mp;
	gradient.block<2,1>(mi0*2,0) += mStiffness*x_p0;
}

void
AttachmentConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
	//Compute H*x
	dg.block<2,1>(mi0*2,0) += mStiffness*dx.block<2,1>(mi0*2,0);
}

void
AttachmentConstraint::
EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.block<2,1>(2*index,0) = mp;
}
void
AttachmentConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0, 2*index, mStiffness));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+1, mStiffness));
}
void
AttachmentConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*mi0+0, mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*mi0+1, mStiffness));
}

int
AttachmentConstraint::
GetDof()
{
	return 1;
}
int
AttachmentConstraint::
GetNumHessianTriplets()
{
	return 2;
}


ConstraintType 
AttachmentConstraint::
GetType()	   
{
	return ConstraintType::ATTACHMENT; 
}
void 
AttachmentConstraint::
AddOffset(const int& offset) 
{
	mi0+=offset;
}

Eigen::Vector2d& 
AttachmentConstraint::
GetP() 
{
	return mp;
}
int&			 
AttachmentConstraint::
GetI0()	   
{
	return mi0;
}