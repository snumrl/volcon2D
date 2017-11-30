#include "SpringConstraint.h"
using namespace FEM;
SpringConstraint::
SpringConstraint(const double& stiffness,int i0,int i1,double l0)
	:Constraint(stiffness),mi0(i0),mi1(i1),ml0(l0)
{

}
double
SpringConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	Eigen::Vector2d x_01 = x.block<2,1>(mi1*2,0)-x.block<2,1>(mi0*2,0);

	return 0.5*mStiffness*x_01.squaredNorm();
}
void
SpringConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	Eigen::Vector2d x_01 = x.block<2,1>(mi1*2,0)-x.block<2,1>(mi0*2,0);	
	double l = x_01.norm();
	gradient.block<2,1>(mi0*2,0) += -mStiffness*(1.0-ml0/l)*x_01;
	gradient.block<2,1>(mi1*2,0) += mStiffness*(1.0-ml0/l)*x_01;
}
// void
// SpringConstraint::
// EvalHessian(const Eigen::VectorXd& x, std::vector<Eigen::Triplet<double>>& hessian_triplets)
// {
// 	Eigen::Vector2d x_01 = x.block<2,1>(mi1*2,0)-x.block<2,1>(mi0*2,0);
// 	double l = x_01.norm();
// 	Eigen::Matrix2d K = mStiffness*(1.0-ml0/l)*Eigen::Matrix2d::Identity() + mStiffness*ml0/(l*l*l)*x_01*x_01.transpose();

// 	for(int i=0;i<2;i++)
// 	{
// 		for(int j=0;j<2;j++)
// 		{
// 			double val = K(i,j);

// 			hessian_triplets.push_back(Eigen::Triplet<double>(mi0*2+i,mi0*2+j,val));
// 			hessian_triplets.push_back(Eigen::Triplet<double>(mi1*2+i,mi0*2+j,-val));
// 			hessian_triplets.push_back(Eigen::Triplet<double>(mi0*2+i,mi1*2+j,-val));
// 			hessian_triplets.push_back(Eigen::Triplet<double>(mi1*2+i,mi1*2+j,val));
// 		}
// 	}

// }

void
SpringConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
	Eigen::Vector2d x_01 = x.block<2,1>(mi1*2,0)-x.block<2,1>(mi0*2,0);	
	Eigen::Vector2d dx_01 = dx.block<2,1>(mi1*2,0)-dx.block<2,1>(mi0*2,0);	
	double l = x_01.norm();
	Eigen::Vector2d Kx = mStiffness*((1.0-ml0/l)*dx_01 + ml0/(l*l*l)*(x_01.dot(dx_01))*x_01);
	
	dg.block<2,1>(mi0*2,0) += -Kx;
	dg.block<2,1>(mi1*2,0) += Kx;
}
void
SpringConstraint::
EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	Eigen::Vector2d x_01 = x.block<2,1>(mi1*2,0)-x.block<2,1>(mi0*2,0);	

	d.block<2,1>(index*2,0) = -ml0*x_01.normalized();
}
void
SpringConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
    J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+0, mStiffness));
    J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+1, mStiffness));
    J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+0, -mStiffness));
    J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+1, -mStiffness));
}
void
SpringConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*mi0+0, mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*mi0+1, mStiffness));

	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*mi0+0, -mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*mi0+1, -mStiffness));

	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*mi1+0, -mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*mi1+1, -mStiffness));

	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*mi1+0, mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*mi1+1, mStiffness));
}



int
SpringConstraint::
GetDof()
{
	return 1;
}
int
SpringConstraint::
GetNumHessianTriplets()
{
	return 16;
}
ConstraintType
SpringConstraint::
GetType()	   
{
	return ConstraintType::SPRING; 
}
void
SpringConstraint::
AddOffset(const int& offset) 
{
	mi0 +=offset;
	mi1 +=offset;
}

const double&
SpringConstraint::
GetI0()
{
	return mi0;
}

const double&
SpringConstraint::
GetI1()
{
	return mi1;
}

