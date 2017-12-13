#include "HillTypeMuscleConstraint.h"
#include <iostream>
using namespace FEM;
HillTypeMuscleConstraint::
HillTypeMuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:MuscleConstraint(stiffness,fiber_direction,i0,i1,i2,vol,invDm),l_previous(1.0)
{

}

void
HillTypeMuscleConstraint::
ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
{
}

double 
HillTypeMuscleConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
}
void
HillTypeMuscleConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
}
void
HillTypeMuscleConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
}
double
HillTypeMuscleConstraint::
EvalLambda()
{

	return (mCacheF*mFiberDirection).norm();
}
double
HillTypeMuscleConstraint::
f_l_curve(double l)
{
	return exp(-2.0*(l-1)*(l-1));
}
double
HillTypeMuscleConstraint::
f_l_dot_curve(double l_dot)
{
	double fld;
	double l_dot_max = 10;
	double Af = 0.3;
	double f_len = 1.8;
	if(l_dot<=0)
	{
		fld = (l_dot + l_dot_max)/(l_dot_max - l_dot/Af);
	}
	else
	{
		fld = (f_len*l_dot*(2+2/Af) + l_dot_max*(f_len-1))/(l_dot*(2+2/Af)+l_dot_max*(f_len-1));
	}
	return fld;
}
double
HillTypeMuscleConstraint::
f_l_passive_curve(double l)
{
	// double kpe = 4.0;
	// double emo = 0.6;
	// if(l>0)
	// 	return (exp(kpe*(l-1)/emo)-1)/(exp(kpe)-1);
	// else
		return 0;
}
double
HillTypeMuscleConstraint::
EvalP0(double a)
{
	double l = EvalLambda();

	double l_dot = -cbrt(mVol)*(l-l_previous)*200;
	// l_dot=0.0;
	// std::cout<<l_dot<<std::endl;
	double b = 0.7;
	double f_hill = mActivationLevel*f_l_curve(l)*f_l_dot_curve(l_dot) + f_l_passive_curve(l)-b*l_dot;
	double p0 = 1- 1.0/l*(f_hill);
	return p0;
}
void
HillTypeMuscleConstraint::
GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.block<2,1>(2*index,0) = md;
}
void
HillTypeMuscleConstraint::
EvaluateDVector(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	md = EvalP0(mActivationLevel)*mCacheF*mFiberDirection;
}
void
HillTypeMuscleConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	Eigen::MatrixXd Ai(2,6);

	auto v = mInvDm*mFiberDirection;
	double a = v[0];
	double b = v[1];

	Ai<<
		-(a+b),0,a,0,b,0,
		0,-(a+b),0,a,0,b;

	auto MuAiT = mVol*mStiffness*Ai.transpose();

	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+0, MuAiT(2*0+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+1, MuAiT(2*0+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+0, MuAiT(2*0+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+1, MuAiT(2*0+1,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+0, MuAiT(2*1+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+1, MuAiT(2*1+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+0, MuAiT(2*1+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+1, MuAiT(2*1+1,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+0, MuAiT(2*2+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+1, MuAiT(2*2+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+0, MuAiT(2*2+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+1, MuAiT(2*2+1,2*0+1)));
};
void
HillTypeMuscleConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	Eigen::MatrixXd Ai(2,6);
	auto v = mInvDm*mFiberDirection;
	double a = v[0];
	double b = v[1];

	Ai<<
		-(a+b),0,a,0,b,0,
		0,-(a+b),0,a,0,b;

	auto MuAiTAi = mVol*mStiffness*((Ai.transpose())*Ai);

	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi0+0,MuAiTAi(2*0+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi0+1,MuAiTAi(2*0+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi0+0,MuAiTAi(2*0+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi0+1,MuAiTAi(2*0+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi1+0,MuAiTAi(2*0+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi1+1,MuAiTAi(2*0+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi1+0,MuAiTAi(2*0+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi1+1,MuAiTAi(2*0+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi2+0,MuAiTAi(2*0+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi2+1,MuAiTAi(2*0+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi2+0,MuAiTAi(2*0+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi2+1,MuAiTAi(2*0+1, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi0+0,MuAiTAi(2*1+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi0+1,MuAiTAi(2*1+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi0+0,MuAiTAi(2*1+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi0+1,MuAiTAi(2*1+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi1+0,MuAiTAi(2*1+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi1+1,MuAiTAi(2*1+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi1+0,MuAiTAi(2*1+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi1+1,MuAiTAi(2*1+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi2+0,MuAiTAi(2*1+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi2+1,MuAiTAi(2*1+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi2+0,MuAiTAi(2*1+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi2+1,MuAiTAi(2*1+1, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi0+0,MuAiTAi(2*2+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi0+1,MuAiTAi(2*2+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi0+0,MuAiTAi(2*2+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi0+1,MuAiTAi(2*2+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi1+0,MuAiTAi(2*2+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi1+1,MuAiTAi(2*2+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi1+0,MuAiTAi(2*2+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi1+1,MuAiTAi(2*2+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi2+0,MuAiTAi(2*2+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi2+1,MuAiTAi(2*2+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi2+0,MuAiTAi(2*2+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi2+1,MuAiTAi(2*2+1, 2*2+1)));
};

void
HillTypeMuscleConstraint::
Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b)
{
	ComputeDeformationGradient(x);

	Eigen::Matrix2d P = mStiffness*mCacheF*mCacheddT;

	P = mVol*P*(mInvDm.transpose());

	b.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
	b.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
	b.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
}

int
HillTypeMuscleConstraint::
GetDof()
{
	return 1;
}

int
HillTypeMuscleConstraint::
GetNumHessianTriplets() 
{
	return 36;
}

ConstraintType
HillTypeMuscleConstraint::
GetType()	   
{
	return ConstraintType::HILL_TYPE_MUSCLE; 
}
void
HillTypeMuscleConstraint::
SetPreviousL(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	if(mi2 == 1)
	{
		double l_dot = (EvalLambda()-l_previous)*10;
		// std::cout<<l_dot<<std::endl;
		// std::cout<<f_l_dot_curve(l_dot)<<std::endl;
	}
		// std::cout<<l_previous<<std::endl;
	l_previous = EvalLambda();
}