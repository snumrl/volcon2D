#include "LinearMuscleConstraint.h"
#include <iostream>
using namespace FEM;
LinearMuscleConstraint::
LinearMuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:MuscleConstraint(stiffness,fiber_direction,i0,i1,i2,vol,invDm)
{

}
void
LinearMuscleConstraint::
ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
{
	dP = mStiffness*(dF*mFiberDirection*(mFiberDirection.transpose()));
}

double 
LinearMuscleConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	auto p = (1.0-mActivationLevel)*mCacheF*mFiberDirection;
	
	return 0.5*mVol*mStiffness*((mCacheF*mFiberDirection - p).squaredNorm());
}
void
LinearMuscleConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	ComputeDeformationGradient(x);
	auto p = (1.0-mActivationLevel)*mCacheF*mFiberDirection;
	
	Eigen::Matrix2d P = mStiffness*(mCacheF*mFiberDirection*(mFiberDirection.transpose())-p*(mFiberDirection.transpose()));

	P = mVol*P*(mInvDm.transpose());

	gradient.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
	gradient.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
	gradient.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
}
void
LinearMuscleConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
	ComputeDeformationGradient(x);
	Eigen::Matrix2d dDs,dF,dP;
	Eigen::Vector2d dx0(dx.block<2,1>(mi0*2,0));
	dDs.block<2,1>(0,0) = dx.block<2,1>(mi1*2,0)-dx0;
	dDs.block<2,1>(0,1) = dx.block<2,1>(mi2*2,0)-dx0;
	
	dF = dDs*(mInvDm);
	ComputedP(dF,dP);

	dP = mVol * dP * (mInvDm.transpose());

	dg.block<2,1>(mi0*2,0) += -(dP.block<2,1>(0,0) + dP.block<2,1>(0,1));
	dg.block<2,1>(mi1*2,0) += dP.block<2,1>(0,0);
	dg.block<2,1>(mi2*2,0) += dP.block<2,1>(0,1);
}
void
LinearMuscleConstraint::
GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.block<2,1>(2*index,0) = md;
}
void
LinearMuscleConstraint::
EvaluateDVector(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	md = (1.0-mActivationLevel)*mCacheF*mFiberDirection;
}
void
LinearMuscleConstraint::
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
LinearMuscleConstraint::
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
LinearMuscleConstraint::
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
LinearMuscleConstraint::
GetDof()
{
	return 1;
}

int
LinearMuscleConstraint::
GetNumHessianTriplets() 
{
	return 36;
}

ConstraintType
LinearMuscleConstraint::
GetType()	   
{
	return ConstraintType::LINEAR_MUSCLE; 
}
