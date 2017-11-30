#include "HillTypeMuscleConstraint.h"
#include <iostream>
using namespace FEM;
void
HillTypeMuscleConstraint::
ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
{
	double f_muscle 		   = ComputeMuscleForce(mCacheL);
	double f_muscle_derivative = ComputeMuscleForceDerivative(mCacheL);

	Eigen::Matrix2d Fddt = mCacheF*mCacheddT;
	double dl = (((1.0/mCacheL)*Fddt).cwiseProduct(dF)).sum();
	dP = (-f_muscle/(mCacheL*mCacheL)+f_muscle_derivative/mCacheL)*dl*Fddt + 
		(f_muscle/mCacheL)*dF*mCacheddT;

}

double
HillTypeMuscleConstraint::
ComputeMuscleForceIntegrate(const double& l)
{
	double f_a_dl = mStiffness*fabs(l-1.0);
	double f_p_dl = 0.0*fabs(l-1.0);

	return mActivationLevel*f_a_dl + f_p_dl;
}
double
HillTypeMuscleConstraint::
ComputeMuscleForce(const double& l)
{
	//f = a*f_a + f_p
	double f_a = mStiffness;
	double f_p = 0.0;
	return mActivationLevel*f_a + f_p;
}
double
HillTypeMuscleConstraint::
ComputeMuscleForceDerivative(const double& l)
{
	double d_f_a = 0.0;
	double d_f_p = 0.0;

	return mActivationLevel*d_f_a + d_f_p;
}
HillTypeMuscleConstraint::
HillTypeMuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:MuscleConstraint(stiffness,fiber_direction,i0,i1,i2,vol,invDm)
{
}

double
HillTypeMuscleConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	double f_muscle = ComputeMuscleForce(mCacheL);

	Eigen::Vector2d x0(x.block<2,1>(mi0*2,0));

	auto x_01 = x.block<2,1>(mi1*2,0)-x0;
	auto x_02 = x.block<2,1>(mi2*2,0)-x0;

	Eigen::Matrix2d P = (f_muscle/mCacheL)*mCacheF*mCacheddT;

	
	P = mVol*P*(mInvDm.transpose());
	return ComputeMuscleForceIntegrate(mCacheL);
	std::cout<<mCacheL<<"\t\t"<<ComputeMuscleForceIntegrate(mCacheL)<<std::endl;
	std::cout<<(P.block<2,1>(0,0)).dot(x_01) + (P.block<2,1>(0,1)).dot(x_02)<<std::endl<<std::endl;
	// return (P.block<2,1>(0,0)).dot(x_01) + (P.block<2,1>(0,1)).dot(x_02);
}
void
HillTypeMuscleConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	ComputeDeformationGradient(x);
	double f_muscle = ComputeMuscleForce(mCacheL);
	Eigen::Matrix2d P = (f_muscle/mCacheL)*mCacheF*mCacheddT;

	P = mVol*P*(mInvDm.transpose());

	gradient.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
	gradient.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
	gradient.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
}
void
HillTypeMuscleConstraint::
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
HillTypeMuscleConstraint::
Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b)
{
	ComputeDeformationGradient(x);
	double save_a = mActivationLevel;
	mActivationLevel = 1.0;
	double f_muscle = ComputeMuscleForce(mCacheL);
	mActivationLevel = 0.0;
	f_muscle -=ComputeMuscleForce(mCacheL);
	mActivationLevel = save_a;

	Eigen::Matrix2d P = (f_muscle/mCacheL)*mCacheF*mCacheddT;

	P = mVol*P*(mInvDm.transpose());

	b.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
	b.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
	b.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
}



void
HillTypeMuscleConstraint::
EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	std::cout<<"HillTypeMuscleConstraint not supported"<<std::endl;
}
void
HillTypeMuscleConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	std::cout<<"HillTypeMuscleConstraint not supported"<<std::endl;
}
void
HillTypeMuscleConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	std::cout<<"HillTypeMuscleConstraint not supported"<<std::endl;
}

ConstraintType
HillTypeMuscleConstraint::
GetType()
{
	return ConstraintType::HILL_TYPE_MUSCLE;
}
