// #include "VenantKirchhoffFEMConstraint.h"
// #include <iostream>

// void
// VenantKirchhoffFEMConstraint::
// ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
// {

// }
// bool
// VenantKirchhoffFEMConstraint::
// ComputeDeformationGradient(const Eigen::VectorXd& x)
// {
// 	if(FEMConstraint::ComputeDeformationGradient(x))
// 	{
// 		mCacheE = 0.5*((mCacheF.transpose()*mCacheF)-Eigen::Matrix2d::Identity());
// 		return true;
// 	}
// 	return false;
// }
// VenantKirchhoffFEMConstraint::
// VenantKirchhoffFEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
// 	:FEMConstraint(stiffness,poisson_ratio,i0,i1,i2,vol,invDm),
// 	mCacheE(Eigen::Matrix2d::Zero())
// {
// }
// double
// VenantKirchhoffFEMConstraint::
// EvalPotentialEnergy(const Eigen::VectorXd& x)
// {
// 	ComputeDeformationGradient(x);

// 	double vol_preserve_sqrt = (mCacheE).trace();
// 	return mVol*(0.5*mMu*(mCacheE.squaredNorm())+0.5*mLambda*vol_preserve_sqrt*vol_preserve_sqrt);
// }
// void
// VenantKirchhoffFEMConstraint::
// EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
// {
// 	ComputeDeformationGradient(x);

// 	Eigen::Matrix2d P = mVol*(mCacheF*(mMu*mCacheE + (mLambda*(mCacheE.trace()))*Eigen::Matrix2d::Identity()));
	
// 	gradient.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
// 	gradient.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
// 	gradient.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
// }

// void
// VenantKirchhoffFEMConstraint::
// EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
// {

// }

// void
// VenantKirchhoffFEMConstraint::
// EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
// {
// 	std::cout<<"VenantKirchhoffFEMConstraint not supported."<<std::endl;
// }
// void
// VenantKirchhoffFEMConstraint::
// EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
// {
// 	std::cout<<"VenantKirchhoffFEMConstraint not supported."<<std::endl;
// }
// void
// VenantKirchhoffFEMConstraint::
// EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
// {
// 	std::cout<<"VenantKirchhoffFEMConstraint not supported."<<std::endl;
// }
// int
// VenantKirchhoffFEMConstraint::
// GetNumHessianTriplets()
// {
// 	return 36;
// }
