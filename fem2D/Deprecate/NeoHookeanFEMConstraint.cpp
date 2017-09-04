// #include "NeoHookeanFEMConstraint.h"
// #include <Eigen/LU>
// #include <math.h>
// #include <iostream>
// bool
// NeoHookeanFEMConstraint::
// ComputeDeformationGradient(const Eigen::VectorXd& x)
// {

// 	if(FEMConstraint::ComputeDeformationGradient(x))
// 	{
// 		mCacheFTF = (mCacheF.transpose())*mCacheF;
// 		mCacheInvF = (mCacheF.inverse());
// 		mCacheInvFT = mCacheInvF.transpose();
// 		// if(fabs(mCacheF.determinant())<1E-6)
// 			// std::cout<<mCacheF<<std::endl;
// 		return true;
// 	}
// 	return false;
// }
// void
// NeoHookeanFEMConstraint::
// ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
// {
// 	double I3 = mCacheFTF.determinant();

// 	dP = 
// 		mMu*dF+
// 		(mMu - mLambda*log(I3)*0.5)*mCacheInvFT*(dF.transpose())*mCacheInvFT+
// 		(mLambda*((mCacheInvF*dF).trace()))*mCacheInvFT;
// }
// NeoHookeanFEMConstraint::
// NeoHookeanFEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
// 	:FEMConstraint(stiffness,poisson_ratio,i0,i1,i2,vol,invDm),
// 	mCacheFTF(Eigen::Matrix2d::Zero()),
// 	mCacheInvFT(Eigen::Matrix2d::Zero())
// {
// }
// double
// NeoHookeanFEMConstraint::
// EvalPotentialEnergy(const Eigen::VectorXd& x)
// {
// 	ComputeDeformationGradient(x);

// 	double I1 = mCacheFTF.trace();
// 	double I3 = mCacheFTF.determinant();

// 	return mVol*(0.25*mMu*(I1-log(I3)-3) + 0.125*mLambda*(log(I3)*log(I3)));
// }
// void
// NeoHookeanFEMConstraint::
// EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
// {
// 	ComputeDeformationGradient(x);

// 	double I3 = mCacheFTF.determinant();
// 	// std::cout<<I3<<std::endl;
// 	Eigen::Matrix2d P = (0.5*mMu)*(mCacheF-mCacheInvFT) + (0.5*mLambda*log(I3))*mCacheInvFT;

// 	P = mVol*P*mInvDm;

// 	gradient.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
// 	gradient.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
// 	gradient.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
// }


// void
// NeoHookeanFEMConstraint::
// EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
// {
// 	ComputeDeformationGradient(x);
// 	Eigen::Matrix2d dDs,dF,dP;
// 	Eigen::Vector2d dx0(dx.block<2,1>(mi0*2,0));
// 	dDs.block<2,1>(0,0) = dx.block<2,1>(mi1*2,0)-dx0;
// 	dDs.block<2,1>(0,1) = dx.block<2,1>(mi2*2,0)-dx0;
	
// 	dF = dDs*(mInvDm);
// 	ComputedP(dF,dP);

// 	dP = mVol * dP * (mInvDm.transpose());

// 	dg.block<2,1>(mi0*2,0) += -(dP.block<2,1>(0,0) + dP.block<2,1>(0,1));
// 	dg.block<2,1>(mi1*2,0) += dP.block<2,1>(0,0);
// 	dg.block<2,1>(mi2*2,0) += dP.block<2,1>(0,1);
// }

// void
// NeoHookeanFEMConstraint::
// EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
// {
// 	std::cout<<"NeoHookeanFEMConstraint not supported."<<std::endl;
// }
// void
// NeoHookeanFEMConstraint::
// EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
// {
// 	std::cout<<"NeoHookeanFEMConstraint not supported."<<std::endl;
// }
// void
// NeoHookeanFEMConstraint::
// EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
// {
// 	std::cout<<"NeoHookeanFEMConstraint not supported."<<std::endl;
// }
// int
// NeoHookeanFEMConstraint::
// GetNumHessianTriplets()
// {
// 	return 36;
// }
