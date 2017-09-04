// #ifndef __NEO_HOOKEAN_FEM_CONSTRAINT_H__
// #define __NEO_HOOKEAN_FEM_CONSTRAINT_H__
// #include "FEMConstraint.h"

// class NeoHookeanFEMConstraint : public FEMConstraint
// {
// protected:
// 	Eigen::Matrix2d	mCacheFTF,mCacheInvFT,mCacheInvF;

// 	virtual void ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP);
// 	virtual bool ComputeDeformationGradient(const Eigen::VectorXd& x);
	
// public:
// 	NeoHookeanFEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);

// 	virtual double  EvalPotentialEnergy(const Eigen::VectorXd& x);
// 	virtual void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient);
// 	virtual void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg);

// 	virtual void EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d);
//     virtual void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets);
//     virtual void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets);

//     virtual int GetDof(){return 2;};
// 	virtual int GetNumHessianTriplets();
// 	virtual ConstraintType GetType()	   {return ConstraintType::NEO_HOOKEAN; };
// };

// #endif
