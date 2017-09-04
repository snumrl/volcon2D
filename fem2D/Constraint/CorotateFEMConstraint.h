#ifndef __COROTATE_FEM_CONSTRAINT_H__
#define	__COROTATE_FEM_CONSTRAINT_H__
#include "FEMConstraint.h"

namespace FEM
{
class CorotateFEMConstraint : public FEMConstraint
{
protected:
	Eigen::Matrix2d			mCacheU,mCacheV,mCacheR,mCacheD;
	
	void ComputeSVD(const Eigen::Matrix2d& F);
	void ComputedR(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dR);

	void ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP) override;
	bool ComputeDeformationGradient(const Eigen::VectorXd& x) override;
public:
	CorotateFEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);

	double  EvalPotentialEnergy(const Eigen::VectorXd& x) override;
	void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) override;
	void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) override;

	void EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) override;
    void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
    void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;

	ConstraintType GetType() override;
};

};
#endif