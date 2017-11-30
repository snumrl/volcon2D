#ifndef __LINEAR_MUSCLE_CONSTRAINT_H__
#define __LINEAR_MUSCLE_CONSTRAINT_H__
#include "MuscleConstraint.h"
#include <iostream>
namespace FEM
{
class LinearMuscleConstraint : public MuscleConstraint
{
protected:
	Eigen::Vector2d			md;
	virtual void ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP);
public:
	LinearMuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);

	double  EvalPotentialEnergy(const Eigen::VectorXd& x) override;
	void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) override;
	void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) override;

	void GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) override;
    void EvaluateDVector(const Eigen::VectorXd& x) override;
    void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
    void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;

	void Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b) override;
    int GetDof() override;
	int GetNumHessianTriplets() override;
	ConstraintType GetType() override;
};
};
#endif