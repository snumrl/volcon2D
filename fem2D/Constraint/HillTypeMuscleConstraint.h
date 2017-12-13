#ifndef __HILL_TYPE_MUSCLE_CONSTRAINT_H__
#define __HILL_TYPE_MUSCLE_CONSTRAINT_H__
#include "MuscleConstraint.h"
#include <iostream>
namespace FEM
{
class HillTypeMuscleConstraint : public MuscleConstraint
{
protected:
	Eigen::Vector2d			md;
	double l_previous;
	virtual void ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP);

	double EvalLambda();
	double f_l_curve(double l);
	double f_l_dot_curve(double l_dot);
	double f_l_passive_curve(double l);
	double EvalP0(double a);
public:
	HillTypeMuscleConstraint(const double& stiffness,const Eigen::Vector2d& fiber_direction,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm);

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

	void SetPreviousL(const Eigen::VectorXd& x);
};
};
#endif