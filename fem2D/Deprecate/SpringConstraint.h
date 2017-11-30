#ifndef __SPRING_CONSTRAINT_H__
#define __SPRING_CONSTRAINT_H__
#include "Constraint.h"
namespace FEM
{
class SpringConstraint : public Constraint
{
protected:
	int mi0,mi1;
	double ml0;
public:
	SpringConstraint(const double& stiffness,int i0,int i1,double l0);
	double  EvalPotentialEnergy(const Eigen::VectorXd& x) override;
	void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) override;
	void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) override;

	void GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) override;
    void EvaluateDVector(const Eigen::VectorXd& x) override;
    void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
    void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;

    int GetDof() override;
	int GetNumHessianTriplets() override;
	ConstraintType GetType() override;
	void AddOffset(const int& offset) override;

	const double& GetI0();
	const double& GetI1();
};
}
#endif