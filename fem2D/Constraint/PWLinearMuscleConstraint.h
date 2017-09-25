#ifndef __PW_LINEAR_MUSCLE_CONSTRAINT_H__
#define __PW_LINEAR_MUSCLE_CONSTRAINT_H__
#include "Constraint.h"


namespace FEM
{
class PWLinearMuscleConstraint : public Constraint
{
	std::vector<int> mIndex;
	double mActivationLevel;
	double ml0;

	double ComputeL(const Eigen::VectorXd& x);
	double ComputeDesiredL();
	double ComputeU(const double& tau);
	double ComputedU(const double& tau);
	double ComputeddU(const double& tau);
public:
	PWLinearMuscleConstraint(const double& stiffness,const std::vector<int>& index,const double& l0);
	double  EvalPotentialEnergy(const Eigen::VectorXd& x) override;
    void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) override;
    void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) override;

    void EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) override;
    void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
    void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;

    void Evaldg_da(const Eigen::VectorXd& x,Eigen::VectorXd& b);
	int GetNumHessianTriplets() override;
	int GetDof() override;
	ConstraintType GetType() override;
	void AddOffset(const int& offset) override;


	void SetActivationLevel(const double& a);
	const double& GetActivationLevel();
	
	const std::vector<int>& GetIndex();	
};
};

#endif