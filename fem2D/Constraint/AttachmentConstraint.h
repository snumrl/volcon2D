#ifndef __ATTACHMENT_CONSTRAINT_H__
#define __ATTACHMENT_CONSTRAINT_H__	
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace FEM
{
class Constraint;
enum ConstraintType;

class AttachmentConstraint : public Constraint
{
protected:
	int mi0;
	Eigen::Vector2d mp;
	Eigen::Vector2d md;
public:
	AttachmentConstraint(const double& stiffness,int i0,const Eigen::Vector2d& p);

	double  EvalPotentialEnergy(const Eigen::VectorXd& x) override;
	void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) override;
	void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) override;

	void GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) override;
    void EvaluateDVector(const Eigen::VectorXd& x) override;
    void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) override;
    void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) override;

	int GetNumHessianTriplets() override;
	int GetDof() override;

	ConstraintType GetType() override;
	void AddOffset(const int& offset) override;

	Eigen::Vector2d& GetP();
	int&			 GetI0();
};

};
#endif
