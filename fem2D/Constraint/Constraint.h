#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__	
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Core>
namespace FEM
{
enum ConstraintType
{
	ATTACHMENT,
	SPRING,
	COROTATE,
	NEO_HOOKEAN,
	VENANT_KIRCHHOFF,
	HILL_TYPE_MUSCLE,
	LINEAR_MUSCLE,
	PW_LINEAR_MUSCLE,
	TENDON
};
class Constraint
{
protected:
	double mStiffness;
public:
	Constraint(const double& stiffness);
	virtual double  EvalPotentialEnergy(const Eigen::VectorXd& x) = 0;
    virtual void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient) = 0;
    virtual void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg) = 0;

    virtual void GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d) = 0;
    virtual void EvaluateDVector(const Eigen::VectorXd& x) = 0;
    virtual void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets) = 0;
    virtual void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets) = 0;

	virtual int GetNumHessianTriplets() = 0;
	virtual int GetDof() = 0;
	virtual ConstraintType GetType() = 0;
	virtual void AddOffset(const int& offset) = 0;
};
};
#endif