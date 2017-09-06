#ifndef __FEM_WORLD_H__
#define __FEM_WORLD_H__
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <vector>

namespace FEM
{
enum IntegrationMethod
{
	SYMPLECTIC_EXPLICIT,
	SEMI_IMPLICIT,
	IMPLICIT_NEWTON_METHOD,
	QUASI_STATIC,
	PROJECTIVE_DYNAMICS,
	PROJECTIVE_QUASI_STATIC
};



class Constraint;
class AttachmentConstraint;
class MuscleConstraint;

class World
{
public:
	World(
		IntegrationMethod integration_method = SYMPLECTIC_EXPLICIT,
		double time_step = 1.0/120.0,
		int max_iteration = 100,
		double damping_coeff = 0.999
		);
	void 								Initialize();
	void 								TimeStepping(bool isIntegrated = true);
	Eigen::VectorXd						ComputeJacobian(std::vector<MuscleConstraint*>& mc,std::vector<AttachmentConstraint*>& attachment_vector);

	void 								AddBody(const Eigen::VectorXd& x0,const std::vector<Constraint*>& c,const double& mass = 1.0);
	void 								AddConstraint(Constraint* c);
	void 								RemoveConstraint(Constraint* c);

	int 								GetClosestNode(const Eigen::Vector2d& x);
	double 								GetTimeStep();
	double 								GetTime();
	int    								GetNumVertices();
	const Eigen::VectorXd& 				GetPositions();
	void 								SetPositions(const Eigen::VectorXd& p);
	std::vector<Constraint*>&			GetConstraints();

private:	
	Eigen::VectorXd 					IntegrateSymplecticExplicit();
	Eigen::VectorXd 					IntegrateSemiImplicit();
	Eigen::VectorXd 					IntegrateNewtonMethod();
	Eigen::VectorXd 					IntegrateQuasiStatic();
	Eigen::VectorXd 					IntegrateProjectiveDynamics();
	Eigen::VectorXd 					IntegrateProjectiveQuasiStatic();
	
	void 								FactorizeLLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& lltSolver);
	void 								FactorizeLDLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldltSolver);
	Eigen::VectorXd 					ConjugateGradient(const Eigen::VectorXd& b,const Eigen::VectorXd& x0);
	double 								ComputeStepSize(const Eigen::VectorXd& x,const Eigen::VectorXd& gradient,const Eigen::VectorXd& descent);
	void 								PreComputeProjectiveDynamics();
	void 								InversionFree(Eigen::VectorXd& x);

	double 								EvalObjectives(const Eigen::VectorXd& x);
	void 								EvalGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g);
	void 								EvalHessian(const Eigen::VectorXd& x,const Eigen::VectorXd& dx,Eigen::VectorXd& dg);

	double 								EvalConstraintsObjectives(const Eigen::VectorXd& x);
	void 								EvalConstraintsGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g);
	void 								EvalConstraintsHessian(const Eigen::VectorXd& x,const Eigen::VectorXd& dx,Eigen::VectorXd& dg);

	void 								EvaluateDVector(const Eigen::VectorXd& x,Eigen::VectorXd& d);
    void 								EvaluateJMatrix(Eigen::SparseMatrix<double>& J);
    void 								EvaluateLMatrix(Eigen::SparseMatrix<double>& L);

	void 								ComputeExternalForces();
	void 								UpdatePositionsAndVelocities(const Eigen::VectorXd& x_n1);
private:
	bool 								mIsInitialized;
	int 								mNumVertices;
	int 								mConstraintDofs;
	int 								mMaxIteration;
	int 								mFrame;

	double 								mTimeStep,mTime;
	double 								mDampingCoefficinent;

	IntegrationMethod 					mIntegrationMethod;
	
	std::vector<double>			 		mUnitMass;
	std::vector<Constraint*>			mConstraints;

	Eigen::VectorXd						mX,mV;
	Eigen::VectorXd						mExternalForces;

	Eigen::SparseMatrix<double> 		mMassMatrix;
	Eigen::SparseMatrix<double> 		mInvMassMatrix;
	Eigen::SparseMatrix<double> 		mIdentityMatrix;

	Eigen::VectorXd						mQn;
	Eigen::SparseMatrix<double> 		mJ,mL;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> mDynamicSolver;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> mQuasiStaticSolver;
};

};
#endif
