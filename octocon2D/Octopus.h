#ifndef __OCTOPUS_H__
#define __OCTOPUS_H__
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include <vector>

struct Muscle
{
	double											activationLevel;
	std::vector<FEM::LinearMuscleConstraint*>		muscleConstraints;	
};

struct Target {
	int idx;
	Eigen::Vector2d coord;
};

class Octopus
{
private:
	FEM::Mesh*											mMesh;
	std::vector<Muscle*>								mMuscles;
	Eigen::VectorXd										mActivationLevel;
	std::vector<FEM::Constraint*>						mConstraints;
	std::vector<Target> 								mTarget;

	double								mMuscleStiffness;
	double								mYoungsModulus;
	double								mPoissonRatio;

	std::vector<FEM::AttachmentConstraint*>	mAttachementConstraintVector;
public:
	Octopus();

	/* Solve IK */
	void AddTarget(const Target& target);
	std::vector<Target> GetTarget();	
	void SolveSoftIK(FEM::World* world);

	int GetNumMuscles() { return mMuscles.size();};
	std::vector<Muscle*>& GetMuscles(){return mMuscles;};
	void AddMuscle(
		const std::vector<Eigen::Vector3i> indexList,
		const Eigen::Vector2d& fiber_direction
		);

	void Initialize(FEM::World* world);

	void TransformAttachmentPoints();
	void SetActivationLevel(const Eigen::VectorXd& a);
	void SoftSolveIK(FEM::World* world);
	Eigen::VectorXd& GetActivationLevel(){return mActivationLevel;};
	Eigen::MatrixXd ComputeForceDerivative(FEM::World* world);
	Eigen::VectorXd ComputeForce(FEM::World* world);
};
void MakeMuscles(const std::string& path,Octopus* ms);



#endif