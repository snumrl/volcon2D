#ifndef __OBJECT_H__
#define __OBJECT_H__
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include <vector>

struct Muscle
{
	double											activationLevel;
	std::vector<FEM::MuscleConstraint*>				muscleConstraints;	
	void SetActivationLevel(double a)
	{
		activationLevel = a;
		for(auto& mc: muscleConstraints)
			mc->SetActivationLevel(a);
	};
};

class Object
{
public:
	FEM::Mesh*											mMesh;
	std::vector<Muscle*>								mMuscles;
	Eigen::VectorXd										mActivationLevel;
	std::vector<FEM::Constraint*>						mConstraints;

	double												mMuscleStiffness;
	double												mYoungsModulus;
	double												mPoissonRatio;

	std::vector<FEM::AttachmentConstraint*>	mAttachementConstraintVector;
public:
	Object(const Eigen::Vector2d& trans);

	int GetNumMuscles() { return mMuscles.size();};
	std::vector<Muscle*>& GetMuscles(){return mMuscles;};
	void AddMuscle(
		const std::vector<Eigen::Vector3i>& indexList,
		const Eigen::Vector2d& fiber_direction
		);

	void Initialize(FEM::World* world);

	// void TransformAttachmentPoints();
	void SetActivationLevel(const Eigen::VectorXd& a);
	Eigen::VectorXd& GetActivationLevel(){return mActivationLevel;};
	// Eigen::MatrixXd ComputeForceDerivative(FEM::World* world);
	double ComputeForce(FEM::World* world);
};
void MakeMuscles(const std::string& path,Object* ms);



#endif