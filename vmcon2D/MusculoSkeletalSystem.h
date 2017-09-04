#ifndef __MUSCULOSKELETAL_SYSTEM_H__
#define __MUSCULOSKELETAL_SYSTEM_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include <vector>


typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
Eigen::Vector2d GetPoint(const AnchorPoint& p);
struct Muscle
{
	FEM::Mesh*								mesh;
	std::vector<AnchorPoint> 				originWayPoints,insertionWayPoints;

	FEM::AttachmentConstraint*				origin;
	FEM::AttachmentConstraint*				insertion;
	std::vector<FEM::MuscleConstraint*>		muscleConstraints;
	std::vector<FEM::Constraint*>			constraints;
	double									activationLevel;
	Eigen::Vector2d							force_origin,force_insertion;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void TransferForce(Eigen::Vector2d& f_origin,Eigen::Vector2d& f_insertion);
};

class MusculoSkeletalSystem
{
private:
	std::vector<Muscle*>				mMuscles;
	dart::dynamics::SkeletonPtr			mSkeleton;
	Eigen::VectorXd						mActivationLevel;
	double								mTendonStiffness;
	double								mMuscleStiffness;
	double								mYoungsModulus;
	double								mPoissonRatio;
	std::vector<FEM::AttachmentConstraint*>	mAttachementConstraintVector;
public:
	MusculoSkeletalSystem();

	int GetNumMuscles() { return mMuscles.size();};
	std::vector<Muscle*>& GetMuscles(){return mMuscles;};
	void AddMuscle(
		const std::vector<AnchorPoint>& origin,
		const int&			origin_i,
		const std::vector<AnchorPoint>& insertion,
		const int&			insertion_i,
		const Eigen::Vector2d& fiber_direction,
		FEM::Mesh* mesh
		);

	void Initialize(FEM::World* world);

	void SetSkeleton(dart::dynamics::SkeletonPtr& skel){mSkeleton = skel;};
	dart::dynamics::SkeletonPtr& GetSkeleton(){return mSkeleton;};

	void TransformAttachmentPoints();
	void SetActivationLevel(const Eigen::VectorXd& a);
	Eigen::VectorXd& GetActivationLevel(){return mActivationLevel;};
	void ApplyForcesToSkeletons(FEM::World* world);
	Eigen::MatrixXd ComputeForceDerivative(FEM::World* world);
	Eigen::VectorXd ComputeForce(FEM::World* world);
};
void MakeMuscles(const std::string& path,MusculoSkeletalSystem* ms);
void MakeSkeleton(MusculoSkeletalSystem* ms);
#endif
