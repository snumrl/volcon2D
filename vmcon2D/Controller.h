#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>	
class Controller;
class MuscleOptimization;
class MusculoSkeletalSystem;
class Machine;
class IKOptimization;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
class Controller
{
private:
	Machine*		mFSM;
	FEM::World*					mSoftWorld;
	dart::simulation::WorldPtr  mRigidWorld;
	MusculoSkeletalSystem*		mMusculoSkeletalSystem;
	Eigen::VectorXd				mTargetPositions;

	Eigen::VectorXd 			mKp,mKv;

	Eigen::VectorXd							mRestPose;
	Eigen::VectorXd							mPreviousPose;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> mMuscleOptimizationSolver;

	Ipopt::SmartPtr<Ipopt::TNLP>			 mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> mIKSolver;
public:
	Controller();
	void Initialize(FEM::World* soft_world,const dart::simulation::WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system);
	const Eigen::VectorXd& GetTargetPositions();
	void SetTargetPositions(const Eigen::VectorXd& tp);
	Eigen::VectorXd ComputeActivationLevels();
	Eigen::VectorXd ComputePDForces();
	Eigen::VectorXd SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap);

	Machine* GetMachine();
};

#endif