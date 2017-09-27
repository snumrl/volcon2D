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
#include <utility>

class Controller;
class MuscleOptimization;
class MusculoSkeletalSystem;
class Machine;
class IKOptimization;
class BallInfo;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;






class Controller
{
private:
	Machine*									mFSM;
	FEM::World*									mSoftWorld;
	dart::simulation::WorldPtr  				mRigidWorld;
	MusculoSkeletalSystem*						mMusculoSkeletalSystem;
	std::vector<BallInfo*>						mBalls;

	Eigen::VectorXd								mTargetPositions;
	Eigen::VectorXd								mTargetVelocities;
	Eigen::VectorXd 							mKp,mKv;

	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mMuscleOptimizationSolver;

public:
	Controller();
	void Initialize(FEM::World* soft_world,const dart::simulation::WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system,const std::vector<dart::dynamics::SkeletonPtr>& balls);
	Eigen::VectorXd Compute();
	Eigen::VectorXd ComputePDForces();
	const Eigen::VectorXd& GetTargetPositions(){return mTargetPositions;};
	// Eigen::VectorXd SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap);

	Machine* GetMachine();
};

#endif