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
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;


class BallInfo
{
public:
	dart::dynamics::SkeletonPtr 				skeleton;
	dart::constraint::WeldJointConstraintPtr 	constraint;

	bool 										isReleased;

	Eigen::Vector3d								releasedPoint;
	Eigen::Vector3d								releasedVelocity;
	Eigen::Vector6d								releasedImpulse;
	int 										releaseFrameCount;
public:
	BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel);

	void ComputeFallingPosition(double h,Eigen::Vector3d& fp);
	void Release(const dart::simulation::WorldPtr& world);
	void GetPosition(Eigen::Vector3d& p);
	void Attach(const dart::simulation::WorldPtr& world,dart::dynamics::BodyNode* bn);
	void TimeStepping();
};
class Controller
{
private:
	Machine*									mFSM;
	FEM::World*									mSoftWorld;
	dart::simulation::WorldPtr  				mRigidWorld;
	MusculoSkeletalSystem*						mMusculoSkeletalSystem;
	std::vector<BallInfo*>						mBalls;
	Eigen::VectorXd								mTargetPositions;

	Eigen::VectorXd 							mKp,mKv;

	Eigen::VectorXd								mRestPose;
	Eigen::VectorXd								mPreviousPose;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mMuscleOptimizationSolver;

	Ipopt::SmartPtr<Ipopt::TNLP>			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;
public:
	Controller();
	void Initialize(FEM::World* soft_world,const dart::simulation::WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system,const std::vector<dart::dynamics::SkeletonPtr>& balls);
	const Eigen::VectorXd& GetTargetPositions();
	void SetTargetPositions(const Eigen::VectorXd& tp);
	Eigen::VectorXd Compute();
	Eigen::VectorXd ComputePDForces();
	Eigen::VectorXd SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap);

	Machine* GetMachine();
};

#endif