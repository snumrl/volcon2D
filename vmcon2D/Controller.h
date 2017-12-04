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
class VelocityControlDDP;
class IKOptimization;
class BallInfo;
class BezierCurve;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;

class Controller
{
private:
	FEM::World*									mSoftWorld;
	dart::simulation::WorldPtr  				mRigidWorld;
	MusculoSkeletalSystem*						mMusculoSkeletalSystem;
	std::vector<BallInfo*>						mBalls;

	Eigen::VectorXd								mTargetPositions;
	Eigen::VectorXd								mTargetVelocities;
	Eigen::VectorXd 							mKp,mKv;

	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mMuscleOptimizationSolver;


	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;

	FEM::World*					mDDPSoftWorld;
	dart::simulation::WorldPtr  mDDPRigidWorld;
	MusculoSkeletalSystem*		mDDPMusculoSkeletalSystem;
	std::vector<BallInfo*>						mDDPBalls;
	std::vector<dart::dynamics::SkeletonPtr>	mDDPBallSkeleletons;

	VelocityControlDDP*			mDDP;
	std::vector<Eigen::VectorXd> mU;
	int u_index;

	BezierCurve* mBezierCurve;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mDDPMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mDDPMuscleOptimizationSolver;
public:
	Controller();
	void Initialize(FEM::World* soft_world,const dart::simulation::WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system,const std::vector<dart::dynamics::SkeletonPtr>& balls);
	void ComputeInitialU0(std::vector<Eigen::VectorXd>& u0);
	Eigen::VectorXd Compute(bool test = true);
	Eigen::VectorXd ComputePDForces(bool test = true);
	const Eigen::VectorXd& GetTargetPositions(){return mTargetPositions;};
};

#endif