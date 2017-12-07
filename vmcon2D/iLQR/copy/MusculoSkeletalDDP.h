#ifndef __MUSCULO_SKELETAL_DDP_H__
#define __MUSCULO_SKELETAL_DDP_H__
#include "DDP.h"
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>	
#include "../MusculoSkeletalSystem.h"
class MusculoSkeletalSystem;
class MusculoSkeletalDDP : public DDP
{
public:
	MusculoSkeletalDDP(const dart::simulation::WorldPtr& rigid_world,
		FEM::World* soft_world,
		MusculoSkeletalSystem* mMusculoSkeletalSystem,int n,int max_iteration);

protected:
	void SetState(const Eigen::VectorXd& x);
	void SetControl(const Eigen::VectorXd& u);
	void GetState(Eigen::VectorXd& x);
	void Step();
	int 										mDofs;
	dart::simulation::WorldPtr 					mRigidWorld;
	FEM::World*									mSoftWorld;
	Eigen::VectorXd								mSoftWorldPositions;
	MusculoSkeletalSystem* 						mMusculoSkeletalSystem;

	Eigen::VectorXd								mTargetPositions,mTargetVelocities;
	Eigen::VectorXd 							mKp,mKv;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mMuscleOptimizationSolver;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;
};

#endif