#ifndef __MUSCULO_SKELETAL_DDP_H__
#define __MUSCULO_SKELETAL_DDP_H__
#include "DDP.h"

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
	void Step(bool fem_update = true);
	int mDofs;
	dart::simulation::WorldPtr mRigidWorld;
	FEM::World*					mSoftWorld;
	MusculoSkeletalSystem* mMusculoSkeletalSystem;
};

#endif