#include "MusculoSkeletalDDP.h"

MusculoSkeletalDDP::
MusculoSkeletalDDP(	const dart::simulation::WorldPtr& rigid_world,
					FEM::World* soft_world,
					MusculoSkeletalSystem* musculo_skeletal_system,int n,int max_iteration)
					:DDP(musculo_skeletal_system->GetSkeleton()->getNumDofs()*2,musculo_skeletal_system->GetNumMuscles(),n,max_iteration),
					mRigidWorld(rigid_world),mSoftWorld(soft_world),
					mMusculoSkeletalSystem(musculo_skeletal_system),mDofs(musculo_skeletal_system->GetSkeleton()->getNumDofs())
{

}

void
MusculoSkeletalDDP::
SetState(const Eigen::VectorXd& x)
{
	mMusculoSkeletalSystem->GetSkeleton()->setPositions(x.head(mDofs));
	mMusculoSkeletalSystem->GetSkeleton()->setVelocities(x.tail(mDofs));
	mMusculoSkeletalSystem->GetSkeleton()->computeForwardKinematics(true,false,false);
}
void
MusculoSkeletalDDP::
SetControl(const Eigen::VectorXd& u)
{
	mMusculoSkeletalSystem->SetActivationLevel(u);
}
void
MusculoSkeletalDDP::
GetState(Eigen::VectorXd& x)
{
	x.head(mDofs) = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	x.tail(mDofs) = mMusculoSkeletalSystem->GetSkeleton()->getVelocities();
}

void
MusculoSkeletalDDP::
Step(bool fem_update)
{
	if(fem_update)
	{
		mMusculoSkeletalSystem->TransformAttachmentPoints();
		mSoftWorld->TimeStepping(false);	
	}
	

	double nn = mSoftWorld->GetTimeStep() / mRigidWorld->getTimeStep();
	for(int i =0; i<nn;i++)
	{
		if(fem_update)
			mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);
		mRigidWorld->step();
	}
}