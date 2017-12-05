#include "MusculoSkeletalDDP.h"
#include "../MuscleOptimization.h"
#include "../BallInfo.h"
using namespace Ipopt;
MusculoSkeletalDDP::
MusculoSkeletalDDP(	const dart::simulation::WorldPtr& rigid_world,
					FEM::World* soft_world,
					MusculoSkeletalSystem* musculo_skeletal_system,int n,int max_iteration)
					:DDP(musculo_skeletal_system->GetSkeleton()->getNumDofs()*2,
						musculo_skeletal_system->GetSkeleton()->getNumDofs()*2,n,max_iteration),
					mRigidWorld(rigid_world),mSoftWorld(soft_world),
					mMusculoSkeletalSystem(musculo_skeletal_system),mDofs(musculo_skeletal_system->GetSkeleton()->getNumDofs())
{
	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	mMuscleOptimizationSolver->Initialize();
	mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);

	double kp = 500.0;
	double kv = 2*sqrt(kp);
	int nn = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(nn,kp);
	mKv = Eigen::VectorXd::Constant(nn,kv);
	
	mTargetPositions = Eigen::VectorXd::Constant(nn,0.0);
	mTargetVelocities = Eigen::VectorXd::Constant(nn,0.0);
}

void
MusculoSkeletalDDP::
SetState(const Eigen::VectorXd& x)
{
	mMusculoSkeletalSystem->GetSkeleton()->setPositions(x.head(mDofs));
	mMusculoSkeletalSystem->GetSkeleton()->setVelocities(x.tail(mDofs));
	mMusculoSkeletalSystem->GetSkeleton()->computeForwardKinematics(true,true,false);
}
void
MusculoSkeletalDDP::
SetControl(const Eigen::VectorXd& u)
{
	mTargetPositions = u.head(mDofs);
	mTargetVelocities = u.tail(mDofs);
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
Step()
{
	auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	Eigen::VectorXd pos_diff = skel->getPositionDifferences(mTargetPositions,skel->getPositions());
	for(int i = 0;i<pos_diff.rows();i++)
		pos_diff[i] = dart::math::wrapToPi(pos_diff[i]);
	Eigen::VectorXd qdd_desired = pos_diff.cwiseProduct(mKp) + (mTargetVelocities - skel->getVelocities()).cwiseProduct(mKv);
	// std::cout<<qdd_desired.transpose()<<std::endl;
	// std::cout<<(skel->getMassMatrix()*qdd_desired+skel->getCoriolisAndGravityForces()).transpose()<<std::endl;
	// static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);
	// mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);

	// auto solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();

	// mMusculoSkeletalSystem->SetActivationLevel(solution.tail(mMusculoSkeletalSystem->GetNumMuscles()));
	// mMusculoSkeletalSystem->TransformAttachmentPoints();
	// mSoftWorld->TimeStepping(false);


	double nn = mSoftWorld->GetTimeStep() / mRigidWorld->getTimeStep();
	for(int i =0; i<nn;i++)
	{
		// mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);
		
		// std::cout<<"pos : "<<skel->getPositions().transpose()<<std::endl;
		// std::cout<<"vel : "<<skel->getVelocities().transpose()<<std::endl;
		// std::cout<<"qdd"<<qdd_desired.transpose()<<std::endl;
		// std::cout<<"getMassMatrix : "<<skel->getMassMatrix()<<std::endl;
		// std::cout<<"getCoriolisAndGravityForces : "<<skel->getCoriolisAndGravityForces()<<std::endl;
		// std::cout<<"setForce"<<(skel->getMassMatrix()*qdd_desired+skel->getCoriolisAndGravityForces()).transpose()<<std::endl;
		// std::cout<<"qdd"<<qdd_desired.transpose()<<std::endl;

		// std::cout<<std::endl;

		skel->setForces(skel->getMassMatrix()*qdd_desired+skel->getCoriolisAndGravityForces());

		mRigidWorld->step();
		// std::cout<<skel->getPositions().transpose()<<std::endl;

	}
}