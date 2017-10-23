#include "Controller.h"
#include "MusculoSkeletalSystem.h"
#include "IKOptimization.h"
#include "FSM_Interface.h"
#include "FSM.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;



Controller::
Controller()
	:mFSM(nullptr)
{
	
}

void
Controller::
Initialize(FEM::World* soft_world,const WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system,const std::vector<SkeletonPtr>& balls)
{
	mSoftWorld = soft_world;
	mRigidWorld = rigid_world;
	mMusculoSkeletalSystem = musculo_skeletal_system;
	for(int i =0;i<balls.size();i++)
	{
		bool is_left_hand = i%2;
		if(is_left_hand)
		{
			auto* abn =mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL");
			Eigen::Vector3d loc = abn->getTransform().translation();
			balls[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

			mBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(balls[i]->getBodyNode(0),abn),balls[i]));
			mBalls.back()->Attach(mRigidWorld,abn);
		}
		else
		{
			auto* abn =mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");
			Eigen::Vector3d loc = abn->getTransform().translation();
			balls[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

			mBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(balls[i]->getBodyNode(0),abn),balls[i]));
			mBalls.back()->Attach(mRigidWorld,abn);	
		}
	
	}

	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	double kp = 4000.0;
	double kv = 2*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);

	mFSM = new Machine(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mBalls,this,mSoftWorld->GetTimeStep());
	MakeMachine("../vmcon2D/export/juggling.xml",mFSM);
	mFSM->Trigger("start");
}

Eigen::VectorXd
Controller::
Compute()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();
	Eigen::VectorXd qdd_desired = ComputePDForces();

	static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

	for(auto& ball : mBalls)
		ball->TimeStepping();




	if(mSoftWorld->GetTime() == 0.0){
		mMuscleOptimizationSolver->Initialize();
		mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	
	}
	else
		mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	

	Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
	Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
	Eigen::VectorXd activation = solution.tail(mMusculoSkeletalSystem->GetNumMuscles()); 

	return activation;
}
Eigen::VectorXd
Controller::
ComputePDForces()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

	mFSM->GetMotion(mTargetPositions,mTargetVelocities);


	Eigen::VectorXd pos_m = mTargetPositions;
	Eigen::VectorXd vel_m = mTargetVelocities;

	Eigen::VectorXd pos = skel->getPositions();
	Eigen::VectorXd vel = skel->getVelocities();


	Eigen::VectorXd pos_diff(pos.rows());

	pos_diff = skel->getPositionDifferences(pos_m,pos);
	for(int i = 0;i<pos_diff.rows();i++)
		pos_diff[i] = dart::math::wrapToPi(pos_diff[i]);
	
	Eigen::VectorXd qdd_desired = 
				pos_diff.cwiseProduct(mKp)+
				(vel_m - vel).cwiseProduct(mKv);

	return qdd_desired;
}



Machine*
Controller::
GetMachine()
{
	return mFSM;
}

