#include "Controller.h"
#include "BallInfo.h"
#include "DART_Interface.h"
#include "MusculoSkeletalSystem.h"
#include "IKOptimization.h"
#include "DDP/VelocityControlDDP.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;



Controller::
Controller()
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

	double kp = 200.0;
	double kv = 2*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);

	mTargetPositions = Eigen::VectorXd::Constant(n,0.0);
	mTargetVelocities = Eigen::VectorXd::Constant(n,0.0);



	mDDPMusculoSkeletalSystem = new MusculoSkeletalSystem();

	mDDPRigidWorld = std::make_shared<World>();
	mDDPRigidWorld->setTimeStep(1.0/1000.0);
	mDDPRigidWorld->setGravity(Eigen::Vector3d(0,-9.81,0));
	mDDPSoftWorld = new FEM::World(
		// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/500.0,										//time_step
		50, 											//max_iteration	
		0.999											//damping_coeff
		);
	MakeSkeleton(mDDPMusculoSkeletalSystem);
	MakeMuscles("../vmcon2D/export/muscle_parameter.xml",mDDPMusculoSkeletalSystem);

	mDDPRigidWorld->addSkeleton(mDDPMusculoSkeletalSystem->GetSkeleton());
	for(int i =0;i<1;i++)
	{
		mDDPBallSkeleletons.push_back(Skeleton::create("Ball_"+std::to_string(i)));
		MakeBall(mDDPBallSkeleletons.back(),0.036,0.13);	
		auto pos = mDDPBallSkeleletons.back()->getPositions();
		pos.tail(3) = Eigen::Vector3d(i*0.1,-0.3,0);
		mDDPBallSkeleletons.back()->setPositions(pos);
	}
	for(auto& ball : mDDPBallSkeleletons)
		mDDPRigidWorld->addSkeleton(ball);
	mDDPMusculoSkeletalSystem->Initialize(mDDPSoftWorld);
	mDDPSoftWorld->Initialize();

	for(int i =0;i<mDDPBallSkeleletons.size();i++)
	{
		bool is_left_hand = i%2;
		if(is_left_hand)
		{
			auto* abn =mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL");
			Eigen::Vector3d loc = abn->getTransform().translation();
			mDDPBallSkeleletons[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

			mDDPBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(mDDPBallSkeleletons[i]->getBodyNode(0),abn),mDDPBallSkeleletons[i]));
			mDDPBalls.back()->Attach(mDDPRigidWorld,abn);
		}
		else
		{
			auto* abn =mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");
			Eigen::Vector3d loc = abn->getTransform().translation();
			mDDPBallSkeleletons[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

			mDDPBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(mDDPBallSkeleletons[i]->getBodyNode(0),abn),mDDPBallSkeleletons[i]));
			mDDPBalls.back()->Attach(mDDPRigidWorld,abn);	
		}
	
	}


	mDDPMuscleOptimization = new MuscleOptimization(mDDPSoftWorld,mDDPRigidWorld,mDDPMusculoSkeletalSystem);
	mDDPMuscleOptimizationSolver = new IpoptApplication();
	mDDPMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mDDPMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mDDPMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mDDPMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mDDPMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mDDPMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mDDPMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	u_index = 0;
	std::vector<Eigen::VectorXd> u0;
	auto X_rigid = mDDPMusculoSkeletalSystem->GetSkeleton()->getPositions();
	auto X_soft = mDDPSoftWorld->GetPositions();

	ComputeInitialU0(u0);
	mDDPMusculoSkeletalSystem->GetSkeleton()->setPositions(X_rigid);
	mDDPMusculoSkeletalSystem->GetSkeleton()->setVelocities(mDDPMusculoSkeletalSystem->GetSkeleton()->getVelocities().setZero());
	mDDPSoftWorld->SetPositions(X_soft);

	mDDP = new VelocityControlDDP(mDDPRigidWorld,mDDPSoftWorld,mDDPMusculoSkeletalSystem,mDDPBalls[0],u0,u0.size(),40);
	mU = mDDP->Solve();
}
void
Controller::
ComputeInitialU0(std::vector<Eigen::VectorXd>& u0)
{
	int n = 100;
	u0.resize(n,Eigen::VectorXd::Zero(mDDPMusculoSkeletalSystem->GetNumMuscles()));

	
	// while(u0.size() != n)
	// {
	// 	bool is_fem_updated = false;
	// 	auto& skel =mDDPMusculoSkeletalSystem->GetSkeleton();


	// 	//Simulation Loop
	// 	if(mDDPSoftWorld->GetTime()<=mDDPRigidWorld->getTime())
	// 	{
	// 		is_fem_updated =true;

	// 		auto ui = Compute(false);
	// 		mDDPMusculoSkeletalSystem->SetActivationLevel(ui);
	// 		u0.push_back(ui);
	// 	}		
		
	// 	if(is_fem_updated)
	// 	{
	// 		mDDPMusculoSkeletalSystem->TransformAttachmentPoints();
	// 		mDDPSoftWorld->TimeStepping();
	// 	}
	// 	mDDPMusculoSkeletalSystem->ApplyForcesToSkeletons(mDDPSoftWorld);

	// 	mDDPRigidWorld->step();
	// }
}
Eigen::VectorXd
Controller::
Compute(bool test)
{
	if(test)
	{
		if(u_index<mU.size())
		{
			return mU[u_index++];
		}
		else if(u_index==mU.size()){
			mBalls[0]->Release(mRigidWorld);
			Eigen::Vector3d vv;
			mBalls[0]->GetVelocity(vv);
			std::cout<<vv.transpose()<<std::endl;
			u_index++;
		}
		
		auto& skel =mMusculoSkeletalSystem->GetSkeleton();
		Eigen::VectorXd qdd_desired = ComputePDForces();
		static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

		// if(mSoftWorld->GetTime() == 0.0){
			mMuscleOptimizationSolver->Initialize();
			mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	
		// }
		// else
			// mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	

		Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
		Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
		Eigen::VectorXd activation = solution.tail(mMusculoSkeletalSystem->GetNumMuscles());
		activation.setZero();
		return activation;
	}
	else
	{
		auto& skel =mDDPMusculoSkeletalSystem->GetSkeleton();
		Eigen::VectorXd qdd_desired = ComputePDForces(false);
		static_cast<MuscleOptimization*>(GetRawPtr(mDDPMuscleOptimization))->Update(qdd_desired);
		std::cout<<qdd_desired.norm()<<std::endl;
		if(mSoftWorld->GetTime() == 0.0){
			mDDPMuscleOptimizationSolver->Initialize();
			mDDPMuscleOptimizationSolver->OptimizeTNLP(mDDPMuscleOptimization);	
		}
		else
			mDDPMuscleOptimizationSolver->ReOptimizeTNLP(mDDPMuscleOptimization);	

		Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mDDPMuscleOptimization))->GetSolution();
		Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
		Eigen::VectorXd activation = solution.tail(mDDPMusculoSkeletalSystem->GetNumMuscles());

		return activation;
	}
	
}
Eigen::VectorXd
Controller::
ComputePDForces(bool test)
{
	if(test)
	{
		auto& skel =mMusculoSkeletalSystem->GetSkeleton();

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
	else
	{
		auto& skel =mDDPMusculoSkeletalSystem->GetSkeleton();
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
	
}

