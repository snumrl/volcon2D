#include "Controller.h"
#include "BezierCurve.h"
#include "BallInfo.h"
#include "DART_Interface.h"
#include "MusculoSkeletalSystem.h"
#include "IKOptimization.h"
#include "iLQR/MusculoSkeletalLQR.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include <fstream>
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

	mMuscleOptimizationSolver->Initialize();
	mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	

	double kp = 1000.0;
	double kv = 2.0*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);

	mTargetPositions = Eigen::VectorXd::Constant(n,0.0);
	mTargetPositions2 = Eigen::VectorXd::Constant(n,0.0);
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
		1.0/250.0,										//time_step
		100, 											//max_iteration	
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
			auto* abn =mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");
			Eigen::Vector3d loc = abn->getTransform().translation();
			mDDPBallSkeleletons[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

			mDDPBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(mDDPBallSkeleletons[i]->getBodyNode(0),abn),mDDPBallSkeleletons[i]));
			mDDPBalls.back()->Attach(mDDPRigidWorld,abn);	
		}
	
	}

	std::vector<Eigen::VectorXd> u0;
	Eigen::VectorXd X_rigid = mDDPMusculoSkeletalSystem->GetSkeleton()->getPositions();
	Eigen::VectorXd X_soft = mDDPSoftWorld->GetPositions();

	double xx,yy,zz,vx,vy,vz;
	std::ifstream target("../vmcon2D/export/target.txt");
	target>>xx>>yy>>zz>>vx>>vy>>vz;
	target.close();
	ComputeInitialU0(u0,Eigen::Vector3d(xx,yy,zz),Eigen::Vector3d(vx,vy,vz));
	// for(int i=0;i<49;i++)
	// 	u0.push_back(X_rigid);
	mDDPMusculoSkeletalSystem->GetSkeleton()->setPositions(X_rigid);
	mDDPMusculoSkeletalSystem->GetSkeleton()->setVelocities(mDDPMusculoSkeletalSystem->GetSkeleton()->getVelocities().setZero());
	mDDPMusculoSkeletalSystem->GetSkeleton()->computeForwardKinematics(true,false,false);
	// mDDPSoftWorld->SetPositions(X_soft);

	Eigen::VectorXd x0(mDDPMusculoSkeletalSystem->GetSkeleton()->getNumDofs()*2);
	x0.head(mDDPMusculoSkeletalSystem->GetSkeleton()->getNumDofs()) = X_rigid;
	x0.tail(mDDPMusculoSkeletalSystem->GetSkeleton()->getNumDofs()).setZero();


	mLQR = new MusculoSkeletalLQR(Eigen::Vector3d(xx,yy,zz),Eigen::Vector3d(vx,vy,vz),mDDPRigidWorld,mDDPSoftWorld,mDDPMusculoSkeletalSystem,u0.size()+1,200);
	mLQR->Initialze(x0,u0);
	mU = mLQR->Solve();
	u_index = 0;
	// mU = u0;
	
	
	// for(int i =0;i<u0.size();i++)
	// 	std::cout<<u0[i]<<" ";
	// std::cout<<std::endl;
	// mU = u0;
	for(int i =0;i<mU.size();i++)
		std::cout<<mU[i].transpose()<<std::endl;
	std::cout<<std::endl;
	// std::cout<<std::endl<<std::endl;
}
void
Controller::
ComputeInitialU0(std::vector<Eigen::VectorXd>& u0,const Eigen::Vector3d& target_p,const Eigen::Vector3d& target_v)
{
	int n = 49;
	// u0.resize(n);
	double t = mDDPSoftWorld->GetTimeStep()*n;
	Eigen::Vector2d p0,p1,p2,v2;
	v2 = Eigen::Vector2d(1,1.0);
	p0 = mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM().block<2,1>(0,0);
	p2 = target_p.block<2,1>(0,0);
	p1 = p2-0.5*target_v.block<2,1>(0,0)*t;
	std::cout<<p0.transpose()<<std::endl;
	std::cout<<p1.transpose()<<std::endl;
	std::cout<<p2.transpose()<<std::endl;
	mBezierCurve = new BezierCurve(p0,p1,p2,t);
	// 
	mIKOptimization = new IKOptimization(mDDPMusculoSkeletalSystem->GetSkeleton());

	mIKSolver = new IpoptApplication();
	mIKSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mIKSolver->Options()->SetStringValue("jac_c_constant", "yes");
	mIKSolver->Options()->SetStringValue("hessian_constant", "yes");
	mIKSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mIKSolver->Options()->SetIntegerValue("print_level", 2);
	mIKSolver->Options()->SetIntegerValue("max_iter", 1000);
	mIKSolver->Options()->SetNumericValue("tol", 1e-4);

	mIKSolver->Initialize();
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));

	// // Eigen::Vector3d l_loc = mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL")->getCOM();
	Eigen::Vector3d r_loc = mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM();
	// // ik->AddTargetPositions(std::make_pair(mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL"),Eigen::Vector3d::Zero()),l_loc);
	ik->AddTargetPositions(std::make_pair(mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d::Zero()),r_loc);

	mIKSolver->OptimizeTNLP(mIKOptimization);
	
	AnchorPoint ap = std::make_pair(mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d(0,0,0));
	std::vector<std::pair<Eigen::VectorXd,double>> motions;
	int num_samples = 10;
	for(int i =0;i<num_samples+3;i++)
	{
		double cur_t = (double)i/(double)num_samples*t;
		Eigen::Vector2d p_ee = mBezierCurve->GetPosition(cur_t);
		Eigen::Vector3d p_ee_3d(p_ee[0],p_ee[1],0);

		ik->AddTargetPositions(ap,p_ee_3d);
		mIKSolver->ReOptimizeTNLP(mIKOptimization);
		motions.push_back(std::make_pair(ik->GetSolution(),cur_t));
	}

	for(int i=0;i<n;i++)
	{
		double cur_t = mDDPSoftWorld->GetTimeStep()*i;
		int k =0,k1 =0;
		for(int i =0;i<motions.size();i++)
		{
			if(motions[i].second<cur_t)
				k=i;
		}
		k1 = k+1;
		double kt = cur_t - motions[k].second;
		double dt = motions[k1].second - motions[k].second;

		kt /= dt;

		auto p = (1.0 - kt)*motions[k].first + kt*motions[k1].first;
		double kt1 = kt + 0.01;
		auto pdp = (1.0 - kt1)*motions[k].first + kt1*motions[k1].first;
		mInitialPositions.push_back(p);

	}

	for(int i =0;i<n;i++)
	{
		Eigen::VectorXd compose(
			mInitialPositions[i].rows() );
		// compose.setZero();
		compose.head(mInitialPositions[i].rows()) = mInitialPositions[i];
		// compose[mInitialPositions[i].rows()] = 1.0;
		// // compose.head(mInitialPositions[i].rows()) = mInitialPositions[i];
		// // compose.block(mInitialPositions[i].rows(),0,mInitialPositions[i].rows(),1) = mInitialVelocities[i];
		// compose.tail(mInitialVelocities[i].rows()) = mKp;
		// u0.push_back(mInitialPositions[i]);
		u0.push_back(compose);
	}
	

}
Eigen::VectorXd
Controller::
Compute()
{	
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

	
	Eigen::VectorXd qdd_desired = ComputePDForces();
	static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

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
	if(u_index == mU.size())
	{
			Eigen::Vector3d vv;
		mBalls[0]->GetVelocity(vv);
		// std::cout<<"Release : "<<vv.transpose()<<std::endl;
		std::cout<<"Release : "<<mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOMLinearVelocity().transpose()<<std::endl;
		mBalls[0]->Release(mRigidWorld);
		u_index++;
	}

	if(u_index<mU.size())
	{
		// IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
		// Eigen::Vector3d new_target;
		// auto tar = ik->GetTargets();
		// new_target = tar[0].first.second;
		// new_target[0] = 0.338854;
		// new_target[1] = mU[u_index][0];
		// ik->AddTargetPositions(std::make_pair(mDDPMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d::Zero()),new_target);
		// mIKSolver->ReOptimizeTNLP(mIKOptimization);
		// mTargetPositions = ik->GetSolution();
		// std::cout<<new_target.transpose()<<std::endl;
		// std::cout<<"mTargetPositions : "<<mTargetPositions.transpose()<<std::endl;
		mTargetPositions = mU[u_index].head(skel->getNumDofs());
		// double kp = 1000.0*mU[u_index][skel->getNumDofs()];
		// double kv = 2.0*sqrt(kp);
		// mKp = Eigen::VectorXd::Constant(skel->getNumDofs(),kp);
		// mKv = Eigen::VectorXd::Constant(skel->getNumDofs(),kv);
		mTargetPositions2 = mInitialPositions[u_index];

		if(u_index == 0)
			mTargetVelocities.setZero();
		else
			mTargetVelocities =  (mU[u_index]-mU[u_index-1]).head(skel->getNumDofs())/mSoftWorld->GetTimeStep();
		// mTargetVelocities = mU[u_index].block(skel->getNumDofs(),0,skel->getNumDofs(),1);
		// mKp = mU[u_index].tail(skel->getNumDofs());
		// mKv = 2*mKp.cwiseSqrt();
		u_index++;
	}
	else
	{
		mTargetPositions = mU.back().head(skel->getNumDofs());
		mTargetPositions2 = mInitialPositions.back();
		mTargetVelocities.setZero();
	}
	

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

