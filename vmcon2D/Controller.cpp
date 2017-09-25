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

BallInfo::
BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel)
	:isReleased(false),constraint(cons),skeleton(skel),releasedPoint(Eigen::Vector3d::Zero()),releasedVelocity(Eigen::Vector3d::Zero()),releaseFrameCount(0)
{
}
void
BallInfo::
ComputeFallingPosition(double h,Eigen::Vector3d& fp)
{
	Eigen::Vector3d p = skeleton->getBodyNode(0)->getCOM();
	Eigen::Vector3d v = skeleton->getBodyNode(0)->getCOMLinearVelocity();

	double dx = h-p[1];
	double g = -9.8;
	double v2_plus_2gdx = v[1]*v[1] + 2.0*g*dx;
	if(v2_plus_2gdx<0)
	{
		std::cout<<"no solution"<<std::endl;
		return;
	}

	double t1 = (-v[1] - sqrt(v2_plus_2gdx))/g;
	double t2 = (-v[1] + sqrt(v2_plus_2gdx))/g;
	double t = (t1>t2?t1:t2);

	fp = p+t*v;
	fp[1] = h;
}
void
BallInfo::
GetPosition(Eigen::Vector3d& p)
{
	p =  skeleton->getBodyNode(0)->getCOM();
}
void
BallInfo::
Release(const WorldPtr& world)
{
	if(!isReleased){
		world->getConstraintSolver()->removeConstraint(constraint);
		isReleased = true;
	}

	releaseFrameCount = 10;
	releasedPoint = skeleton->getBodyNode(0)->getCOM();
	releasedVelocity = skeleton->getBodyNode(0)->getCOMLinearVelocity();
	releasedImpulse = skeleton->getBodyNode(0)->getConstraintImpulse()*1000;
	std::cout<<releasedImpulse.transpose()<<std::endl;
}

void
BallInfo::
Attach(const WorldPtr& world,BodyNode* bn)
{
	constraint.reset();
	constraint = std::make_shared<dart::constraint::WeldJointConstraint>(skeleton->getBodyNode(0),bn);
	isReleased = false;
	world->getConstraintSolver()->addConstraint(constraint);
}
void
BallInfo::
TimeStepping()
{
	if(releaseFrameCount>0)
	{
		skeleton->getBodyNode(0)->addExtForce(releasedImpulse.tail(3),Eigen::Vector3d(0,0,0));
	}
	releaseFrameCount--;
}
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
		auto* abn =mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");
		Eigen::Vector3d loc = abn->getTransform().translation();
		balls[i]->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));

		mBalls.push_back(new BallInfo(std::make_shared<dart::constraint::WeldJointConstraint>(balls[i]->getBodyNode(0),abn),balls[i]));
		mBalls.back()->Attach(mRigidWorld,abn);
	}

	

	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	mIKOptimization = new IKOptimization(mMusculoSkeletalSystem->GetSkeleton());
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	mIKSolver = new IpoptApplication();

	mIKSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mIKSolver->Options()->SetStringValue("jac_c_constant", "yes");
	mIKSolver->Options()->SetStringValue("hessian_constant", "yes");
	mIKSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mIKSolver->Options()->SetIntegerValue("print_level", 2);
	mIKSolver->Options()->SetIntegerValue("max_iter", 1000);
	mIKSolver->Options()->SetNumericValue("tol", 1e-4);

	mRestPose = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	mPreviousPose = mRestPose;
	mTargetPositions = mRestPose;

	double kp = 400.0;
	double kv = 2*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);

	double reinforce_ratio = 2.0;
	for(int i =0;i<n;i++)
		if(!mMusculoSkeletalSystem->GetSkeleton()->getDof(i)->getChildBodyNode()->getName().compare("ElbowR"))
		{
			mKp[i] *=reinforce_ratio;
			mKv[i] *= sqrt(reinforce_ratio);
		}
		else if(!mMusculoSkeletalSystem->GetSkeleton()->getDof(i)->getChildBodyNode()->getName().compare("ElbowL"))
		{
			mKp[i] *=reinforce_ratio;
			mKv[i] *= sqrt(reinforce_ratio);
		}
	mFSM = new Machine();
	MakeMachine("../vmcon2D/export/juggling.xml",mFSM);
	mFSM->Trigger("start_right");
}
const Eigen::VectorXd&
Controller::
GetTargetPositions()
{
	return mTargetPositions;
}
void
Controller::
SetTargetPositions(const Eigen::VectorXd& tp)
{
	mTargetPositions = tp;
}
Eigen::VectorXd
Controller::
Compute()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();
	Eigen::VectorXd qdd_desired = ComputePDForces();

	static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);
	mBalls[0]->TimeStepping();
	if(mFSM->GetStates()["LEFT_CATCH"] == mFSM->GetCurrentState() && mFSM->GetCurrentState()->GetTime()==0.0)
		mBalls[0]->Release(mRigidWorld);
	else if(mFSM->GetStates()["RIGHT_CATCH"] == mFSM->GetCurrentState() && mFSM->GetCurrentState()->GetTime()==0.0)
		mBalls[0]->Release(mRigidWorld);

	if(mFSM->GetStates()["LEFT_CATCH"] == mFSM->GetCurrentState())
	{
		Eigen::Vector3d ball_position,left_hand_position;
		mBalls[0]->GetPosition(ball_position);
		left_hand_position = skel->getBodyNode("HandL")->getTransform().translation();
		if((left_hand_position-ball_position).norm()<5E-2){
			mBalls[0]->Attach(mRigidWorld,skel->getBodyNode("HandL"));
			mFSM->Trigger("catch");
		}
	}
	else if(mFSM->GetStates()["RIGHT_CATCH"] == mFSM->GetCurrentState())
	{
		Eigen::Vector3d ball_position,left_hand_position;
		mBalls[0]->GetPosition(ball_position);
		left_hand_position = skel->getBodyNode("HandR")->getTransform().translation();
		if((left_hand_position-ball_position).norm()<5E-2){
			mBalls[0]->Attach(mRigidWorld,skel->getBodyNode("HandR"));
			mFSM->Trigger("catch");
		}
	}



	if(mSoftWorld->GetTime() == 0.0){
		mMuscleOptimizationSolver->Initialize();
		mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	
	}
	else
		mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	

	Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
	Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
	Eigen::VectorXd activation = solution.tail(mMusculoSkeletalSystem->GetNumMuscles()); 
	// 
	return activation;
}
Eigen::VectorXd
Controller::
ComputePDForces()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

	mFSM->GetCurrentState()->SetTime(mFSM->GetCurrentState()->GetTime()+mSoftWorld->GetTimeStep());
	mFSM->GetMotion(mTargetPositions);

	if(mFSM->GetCurrentState()==mFSM->GetStates()["LEFT_CATCH"])
	{
		Eigen::Vector3d target;
		mBalls[0]->ComputeFallingPosition(mBalls[0]->releasedPoint[1],target);
		AnchorPoint ap;
		ap = std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL"),Eigen::Vector3d(0,0,0));
		const auto& ik_targets = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetTargets();
		bool need_ik_update = true;
		for(auto& t : ik_targets)
		{
			if(!t.first.first->getName().compare("HandL"))
			{
				if((t.second-target).norm()<1E-2)
					need_ik_update = false;
			}
		}
		if(need_ik_update)
			mTargetPositions = SolveIK(target,ap);
		else
			mTargetPositions = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetSolution();
	}
	else if(mFSM->GetCurrentState()==mFSM->GetStates()["RIGHT_CATCH"])
	{
		Eigen::Vector3d target;
		mBalls[0]->ComputeFallingPosition(mBalls[0]->releasedPoint[1],target);
		AnchorPoint ap;
		ap = std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d(0,0,0));
		const auto& ik_targets = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetTargets();
		bool need_ik_update = true;
		for(auto& t : ik_targets)
		{
			if(!t.first.first->getName().compare("HandR"))
			{
				if((t.second-target).norm()<1E-2)
					need_ik_update = false;
			}
		}
		if(need_ik_update)
			mTargetPositions = SolveIK(target,ap);
		else
			mTargetPositions = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetSolution();
	}	
	Eigen::VectorXd pos_m = mTargetPositions;
	Eigen::VectorXd vel_m(pos_m.rows());
	vel_m.setZero();

	Eigen::VectorXd pos = skel->getPositions();
	Eigen::VectorXd vel = skel->getVelocities();

	// for(int i = 0;i<pos.rows();i++)
	// 	pos[i] = dart::math::wrapToPi(pos[i]);
	// for(int i = 0;i<pos.rows();i++)
	// 	pos_m[i] = dart::math::wrapToPi(pos_m[i]);

	Eigen::VectorXd pos_diff(pos.rows());
	pos_diff.setZero();

	// for(int i =0;i<pos.rows();i++)
	// {
	// 	double phi = pos[i];
	// 	double two_phi = 2*3.141592 + pos[i];

	// 	pos_diff[i] = pos_m[i] -(phi<two_phi?phi:two_phi);
	// }
	// std::cout<<pos.transpose()<<std::endl;
	// std::cout<<pos_m.transpose()<<std::endl;
	// std::cout<<pos_diff.transpose()<<std::endl;
	pos_diff = skel->getPositionDifferences(pos_m,pos);
	for(int i = 0;i<pos_diff.rows();i++)
		pos_diff[i] = dart::math::wrapToPi(pos_diff[i]);
	
	Eigen::VectorXd qdd_desired = 
				pos_diff.cwiseProduct(mKp)+
				(vel_m - vel).cwiseProduct(mKv);

	// std::cout<<skel->getInvMassMatrix()*mMusculoSkeletalSystem->GetSkeleton()->getConstraintForces()<<std::endl;
	// qdd_desired -= 10*skel->getInvMassMatrix()*skel->getConstraintForces();


	// if(mFSM->GetCurrentState()==mFSM->GetStates()["LEFT_THROW"])
	// {

		
	// }
	// if(mFSM->GetCurrentState()==mFSM->GetStates()["RIGHT_SWING"]&&)
	// {
	// 	Eigen::Vector3d target_linear_velocity(0,10,0);
	// 	Eigen::Vector3d target_force = 100*(target_linear_velocity-skel->getBodyNode("HandR")->getCOMLinearVelocity());
	// 	qdd_desired += skel->getInvMassMatrix()*(skel->getLinearJacobian(skel->getBodyNode("HandR")).transpose())*target_force;

	// }
	return qdd_desired;
}

Eigen::VectorXd
Controller::
SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap)
{
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	ik->AddTargetPositions(ap,target_position);

	mIKSolver->Initialize();
	mIKSolver->OptimizeTNLP(mIKOptimization);

	Eigen::VectorXd solution = ik->GetSolution();
	return solution;
}


Machine*
Controller::
GetMachine()
{
	return mFSM;
}

