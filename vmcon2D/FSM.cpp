#include "FSM.h"
#include "BezierCurve.h"
#include "IKOptimization.h"
#include "MusculoSkeletalSystem.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;
BallInfo::
BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel)
	:isReleased(true),constraint(cons),skeleton(skel),releasedPoint(Eigen::Vector3d::Zero()),releasedVelocity(Eigen::Vector3d::Zero()),releaseFrameCount(0)
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
	if(t<0.0)
		t=0.05;

	fp = p+t*v;
	if(fp[1]>h)
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
GetVelocity(Eigen::Vector3d& v)
{
	v = skeleton->getBodyNode(0)->getCOMLinearVelocity();
}
void
BallInfo::
Release(const WorldPtr& world)
{
	if(!isReleased){
		world->getConstraintSolver()->removeConstraint(constraint);
		isReleased = true;
	}

	releaseFrameCount = 0;
	releasedPoint = skeleton->getBodyNode(0)->getCOM();
	releasedVelocity = skeleton->getBodyNode(0)->getCOMLinearVelocity();
	releasedImpulse = skeleton->getBodyNode(0)->getConstraintImpulse()*1000;
	
	std::cout<<"Released Velocity : "<<releasedVelocity.transpose()<<std::endl;
}

void
BallInfo::
Attach(const WorldPtr& world,BodyNode* bn)
{
	if(isReleased)
	{
		constraint.reset();
		constraint = std::make_shared<dart::constraint::WeldJointConstraint>(skeleton->getBodyNode(0),bn);
		isReleased = false;
		// Eigen::Vector3d loc = bn->getCOM();
		// skeleton->setPositions(Eigen::compose(Eigen::Vector3d(0,0,0),loc));
		world->getConstraintSolver()->addConstraint(constraint);

	}
	
	
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




State::
State(const WorldPtr& rigid_world,
	const SkeletonPtr& skeleton,
	BodyNode* bn,
	const std::vector<BallInfo*>& ball_info,
				const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
				const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver)
	:mWorld(rigid_world),
	mSkeleton(skeleton),
	mAnchorPoint(std::make_pair(bn,Eigen::Vector3d(0,0,0))),
	mBallInfo(ball_info),
	mBallIndex(-1),mV(0),
	mIKOptimization(ik_optimization),
	mIKSolver(ik_solver),
	mTimeElapsed(0.0)
{

}
void
State::
AddEvent(const std::string& name,State* next_state)
{
	mEvents.insert(std::make_pair(name,next_state));	
}

State*
State::
GetNextState(const std::string& event_name)
{
	if(mEvents.find(event_name)==mEvents.end())
	{
		// std::cout<<"NO event name : "<<event_name<<std::endl;
		return nullptr;
	}
	return mEvents.at(event_name);
}
std::map<std::string,State*>&
State::
GetEvents()
{
	return mEvents;
}
void
State::
TimeStepping(double time_step)
{
	mTimeElapsed += time_step;
}


void
IKState::
Solve()
{
	bool need_ik_update= true;
	
	Eigen::Vector3d target;
	mBallInfo[mBallIndex]->ComputeFallingPosition(mBallInfo[mBallIndex]->releasedPoint[1],target);
	
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	
	const auto& ik_targets = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetTargets();
	for(auto& t : ik_targets)
	{
		if(!t.first.first->getName().compare(mAnchorPoint.first->getName()))
		{
			if((t.second-target).norm()<5E-3)
				need_ik_update = false;
		}
	}
	if(need_ik_update)
	{
		ik->AddTargetPositions(mAnchorPoint,target);
		mIKSolver->ReOptimizeTNLP(mIKOptimization);
	}
}
IKState::
IKState(const WorldPtr& rigid_world,
	const SkeletonPtr& skeleton,
	BodyNode* bn,
			const std::vector<BallInfo*>& ball_info,
			const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
			const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver)
	:State(rigid_world,skeleton,bn,ball_info,ik_optimization,ik_solver)
{

}

void 
IKState::
Initialize(int ball_index,int V)
{
	
	mBallIndex = ball_index;
	std::cout<<"Catch "<<mBallIndex<<std::endl;
	mV = V;
	Solve();
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	mTimeElapsed=0.0;
}
std::string 
IKState::
GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v)
{
	std::string event("no_event");
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));

	Solve();
	Eigen::VectorXd sol = ik->GetSolution();
	p = sol;
	v.resize(p.rows());
	v.setZero();

	Eigen::Vector3d body_COM = mAnchorPoint.first->getCOM();
	Eigen::Vector3d ball_COM;
	mBallInfo[mBallIndex]->GetPosition(ball_COM);
	if((body_COM-ball_COM).norm()<5E-2){
		std::cout<<"Attach "<<mBallIndex<<std::endl;
		mBallInfo[mBallIndex]->Attach(mWorld,mAnchorPoint.first);
		event = "catch";
	}

	return event;
}

Eigen::Vector2d
IKState::
GetTarget()
{
	const auto& ik_targets = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->GetTargets();
	for(auto& t : ik_targets)
	{
		if(!t.first.first->getName().compare(mAnchorPoint.first->getName()))
		{
			Eigen::Vector2d ret(t.second[0],t.second[1]);
			return ret;
		}
	}
	return Eigen::Vector2d(0,0);
}
BezierCurveState::
BezierCurveState(const WorldPtr& rigid_world,
	const SkeletonPtr& skeleton,
	BodyNode* bn,
	const std::vector<BallInfo*>& ball_info,
			const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
			const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver,
			double D,double T,
			int num_curve_sample)
	:State(rigid_world,skeleton,bn,ball_info,ik_optimization,ik_solver),mCurve(new BezierCurve())
	,mD(D),mT(T),mNumCurveSample(num_curve_sample)
{
}

void 
BezierCurveState::
Initialize(int ball_index,int V)
{
	mBallIndex = ball_index;
	std::cout<<"Swing "<<mBallIndex<<std::endl;
	mV = V;
	mMotions.clear();

	double t = ((double)mV-2*mD)*mT;
	Eigen::Vector3d p_3d;
	Eigen::Vector2d p0,p1,p2,v2;

	mBallInfo[mBallIndex]->GetPosition(p_3d);
	bool isleft = false;
	if(mAnchorPoint.first->getName().find("L")!=std::string::npos)
		isleft = true;

	
	// p0[0] = mAnchorPoint.first->getCOM()[0];
	// p0[1] = mAnchorPoint.first->getCOM()[1];
	// p0[1] = 0.2;
	// if(isleft)
	// 	p0[0] = -0.4;
	// else
	// 	p0[0] = 0.4;

	
	p0[1] = 0.2;
	if(isleft)
		p2[0] = -0.2;
	else
		p2[0] = 0.2;
	p2[1] = 0.2;
	v2[0] = -2.0*p2[0]/t;
	// v2[0] = 0;
	v2[1] = 0.5*9.81*t;
	p0[0] = p2[0] - mD*mT*v2[0];
	p1[0] = p2[0] - 0.5*mD*mT*v2[0];
	p1[1] = p2[1] - 0.5*mD*mT*v2[1];
	mTargetVelocity.setZero();
	mTargetVelocity[0] = v2[0];
	mTargetVelocity[1] = v2[1];
	std::cout<<"Target Velocity : "<<v2.transpose()<<std::endl;
	mCurve->Initialize(p0,p1,p2,mD*mT);

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	auto save_target = ik->GetTargets();
	for(int i =0;i<mNumCurveSample+2;i++)
	{
		double tt = ((double)i)/((double)mNumCurveSample) * mD*mT;
		
		Eigen::Vector2d p_ee = mCurve->GetPosition(tt);
		Eigen::Vector3d p_ee_3d(p_ee[0],p_ee[1],0);

		ik->AddTargetPositions(mAnchorPoint,p_ee_3d);
		mIKSolver->ReOptimizeTNLP(mIKOptimization);
		Eigen::VectorXd sol = ik->GetSolution();
		mMotions.push_back(std::make_pair(sol,tt));
	}
	for(auto& target : save_target)
		ik->AddTargetPositions(target.first,target.second);
	mTimeElapsed=0.0;
}
std::string 
BezierCurveState::
GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v)
{
	int k =0,k1 =0;
	for(int i =0;i<mMotions.size();i++)
	{
		if(mMotions[i].second<mTimeElapsed)
			k=i;
	}

	// if(k>=mMotions.size()-2)
	// {
	// 	mBallInfo[mBallIndex]->Release(mWorld);
	// }
	Eigen::Vector3d vv;
	mBallInfo[mBallIndex]->GetVelocity(vv);
	
	std::cout<<"Current Velocity : "<<vv.transpose()<<std::endl;
	if(k==mMotions.size()-1)// || ((mTargetVelocity-vv).norm()<5E-1)&&k>5){
	{
		std::cout<<"Release "<<mBallIndex<<std::endl;
		mBallInfo[mBallIndex]->Release(mWorld);
		// mBallInfo[mBallIndex]->skeleton->resetAccelerations();
		// mBallInfo[mBallIndex]->skeleton->setVelocities(Eigen::compose(Eigen::Vector3d(0,0,0),mTargetVelocity));
		// std::cout<<mBallInfo[mBallIndex]->skeleton->getBodyNode(0)->getCOMLinearVelocity()<<std::endl;
		return std::string("end");
	}

	k1 = k+1;
	double t = mTimeElapsed-mMotions[k].second;
	double dt = mMotions[k1].second-mMotions[k].second;

	t/= dt;
	p = (1.0-t)*(mMotions[k].first) + (t)*(mMotions[k1].first);

	double t1 = t+0.01;
	auto pdp =(1.0-t)*(mMotions[k].first) + (t)*(mMotions[k1].first);

	v = (pdp-p)/0.01;

	return std::string("no_event");
}




// void ReadMotion(const std::string path,std::vector<std::pair<double,Eigen::VectorXd>>& motion)
// {
// 	std::ifstream ifs(path);
// 	std::string str;
// 	std::stringstream ss;
// 	std::string t;
// 	std::string token;
// 	std::vector<std::string> token_vector;
// 	while(!ifs.eof())
// 	{
// 		str.clear();
// 		ss.clear();
// 		token_vector.clear();

// 		std::getline(ifs,str);
// 		ss.str(str);
// 		if(ss>>t){
// 			while(!ss.eof())
// 			{
// 				ss>>token;
// 				token_vector.push_back(token);
// 			}
// 			Eigen::VectorXd vec(token_vector.size());
// 			for(int i =0;i<token_vector.size();i++)
// 				vec[i] = std::stod(token_vector[i]);

// 			motion.push_back(std::make_pair(std::stod(t),vec));
// 		}
// 	}
// 	// std::cout<<path<<std::endl;
// 	// for(int i =0;i<motion.size();i++)
// 	// {
// 	// 	std::cout<<"time : "<<motion[i].first<<std::endl<<motion[i].second.transpose()<<std::endl;
// 	// }
// }


// MotionState::
// MotionState(const std::string& motion_path)
// {
// 	std::string path_to_motion_dir = "../vmcon2D/export/motions/";
	
// 	ReadMotion(path_to_motion_dir+motion_path,mMotion);
// }



// std::string
// MotionState::
// GetMotion(Eigen::VectorXd& motion)
// {
// 	int k =0,k1 =0;
// 	for(int i =0;i<mMotion.size();i++)
// 	{
// 		if(mMotion[i].first<mTimeElapsed)
// 			k=i;
// 	}

// 	if(k==mMotion.size()-1)
// 		return std::string("end");

// 	k1 = k+1;
// 	double t = mTimeElapsed-mMotion[k].first;
// 	double dt = mMotion[k1].first-mMotion[k].first;

// 	t/= dt;
// 	motion = (1.0-t)*(mMotion[k].second) + (t)*(mMotion[k1].second);

// 	return std::string("no_event");
// }

void
Machine::
InitializeJugglingState(const std::vector<int>& V)
{
	int max_V = 0;
	for(int i=0;i<V.size();i++)
		if(max_V<V[i])
			max_V = V[i];

	mJugglingStates.resize(V.size() + max_V);

	for(int i =0;i<V.size();i++)
	{
		if(mJugglingStates[i].ball_index==-1)
		{
			mJugglingStates[i].isLeftHand = i%2;
			mJugglingStates[i].ball_index = i;
			mJugglingStates[i].V = V[i];

			mJugglingStates[i + V[i]].ball_index = i;
			mJugglingStates[i + V[i]].isLeftHand = !mJugglingStates[i].isLeftHand;
		}
		else
		{
			mJugglingStates[i].V = V[i];
			mJugglingStates[i + V[i]].ball_index = mJugglingStates[i].ball_index;
			mJugglingStates[i + V[i]].isLeftHand = !mJugglingStates[i].isLeftHand;
		}
	}

	for(auto& js : mJugglingStates)
	{
		std::cout<<(js.isLeftHand?"LEFT ":"RIGHT ");
		std::cout<<js.ball_index<<" ";
		std::cout<<js.V<<std::endl;
	}
	mJugglingFrame = 0;
}
Machine::
Machine(const WorldPtr& rigid_world,
	MusculoSkeletalSystem* musculo_skeletal_system,
	const std::vector<BallInfo*>& balls,const double& dt)
	:mCurrentState(nullptr),mRigidWorld(rigid_world),mMusculoSkeletalSystem(musculo_skeletal_system),mBallInfo(balls),mdt(dt)
{
	mIKOptimization = new IKOptimization(mMusculoSkeletalSystem->GetSkeleton());

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

	Eigen::Vector3d l_loc = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL")->getCOM();
	Eigen::Vector3d r_loc = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM();

	ik->AddTargetPositions(std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL"),Eigen::Vector3d::Zero()),l_loc);
	ik->AddTargetPositions(std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d::Zero()),r_loc);

	mIKSolver->OptimizeTNLP(mIKOptimization);


	// std::vector<int> V_list{5,3,1,5,3,1,5,3,1,5,3,1,5,3,1};
	std::vector<int> V_list{3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
	// std::vector<int> V_list{5,5,5,5,5,5,5,5,5};
	// std::vector<int> V_list{7,7,7,7,7,7,7,7,7};
	InitializeJugglingState(V_list);
}


State*
Machine::
AddState(const std::string& name)
{
	bool is_left =false;
	bool is_catch = false;
	if(name.find("LEFT")!=std::string::npos)
		is_left = true;
	if(name.find("CATCH")!=std::string::npos)
		is_catch = true;
	if(is_catch)
	mStates.insert(std::make_pair(name,
		new IKState(
			mRigidWorld,
			mMusculoSkeletalSystem->GetSkeleton(),
			mMusculoSkeletalSystem->GetSkeleton()->getBodyNode((is_left?"HandL":"HandR")),
			mBallInfo,
			mIKOptimization,
			mIKSolver)));
	else
	mStates.insert(std::make_pair(name,
		new BezierCurveState(
			mRigidWorld,
			mMusculoSkeletalSystem->GetSkeleton(),
			mMusculoSkeletalSystem->GetSkeleton()->getBodyNode((is_left?"HandL":"HandR")),
			mBallInfo,
			mIKOptimization,
			mIKSolver)));
	if(mStates.size()==1)
		mCurrentState = mStates.at(name);
	
	return mStates.at(name);
}
void 	
Machine::
AddEvent(State* state_from,State* state_to,const std::string& name)
{
	state_from->AddEvent(name,state_to);
}

void
Machine::
SetCurrentState(State* s)
{
	mCurrentState = s;
}
State* 	
Machine::
GetCurrentState()
{
	return mCurrentState;
}
std::map<std::string,State*>&
Machine::
GetStates()
{
	return mStates;
}
void	
Machine::
Trigger(const std::string& name)
{
	std::string re_name = name;
	if(!name.compare("end"))
	{
		bool isPreviousLeft = mJugglingStates[mJugglingFrame].isLeftHand;
		std::cout<<"isPreviousLeft : "<<isPreviousLeft<<std::endl;
		for(int i = mJugglingFrame+1;i<mJugglingStates.size();i++)
		{
			std::cout<<"\nFrame : "<<i<<std::endl;
			std::cout<<"isLeft : "<<mJugglingStates[i].isLeftHand<<std::endl;
			std::cout<<"Ball index "<<mJugglingStates[i].ball_index<<std::endl;
			if(isPreviousLeft == mJugglingStates[i].isLeftHand)
			{
				if(mBallInfo[mJugglingStates[i].ball_index]->isReleased)
				{
					
					IKState* ikstate = dynamic_cast<IKState*>(mCurrentState->GetNextState(name +"_same_hand"));
					ikstate->Initialize(mJugglingStates[i].ball_index,mJugglingStates[i].V);
					
				}
				break;
			}
		}

		mJugglingFrame++;
		bool isCurrentLeft = mJugglingStates[mJugglingFrame].isLeftHand;
		if(isCurrentLeft == isPreviousLeft)
			re_name = name +"_same_hand";

	}
	State* next_state;
	next_state = mCurrentState->GetNextState(re_name);
	if(next_state!=nullptr){
		mCurrentState = next_state;
		mCurrentState->Initialize(mJugglingStates[mJugglingFrame].ball_index,mJugglingStates[mJugglingFrame].V);
	}	
	
}

void
Machine::
GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v)
{
	mCurrentState->TimeStepping(mdt);
	std::string event = mCurrentState->GetMotion(p,v);	
	Trigger(event);
}