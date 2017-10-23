#include "FSM.h"
#include "Record.h"
#include "BezierCurve.h"
#include "IKOptimization.h"
#include "MusculoSkeletalSystem.h"
#include "Controller.h"
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
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
State(const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
				BodyNode* bn,
				const std::vector<BallInfo*>& ball_info,
				const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
				const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver)
	:mRigidWorld(rigid_world),
	mSoftWorld(soft_world),
	mMusculoSkeletalSystem(musculo_skeletal_system),
	mController(controller),
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
IKState(const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
	BodyNode* bn,
			const std::vector<BallInfo*>& ball_info,
			const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
			const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver)
	:State(rigid_world,soft_world,musculo_skeletal_system,controller,bn,ball_info,ik_optimization,ik_solver)
{

}

void 
IKState::
Initialize(int ball_index,int V)
{
	
	mBallIndex = ball_index;
	// std::cout<<"Catch "<<mBallIndex<<std::endl;
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
		// std::cout<<"Attach "<<mBallIndex<<std::endl;
		mBallInfo[mBallIndex]->Attach(mRigidWorld,mAnchorPoint.first);
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

void
BezierCurveState::
GenerateMotions(BezierCurve* bc,std::vector<std::pair<Eigen::VectorXd,double>>& motions)
{
	motions.clear();

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	Eigen::VectorXd save_positions = ik->GetSolution();

	auto save_target = ik->GetTargets();
	for(int i =0;i<mNumCurveSample+4;i++)
	{
		double tt = ((double)i)/((double)mNumCurveSample) * mD*mT;
		
		Eigen::Vector2d p_ee = bc->GetPosition(tt);
		Eigen::Vector3d p_ee_3d(p_ee[0],p_ee[1],0);

		ik->AddTargetPositions(mAnchorPoint,p_ee_3d);
		mIKSolver->ReOptimizeTNLP(mIKOptimization);
		Eigen::VectorXd sol = ik->GetSolution();
		motions.push_back(std::make_pair(sol,tt));
	}

	for(auto& target : save_target)
		ik->AddTargetPositions(target.first,target.second);

	ik->SetSolution(save_positions);
}
void
BezierCurveState::
OptimizeBezierCurvePoint(int num_samples)
{
	double t = ((double)mV-2*mD+0.3)*mT;

	Eigen::Vector3d p_3d;
	Eigen::Vector2d p0,p1,p2,v2;

	mBallInfo[mBallIndex]->GetPosition(p_3d);
	bool isleft = false;
	if(mAnchorPoint.first->getName().find("L")!=std::string::npos)
		isleft = true;

	if(isleft)
		p2[0] = -0.3;
	else
		p2[0] = 0.3;
	p2[1] = 0.2;

	v2[0] = -2.0*p2[0]/t;
	// v2[0] = 0;
	v2[1] = 0.5*9.81*t;

	std::cout<<"Target Velocity : "<<v2.transpose()<<std::endl;

	p0[1] = 0.2;
	//Free parameter Optimization
	p0[0] = p2[0] - mD*mT*v2[0];
	p1[0] = p2[0] - 0.5*mD*mT*v2[0];
	p1[1] = p2[1] - 0.5*mD*mT*v2[1];

	std::cout<<"LOOK AHEAD start"<<std::endl;

	double vel_diff_norm = 1E6;
	int sample_check = -1;
	int result_release_count = -1;
	Eigen::Vector2d result_velocity;
	Eigen::Vector2d rp0,rp1;

	for(int j=0;j<3;j++)
	{
		std::vector<std::pair<Eigen::Vector4d,Eigen::Vector3d>> samples;
		for(int i =0;i<num_samples;i++)
		{
			std::cout<<"\n\tSample "<<i<<std::endl;
			Eigen::Vector2d new_p0,new_p1;
			new_p0 = p0;
			new_p1 = p1;
			double deltax0 = dart::math::random(-0.03,0.03);
			double deltay0 = dart::math::random(-0.07,0.07);
			double deltax1 = dart::math::random(-0.1,0.1);
			double deltay1 = dart::math::random(-0.15,0.15);

			for(int k=0;k<j;k++)
			{
				deltax0 *=0.5;
				deltay0 *=0.5;
				deltax1 *=0.5;
				deltay1 *=0.5;
			}

			new_p0[0] += deltax0;
			new_p0[1] += deltay0;
			new_p1[0] += deltax1;
			new_p1[1] += deltay1;

			Eigen::Vector2d v_result;
			int release_count = -1;
			
			std::cout<<"\tp : "<<new_p0.transpose()<<"\t"<<new_p1.transpose()<<std::endl;
			GenerateSample(new_p0,new_p1,p2,v2,v_result,release_count);

			std::pair<Eigen::Vector4d,Eigen::Vector3d> sample_pair;
			sample_pair.first.block<2,1>(0,0) = new_p0;
			sample_pair.first.block<2,1>(2,0) = new_p1;
			sample_pair.second.block<2,1>(0,0) = v_result;
			sample_pair.second[2] = release_count;
			samples.push_back(sample_pair);
			if((v2-v_result).norm()<vel_diff_norm)
			{
				vel_diff_norm = (v2-v_result).norm();
				
				result_velocity = v_result;
				result_release_count = release_count;
				rp0 = new_p0;
				rp1 = new_p1;
				sample_check = i;
			}
		}

		//Linear Regression
		for(int i =0;i<num_samples;i++)
			std::cout<<samples[i].first.transpose()<<"\t"<<samples[i].second.transpose()<<std::endl;

		Eigen::MatrixXd xxt(4,4),yxt(3,4);
		Eigen::VectorXd x_sum(4),y_sum(3);
		xxt.setZero();
		yxt.setZero();
		x_sum.setZero();
		y_sum.setZero();
		for(int i =0;i<num_samples;i++)
		{
			xxt += samples[i].first*samples[i].first.transpose();
			yxt += samples[i].second*samples[i].first.transpose();
			x_sum += samples[i].first;
			y_sum += samples[i].second;
		}
		Eigen::MatrixXd J(3,4);
		Eigen::VectorXd b(3);
		J = yxt*xxt.inverse();
		b = y_sum-J*x_sum;

		for(int i =0;i<num_samples;i++)
			std::cout<<(J*samples[i].first+b).transpose()<<std::endl;
		

		Eigen::Vector3d reg_target;
		Eigen::Vector4d reg_p;
		reg_target.block<2,1>(0,0) = v2;
		reg_target[2] = y_sum[2]/(double)num_samples;

		reg_p = J.transpose()*((J*J.transpose()).inverse())*(reg_target-b);
		std::cout<<"regresion result : "<<reg_p.transpose()<<std::endl;
		// p0 = rp0;
		// p1 = rp1;
		p0 = (0.3*reg_p.block<2,1>(0,0)+0.7*rp0);
		p1 = (0.3*reg_p.block<2,1>(2,0)+0.7*rp1);

		std::cout<<"\n\tNew initial Point : : "<<p0.transpose()<<"\t"<<p1.transpose()<<std::endl;
		std::cout<<"Sample "<<sample_check<<" is choosed."<<std::endl;
		std::cout<<"Residual : "<<vel_diff_norm<<std::endl;
		std::cout<<"velocity : "<<result_velocity.transpose()<<std::endl;
		std::cout<<std::endl;
	}

	
	std::cout<<"LOOK AHEAD end"<<std::endl;
	std::cout<<"Sample "<<sample_check<<" is choosed."<<std::endl;
	
	mReleaseCount = result_release_count;
	mCurve->Initialize(rp0,rp1,p2,mD*mT);
	
	std::cout<<"Residual : "<<vel_diff_norm<<std::endl;
	std::cout<<"velocity : "<<result_velocity.transpose()<<std::endl;
	std::cout<<"target velocity : "<<v2.transpose()<<std::endl;
	std::cout<<"Release count : "<<mReleaseCount<<std::endl;
	
}
void
BezierCurveState::
GenerateSample(const Eigen::Vector2d& p0,const Eigen::Vector2d& p1,const Eigen::Vector2d& p2, const Eigen::Vector2d& v_target,
				Eigen::Vector2d& v_result,int& release_count)
{
	mCount = 0;
	mReleaseCount = 10000;
	mTimeElapsed = 0.0;

	mCurve->Initialize(p0,p1,p2,mD*mT);
	GenerateMotions(mCurve,mMotions);
	Record*	current_record = new Record();
	current_record->Set(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);

	


	int closest_velocity_count = 0;
	double vel_diff_norm = 1E6;
	Eigen::Vector2d vel_closest;
	while(true)
	{
		//Simulation START
		bool is_fem_updated = false;
		if(mSoftWorld->GetTime()<=mRigidWorld->getTime())
		{
			is_fem_updated =true;
			mMusculoSkeletalSystem->SetActivationLevel(mController->Compute());	
		}

		mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);

		if(is_fem_updated)
		{
			mMusculoSkeletalSystem->TransformAttachmentPoints();
			mSoftWorld->TimeStepping();
		}
		mRigidWorld->step();
		//Simulation END

		if(is_fem_updated)
		{
			Eigen::Vector3d velocity;
			mBallInfo[mBallIndex]->GetVelocity(velocity);
			// std::cout<<"\t\t"<<mCount<<"\t"<<velocity.block<2,1>(0,0).transpose()<<std::endl;
			if( (v_target-velocity.block<2,1>(0,0)).norm()<vel_diff_norm)
			{
				vel_diff_norm = (v_target-velocity.block<2,1>(0,0)).norm();
				closest_velocity_count = mCount;
				vel_closest = velocity.block<2,1>(0,0);
				// std::cout<<mCount<<std::endl;
			}
		}
		if(mTimeElapsed>mMotions[mMotions.size()-2].second)
			break;
	}

	
	std::cout<<"\tClosest count : "<<closest_velocity_count<<std::endl;
	std::cout<<"\tres : "<<vel_diff_norm<<std::endl;
	std::cout<<"\tvelocity : "<<vel_closest.transpose()<<std::endl;
	

	current_record->Get(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);

	mTimeElapsed = 0.0;
	mCount = 0;

	release_count = closest_velocity_count;
	v_result = vel_closest;
}



BezierCurveState::
BezierCurveState(const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
	BodyNode* bn,
	const std::vector<BallInfo*>& ball_info,
			const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
			const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver,
			double D,double T,
			int num_curve_sample)
	:State(rigid_world,soft_world,musculo_skeletal_system,controller,bn,ball_info,ik_optimization,ik_solver),mCurve(new BezierCurve())
	,mD(D),mT(T),mNumCurveSample(num_curve_sample)
{
}

void 
BezierCurveState::
Initialize(int ball_index,int V)
{
	// mCount = 0;
	// mReleaseCount = 10000;
	mBallIndex = ball_index;
	mV = V;
	// std::cout<<"Swing "<<mBallIndex<<std::endl;

	OptimizeBezierCurvePoint(5);		//This optimize bezier mCurve control point.

	GenerateMotions(mCurve,mMotions);

	// mTimeElapsed=0.0;
	//Look Ahead step
	// Record*	current_record = new Record();
	// current_record->Set(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);
	// std::cout<<"LOOK AHEAD start"<<std::endl;


	// int closest_velocity_count = 0;
	// double vel_diff_norm = 1E6;
	// while(true)
	// {
	// 	//Simulation START
	// 	bool is_fem_updated = false;
	// 	if(mSoftWorld->GetTime()<=mRigidWorld->getTime())
	// 	{
	// 		is_fem_updated =true;
	// 		mMusculoSkeletalSystem->SetActivationLevel(mController->Compute());	
	// 	}

	// 	mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);

	// 	if(is_fem_updated)
	// 	{
	// 		mMusculoSkeletalSystem->TransformAttachmentPoints();
	// 		mSoftWorld->TimeStepping();
	// 	}
	// 	mRigidWorld->step();
	// 	//Simulation END

	// 	Eigen::Vector3d velocity;
	// 	mBallInfo[mBallIndex]->GetVelocity(velocity);
	// 	if( (mTargetVelocity-velocity).norm()<vel_diff_norm)
	// 	{
	// 		vel_diff_norm = (mTargetVelocity-velocity).norm();
	// 		closest_velocity_count = mCount;
	// 	}
	// 	if(mTimeElapsed>mMotions[mMotions.size()-2].second)
	// 		break;
	// }

	
	// std::cout<<"LOOK AHEAD end"<<std::endl;
	// std::cout<<"result : \n";
	// std::cout<<"Closest count : "<<closest_velocity_count<<std::endl;
	// std::cout<<"difference : "<<vel_diff_norm<<std::endl;



	// current_record->Get(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);

	// mReleaseCount = closest_velocity_count;
	// mTimeElapsed = 0.0;
	// mCount = 0;
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

	
/*	Eigen::Vector3d vv;
	mBallInfo[mBallIndex]->GetVelocity(vv);
	
	std::cout<<mCount<<" Current Velocity : "<<vv.transpose()<<std::endl;
	std::cout<<mCount<<"\t"<<mReleaseCount<<std::endl;*/
	if(mCount == mReleaseCount-1)// || ((mTargetVelocity-vv).norm()<5E-1)&&k>5){
	{
		// std::cout<<"Release "<<mBallIndex<<std::endl;
		mBallInfo[mBallIndex]->Release(mRigidWorld);

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
	mCount++;
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
	FEM::World* soft_world,
	MusculoSkeletalSystem* musculo_skeletal_system,
	const std::vector<BallInfo*>& balls,
	Controller* controller,const double& dt)
	:mCurrentState(nullptr),mRigidWorld(rigid_world),mSoftWorld(soft_world),mMusculoSkeletalSystem(musculo_skeletal_system),mBallInfo(balls),mController(controller),mdt(dt)
{
	dart::math::seedRand();
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
	// std::vector<int> V_list{3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
	// std::vector<int> V_list{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	// std::vector<int> V_list{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	std::vector<int> V_list{5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
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
			mSoftWorld,
			mMusculoSkeletalSystem,
			mController,
			mMusculoSkeletalSystem->GetSkeleton()->getBodyNode((is_left?"HandL":"HandR")),
			mBallInfo,
			mIKOptimization,
			mIKSolver)));
	else
	mStates.insert(std::make_pair(name,
		new BezierCurveState(
			mRigidWorld,
			mSoftWorld,
			mMusculoSkeletalSystem,
			mController,
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
		// std::cout<<"isPreviousLeft : "<<isPreviousLeft<<std::endl;
		for(int i = mJugglingFrame+1;i<mJugglingStates.size();i++)
		{
			// std::cout<<"\nFrame : "<<i<<std::endl;
			// std::cout<<"isLeft : "<<mJugglingStates[i].isLeftHand<<std::endl;
			// std::cout<<"Ball index "<<mJugglingStates[i].ball_index<<std::endl;
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
	if(dynamic_cast<BezierCurveState*>(mCurrentState)!=nullptr)
	{
		if(dynamic_cast<BezierCurveState*>(mCurrentState)->mReleaseCount != 10000)
		{


		bool isCurrentLeft = mJugglingStates[mJugglingFrame].isLeftHand;
		for(int i = mJugglingFrame+1;i<mJugglingStates.size();i++)
		{
			bool isLeft = mJugglingStates[i].isLeftHand;
			if(isCurrentLeft == !isLeft)
			{
				IKState* ikstate = dynamic_cast<IKState*>(mCurrentState->GetNextState("end"));

				Eigen::Vector3d body_COM = ikstate->mAnchorPoint.first->getCOM();
				Eigen::Vector3d ball_COM;
				mBallInfo[mJugglingStates[i].ball_index]->GetPosition(ball_COM);
				// std::cout<<"\t\t\t"<<ball_COM.transpose()<<std::endl;
				// std::cout<<"\t\t\t"<<body_COM.transpose()<<std::endl;
				if((body_COM-ball_COM).norm()<5E-2){
					mBallInfo[mJugglingStates[i].ball_index]->Attach(mRigidWorld,ikstate->mAnchorPoint.first);
				}
				break;
			}
		}
		}
	}
	mCurrentState->TimeStepping(mdt);
	std::string event = mCurrentState->GetMotion(p,v);	
	Trigger(event);

}