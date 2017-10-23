#ifndef __FSM_H__
#define __FSM_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <Eigen/Core>

class BezierCurve;
class IKOptimization;
class MusculoSkeletalSystem;
class Record;
class Controller;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
class BallInfo
{
public:

	dart::dynamics::SkeletonPtr 				skeleton;
	dart::constraint::WeldJointConstraintPtr 	constraint;

	bool 										isReleased;

	Eigen::Vector3d								releasedPoint;
	Eigen::Vector3d								releasedVelocity;
	Eigen::Vector6d								releasedImpulse;
	int 										releaseFrameCount;
public:
	BallInfo(const dart::constraint::WeldJointConstraintPtr& cons,const dart::dynamics::SkeletonPtr& skel);

	void ComputeFallingPosition(double h,Eigen::Vector3d& fp);
	void Release(const dart::simulation::WorldPtr& world);
	void GetPosition(Eigen::Vector3d& p);
	void GetVelocity(Eigen::Vector3d& v);
	void Attach(const dart::simulation::WorldPtr& world,dart::dynamics::BodyNode* bn);
	void TimeStepping();
};
struct JugglingState
{
	bool isLeftHand;  		// true : left hand, false : right hand
	int ball_index;
	int V;					//SiteSwap Value;

	JugglingState()
		:isLeftHand(false),ball_index(-1),V(0){};
};

class State
{
public:
	
	dart::simulation::WorldPtr							mRigidWorld;
	FEM::World*											mSoftWorld;
	MusculoSkeletalSystem*								mMusculoSkeletalSystem;
	Controller*											mController;

	AnchorPoint 										mAnchorPoint;

	std::vector<BallInfo*>								mBallInfo;
	int 												mBallIndex;
	int 												mV;

	Ipopt::SmartPtr<Ipopt::TNLP>			 			mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 			mIKSolver;

	std::map<std::string,State*>						mEvents;
	double												mTimeElapsed;

public:
	State(		const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
				dart::dynamics::BodyNode* bn,
				const std::vector<BallInfo*>& ball_info,
				const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
				const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver);

	void 	AddEvent(const std::string& name,State* next_state);
	State*	GetNextState(const std::string& event_name);
	std::map<std::string,State*>& GetEvents();

	
	void TimeStepping(double time_step);
	virtual void Initialize(int ball_index,int V) = 0; 
	double GetTime() {return mTimeElapsed;};
	virtual std::string GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v) = 0;
};

class IKState : public State
{
protected:
	void Solve();
public:
	IKState(	const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
				dart::dynamics::BodyNode* bn,
				const std::vector<BallInfo*>& ball_info,
				const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
				const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver
				);	

	void Initialize(int ball_index,int V) override;
	std::string GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v) override;
	Eigen::Vector2d GetTarget();
};

class BezierCurveState : public State
{
public:

	BezierCurve* mCurve;
	// Eigen::Vector3d								   mTargetVelocity;
	std::vector<std::pair<Eigen::VectorXd,double>> mMotions;
	double mD,mT;
	int mNumCurveSample;
	int mReleaseCount;
	int mCount;

	
	void GenerateMotions(BezierCurve* bc,std::vector<std::pair<Eigen::VectorXd,double>>& motions);
	void OptimizeBezierCurvePoint(int num_samples);
	void GenerateSample(
		const Eigen::Vector2d& p0,const Eigen::Vector2d& p1,const Eigen::Vector2d& p2, const Eigen::Vector2d& v_target,
		Eigen::Vector2d& v_result,int& release_count);
public:
	BezierCurveState(const dart::simulation::WorldPtr& rigid_world,
				FEM::World* soft_world,
				MusculoSkeletalSystem* musculo_skeletal_system,
				Controller* controller,
				dart::dynamics::BodyNode* bn,
				const std::vector<BallInfo*>& ball_info,
				const Ipopt::SmartPtr<Ipopt::TNLP>& ik_optimization,
				const Ipopt::SmartPtr<Ipopt::IpoptApplication>& ik_solver,
				double D = 0.7,
				double T = 0.3,
				int num_curve_sample = 10
				);

	void Initialize(int ball_index,int V) override;
	std::string GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v) override;
	BezierCurve* GetCurve() {return mCurve;};
};




class Machine
{
private:
	State*										mCurrentState;

	std::map<std::string,State*>				mStates;

	dart::simulation::WorldPtr  				mRigidWorld;
	FEM::World*									mSoftWorld;
	MusculoSkeletalSystem*						mMusculoSkeletalSystem;
	
	std::vector<BallInfo*>						mBallInfo;
	Controller*									mController;

	Ipopt::SmartPtr<Ipopt::TNLP>			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;

	std::vector<JugglingState>					mJugglingStates;
	int 										mJugglingFrame;

	double										mdt;
	void InitializeJugglingState(const std::vector<int>& V);
public:
	Machine(
		const dart::simulation::WorldPtr& rigid_world,
		FEM::World* soft_world,
		MusculoSkeletalSystem* musculo_skeletal_system,
		const std::vector<BallInfo*>& balls,
		Controller* controller,const double& dt);
	
	State* 	AddState(const std::string& name);
	void 	AddEvent(State* state_from,State* state_to,const std::string& name);

	void 	SetCurrentState(State* s);
	State* 	GetCurrentState();
	std::map<std::string,State*>& GetStates();
	void	Trigger(const std::string& event_name);
	void GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v);

};



#endif