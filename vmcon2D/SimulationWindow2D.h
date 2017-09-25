#ifndef __SIMULATION_WINDOW_H__
#define __SIMULATION_WINDOW_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "GUI/Window2D.h"
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <vector>
class Controller;
class MuscleOptimization;
class MusculoSkeletalSystem;
class Machine;
class State;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
enum MOUSE_MODE
{
	CAMERA_CONTROL,
	CONSTRAINT_CONTROL
};
struct Record
{
	double			time;
	std::vector<Eigen::VectorXd> rigid_body_positions;
	Eigen::VectorXd soft_body_positions;
	Eigen::VectorXd activation_levels;
	std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> muscle_forces;
	State*	state;
};
class SimulationWindow2D : public Window2D
{
protected:
	MOUSE_MODE					mMouseMode;
	FEM::AttachmentConstraint* 	mDragConstraint;
	AnchorPoint					mDragAnchorPoint;
	
	FEM::World*					mSoftWorld;
	dart::simulation::WorldPtr  mRigidWorld;
	MusculoSkeletalSystem*		mMusculoSkeletalSystem;

	std::vector<dart::dynamics::SkeletonPtr> mBalls;

	Controller*					mController;
	bool 						mIsPlay;
	bool 						mIsReplay;
	bool 						mIsPaused;
	
	std::vector<Record*>		mRecords;
	int 						mRecordFrame;
	double						mTime,mTimeStep,mSimTime;
public:
	SimulationWindow2D();
	void Initialize();

	bool TimeStepping();  //return true if soft simulation is updated.
	void SetRecord(Record* rec);
public:
	
	void Display() override;
	void Keyboard(unsigned char key,int x,int y) override;
	void Mouse(int button, int state, int x, int y) override;
	void Motion(int x, int y) override;
	void Reshape(int w, int h) override;
	void Timer(int value) override;
};

#endif