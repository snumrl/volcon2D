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
class VelocityControlDDP;
class Record;
class Controller;
class MuscleOptimization;
class MusculoSkeletalSystem;
class BallInfo;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
enum MOUSE_MODE
{
	CAMERA_CONTROL,
	CONSTRAINT_CONTROL
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

	FEM::World*					mDDPSoftWorld;
	dart::simulation::WorldPtr  mDDPRigidWorld;
	MusculoSkeletalSystem*		mDDPMusculoSkeletalSystem;
	std::vector<dart::dynamics::SkeletonPtr> mDDPBalls;

	Controller*					mController;
	VelocityControlDDP*			mDDP;
	std::vector<Eigen::VectorXd> mU;
	int u_index;
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