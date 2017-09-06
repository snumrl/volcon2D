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

class MuscleOptimization;
class MusculoSkeletalSystem;
enum MOUSE_MODE
{
	CAMERA_CONTROL,
	CONSTRAINT_CONTROL
};
struct Record
{
	double			time;
	Eigen::VectorXd rigid_body_positions;
	Eigen::VectorXd soft_body_positions;
	Eigen::VectorXd activation_levels;
	std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> muscle_forces;

};

class SimulationWindow2D : public Window2D
{
protected:
	MOUSE_MODE					mMouseMode;
	FEM::AttachmentConstraint* 	mDragConstraint;
	FEM::World*					mSoftWorld;
	dart::simulation::WorldPtr  mRigidWorld;
	MusculoSkeletalSystem*		mMusculoSkeletalSystem;

	Eigen::VectorXd							mRestPose;
	Eigen::VectorXd							mPreviousPose;
	Ipopt::SmartPtr<Ipopt::TNLP> 			 mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> mMuscleOptimizationSolver;

	bool 						mIsPlay;
	bool 						mIsReplay;
	
	std::vector<Record*>		mRecords;
	int 						mRecordFrame;
	double						mTime,mTimeStep,mSimTime;
public:
	SimulationWindow2D();
	void Initialize();

	void TimeStepping();
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