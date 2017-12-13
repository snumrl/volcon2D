#ifndef __SIMULATION_WINDOW_H__
#define __SIMULATION_WINDOW_H__
#include "GUI/Window2D.h"
#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"

class Object;
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
	FEM::World*					mSoftWorld;
	std::vector<Object*> 		mObjects;
	

	bool 						mIsPlay;

public:
	SimulationWindow2D();
	void Initialize();
	void TimeStepping();

	void Plotting(const Eigen::VectorXd& a, Eigen::VectorXd& f0,Eigen::VectorXd& f1);
public:
	void Display() override;
	void Keyboard(unsigned char key,int x,int y) override;
	void Mouse(int button, int state, int x, int y) override;
	void Motion(int x, int y) override;
	void Reshape(int w, int h) override;
	void Timer(int value) override;
};



#endif