#include "SimulationWindow2D.h"
#include "Object.h"
#include "GUI/Camera2D.h"
#include "GUI/GL_function.h"
#include "FEM2D_Interface.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GL/glut.h"
#include <fstream>

using namespace FEM;

SimulationWindow2D::
SimulationWindow2D()
	:mMouseMode(MOUSE_MODE::CAMERA_CONTROL),
	mDragConstraint(new AttachmentConstraint(50000.0,0,Eigen::Vector2d(0,0))),
	mIsPlay(false)
{
	mObjects.push_back(new Object(Eigen::Vector2d(0,0)));
	mObjects.push_back(new Object(Eigen::Vector2d(0,0.5)));
	Initialize();
	mDisplayTimeout = mSoftWorld->GetTimeStep()*1000;
}

void
SimulationWindow2D::
Initialize()
{
	mSoftWorld = new FEM::World(
		// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/200.0,										//time_step
		100, 											//max_iteration	
		0.999											//damping_coeff
		);

	// MakeMuscles("../fem_test/export/muscle.xml",mObject);
	for(auto& obj : mObjects)
		obj->Initialize(mSoftWorld);
	mSoftWorld->Initialize();
}

void
SimulationWindow2D::
Plotting(const Eigen::VectorXd& a, Eigen::VectorXd& f0,Eigen::VectorXd& f1)
{
	Eigen::VectorXd act(1);
	act[0] = 1.0;
	for(int i =0;i<a.rows();i++)
	{
		auto prev_p0 = mObjects[0]->mAttachementConstraintVector[0]->GetP();
		prev_p0[0] += -0.01;
		mObjects[0]->mAttachementConstraintVector[0]->GetP() = prev_p0;
		mObjects[0]->SetActivationLevel(act);
		mSoftWorld->TimeStepping();	
		f0[i] = mObjects[0]->ComputeForce(mSoftWorld);
	}

	// for(int i =0;i<a.rows();i++)
	// {
	// 	Eigen::VectorXd act(1);
	// 	act[0] = a[i];
	// 	mObjects[0]->SetActivationLevel(act);
	// 	mObjects[1]->SetActivationLevel(act);
	// 	f0[i] = mObjects[0]->ComputeForce(mSoftWorld);
	// 	f1[i] = mObjects[1]->ComputeForce(mSoftWorld);
	// 	mSoftWorld->TimeStepping();	
	// }
	// std::cout<<a.transpose()<<std::endl;
	std::cout<<f0.transpose()<<std::endl;
	// std::cout<<f1.transpose()<<std::endl;
	
}
void
SimulationWindow2D::
Display()
{
	glClearColor(0.95, 0.95, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_DEPTH_TEST);
	mCamera->Apply();
	
	glColor3f(0,0,0);
	glLineWidth(1.0);
	glBegin(GL_LINES);
	for(double x=-10.0;x<=10.0;x+=0.5)
	{
		glVertex2f(x,-10.0);
		glVertex2f(x,10.0);
	}
	for(double y=-10.0;y<=10.0;y+=0.5)
	{
		glVertex2f(-10.0,y);
		glVertex2f(10.0,y);
	}
	glEnd();
    glColor3f(0,0,0);

	glColor3f(0,0,1);
	glLineWidth(5.0);
	glBegin(GL_LINES);
	glVertex2f(-40.0,-8);
	glVertex2f(40.0,-8);
	glEnd();

	glColor3f(0,0,0);
    DrawStringOnScreen(0.8,0.2,std::to_string(mSoftWorld->GetTime()),true);
	glLineWidth(1.0);

	const auto& x = mSoftWorld->GetPositions();
	
	for(auto& c : mSoftWorld->GetConstraints())
	{
		DrawConstraint(c,x);
	}
	if(mIsDrag)
		DrawConstraint(mDragConstraint,x);
	glEnable(GL_DEPTH_TEST);
	glutSwapBuffers();
}
void
SimulationWindow2D::
Keyboard(unsigned char key,int x,int y)
{
	Eigen::VectorXd act(1);
	if (key == 'a')
	{
		Eigen::Vector2d prev_p0,prev_p1;
		prev_p0 = mObjects[0]->mAttachementConstraintVector[0]->GetP();
		prev_p1 = mObjects[1]->mAttachementConstraintVector[0]->GetP();
		prev_p0[0] += -0.05;
		prev_p1[0] += -0.05;
		mObjects[0]->mAttachementConstraintVector[0]->GetP() = prev_p0;
		mObjects[1]->mAttachementConstraintVector[0]->GetP() = prev_p1;
	}
	else if(key == 'd')
	{
		Eigen::Vector2d prev_p0,prev_p1;
		prev_p0 = mObjects[0]->mAttachementConstraintVector[0]->GetP();
		prev_p1 = mObjects[1]->mAttachementConstraintVector[0]->GetP();
		prev_p0[0] += 0.05;
		prev_p1[0] += 0.05;
		mObjects[0]->mAttachementConstraintVector[0]->GetP() = prev_p0;
		mObjects[1]->mAttachementConstraintVector[0]->GetP() = prev_p1;
	}
	Eigen::VectorXd act_set = Eigen::VectorXd::LinSpaced(20, 0, 1.0);
	Eigen::VectorXd f0(act_set.rows()),f1(act_set.rows());
	switch(key)
	{
		case 't' : 
		if(mMouseMode== CAMERA_CONTROL){
			mMouseMode = CONSTRAINT_CONTROL;
			std::cout<<"CONSTRAINT CONTROL mode"<<std::endl;
		}
		else
		{
			mMouseMode = CAMERA_CONTROL;
			std::cout<<"CAMERA_CONTROL mode"<<std::endl;
		}
		break;

		case 'z' : Plotting(act_set,f0,f1);break;
		case ' ' : mIsPlay = !mIsPlay; break;	
		case 'r' : 
		act = mObjects[0]->GetActivationLevel();
		act.setZero();
		mObjects[0]->SetActivationLevel(act);
		mObjects[1]->SetActivationLevel(act);
		break;break;
		case '1' :
		act = mObjects[0]->GetActivationLevel();
		act[0] +=0.05;
		if(act[0]>1.0)
			act[0]=1.0;
		mObjects[0]->SetActivationLevel(act);
		// break;
		case '2' :
		act = mObjects[1]->GetActivationLevel();
		act[0] +=0.05;
		if(act[0]>1.0)
			act[0]=1.0;
		mObjects[1]->SetActivationLevel(act);
		break;
		
		case 27: exit(0);break;
		default : break;
	}

	glutPostRedisplay();
}

void
SimulationWindow2D::
Mouse(int button, int state, int x, int y)
{
	
	if (state == GLUT_DOWN)
	{	
		auto mouse_world = mCamera->GetWorldPosition(x,y);
		mIsDrag = true;
		mMouseType = button;
		mPrevX = x;
		mPrevY = y;


		if(mMouseMode==CONSTRAINT_CONTROL)
		{
		}
	}
	else
	{
		if(mMouseMode==CONSTRAINT_CONTROL) {
		}

		mIsDrag = false;
		mMouseType = 0;
	}

	glutPostRedisplay();
}
void
SimulationWindow2D::
Motion(int x, int y)
{
	if (!mIsDrag)
		return;

	if(mMouseMode==CAMERA_CONTROL)
	{
		if (mMouseType == GLUT_LEFT_BUTTON)
			mCamera->Translate(x,y,mPrevX,mPrevY);
		else if (mMouseType == GLUT_RIGHT_BUTTON)
			mCamera->Pan(x,y,mPrevX,mPrevY);
	}
	else if(mMouseMode==CONSTRAINT_CONTROL)
	{
		auto mouse_world = mCamera->GetWorldPosition(x,y);
	}
	mPrevX = x;
	mPrevY = y;
	glutPostRedisplay();
}
void
SimulationWindow2D::
Reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	mCamera->Apply();
	glutPostRedisplay();
}

void
SimulationWindow2D::
Timer(int value)
{
	if(mIsPlay)
		mSoftWorld->TimeStepping();
	glutPostRedisplay();
	glutTimerFunc(mDisplayTimeout, TimerEvent,1);
}

