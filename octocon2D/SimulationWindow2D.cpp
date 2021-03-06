#include "SimulationWindow2D.h"
#include "Octopus.h"
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
	mOctopus(new Octopus()),
	mIsPlay(false)
{
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
		// FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/500.0,										//time_step
		50, 											//max_iteration	
		0.999											//damping_coeff
		);

	MakeMuscles("../octocon2D/export/muscle.xml",mOctopus);
	mOctopus->Initialize(mSoftWorld);
	mSoftWorld->Initialize();
}

void
SimulationWindow2D::
PrintPosition() {
	std::cout << "Print Target" << std::endl;
	auto target = mSoftWorld->GetPositions();
	
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
	auto act = mOctopus->GetActivationLevel();
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
	
		case ' ' : mIsPlay = !mIsPlay; break;	
		case 'r' : act.setZero();break;
		case '1' : act[0] += 0.05;break;
		case '2' : act[1] += 0.05;break;
		case '3' : act[2] += 0.05;break;
		case '4' : act[3] += 0.05;break;
		case '5' : act[4] += 0.05;break;
		case '6' : act[5] += 0.05;break;
		case '7' : act[6] += 0.05;break;
		case '8' : act[7] += 0.05;break;
		case '9' : act[8] += 0.05;break;
		case '0' : act[9] += 0.05;break;
	
		case 'p' : PrintPosition();break;

		case 27: exit(0);break;
		default : break;
	}

	for(int i=0;i<act.rows();i++)
		if(act[i]>1.0)
			act[i] =1.0;
	mOctopus->SetActivationLevel(act);
	glutPostRedisplay();
}

Target target;
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

		target.idx = mSoftWorld->GetClosestNode(mouse_world);
		// std::cout << "DOWN: " << target.idx << std::endl;

		if(mMouseMode==CONSTRAINT_CONTROL)
		{
			/* Simulation */
			// int closest_node = mSoftWorld->GetClosestNode(mouse_world);
			// if(closest_node>=0)
			// 
				// mDragConstraint->GetI0() = closest_node;
				// mDragConstraint->GetP() = mouse_world;
				// mSoftWorld->AddConstraint(mDragConstraint);
				// std::cout<<closest_node<<std::endl;
			// }	
		}
	}
	else
	{
		if(mMouseMode==CONSTRAINT_CONTROL) {
			auto mouse_world = mCamera->GetWorldPosition(x,y);
			
			target.coord = mouse_world;
			// std::cout << "UP: " << target.idx << std::endl;
			mOctopus->AddTarget(target);
			mOctopus->SolveSoftIK(mSoftWorld);
		}

		mIsDrag = false;
		mMouseType = 0;
		mSoftWorld->RemoveConstraint(mDragConstraint);

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
		mDragConstraint->GetP() = mouse_world;
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

