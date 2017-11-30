#include "SimulationWindow2D.h"
#include "MusculoSkeletalSystem.h"
#include "GUI/Camera2D.h"
#include "GUI/GL_function.h"
#include "FEM2D_Interface.h"
#include "DART_Interface.h"
#include "Record.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "Controller.h"
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;

SimulationWindow2D::
SimulationWindow2D()
	:mMouseMode(MOUSE_MODE::CONSTRAINT_CONTROL),
	mDragConstraint(new FEM::AttachmentConstraint(50000.0,0,Eigen::Vector2d(0,0))),
	mMusculoSkeletalSystem(new MusculoSkeletalSystem()),mController(new Controller()),
	mIsPlay(false),mIsReplay(false),mIsPaused(false),mSimTime(0.0),mRecordFrame(0)
{
	Initialize();
	mDisplayTimeout = mSoftWorld->GetTimeStep()*1000;
}

void
SimulationWindow2D::
Initialize()
{
	mRigidWorld = std::make_shared<World>();
	mRigidWorld->setGravity(Eigen::Vector3d(0,-9.81,0));
	mSoftWorld = new FEM::World(
		// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/200.0,										//time_step
		50, 											//max_iteration	
		0.999											//damping_coeff
		);
	MakeSkeleton(mMusculoSkeletalSystem);
	MakeMuscles("../vmcon2D/export/muscle_parameter.xml",mMusculoSkeletalSystem);

	mRigidWorld->addSkeleton(mMusculoSkeletalSystem->GetSkeleton());
	for(int i =0;i<5;i++)
	{
		mBalls.push_back(Skeleton::create("Ball_"+std::to_string(i)));
		MakeBall(mBalls.back(),0.036,0.13);	
		auto pos = mBalls.back()->getPositions();
		// std::cout<<pos.transpose()<<std::endl;
		pos.tail(3) = Eigen::Vector3d(i*0.1,-0.3,0);
		mBalls.back()->setPositions(pos);
	}
	for(auto& ball : mBalls)
		mRigidWorld->addSkeleton(ball);
	mMusculoSkeletalSystem->Initialize(mSoftWorld);
	mSoftWorld->Initialize();

	mController->Initialize(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem,mBalls);
	mDragAnchorPoint = std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode(0),Eigen::Vector3d(0,0,0));

}
bool
SimulationWindow2D::
TimeStepping()
{
	bool is_fem_updated = false;
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();


	//Simulation Loop
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

		
	//Record Loop
	mRecords.push_back(new Record());
	auto rec = mRecords.back();
	rec->Set(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);
	// rec->time = mTime;
	// for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
	// 	rec->rigid_body_positions.push_back(mRigidWorld->getSkeleton(i)->getPositions());

	// for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
	// 	rec->rigid_body_velocities.push_back(mRigidWorld->getSkeleton(i)->getVelocities());

	// rec->soft_body_positions = mSoftWorld->GetPositions();
	// rec->activation_levels = mMusculoSkeletalSystem->GetActivationLevel();
	// for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
	// 	rec->muscle_forces.push_back(std::make_pair(muscle->force_origin,muscle->force_insertion));

	// rec->state = mController->GetMachine()->GetCurrentState();
	
	// std::cout<<std::endl<<std::endl;
	return is_fem_updated;
}
void
SimulationWindow2D::
SetRecord(Record* rec)
{
	rec->Get(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem,mController);
	
	// for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
	// 	mRigidWorld->getSkeleton(i)->setPositions(rec->rigid_body_positions[i]);
	// for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
	// 	mRigidWorld->getSkeleton(i)->setVelocities(rec->rigid_body_velocities[i]);

	// mSoftWorld->SetPositions(rec->soft_body_positions);
	// int count = 0;
 //    for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
 //    {
 //        muscle->activationLevel = rec->activation_levels[count];
 //        for(auto& mc : muscle->muscleConstraints)
 //            mc->SetActivationLevel(muscle->activationLevel);
 //        muscle->force_origin = rec->muscle_forces[count].first;
 //        muscle->force_insertion = rec->muscle_forces[count].second;
 //        muscle->origin->GetP() = GetPoint(muscle->originWayPoints[0]);
 //        muscle->insertion->GetP() = GetPoint(muscle->insertionWayPoints[0]);
 //        count++;
 //    }
 //    mController->GetMachine()->SetCurrentState(rec->state);
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
    DrawStringOnScreen(0.8,0.2,std::to_string(mRigidWorld->getTime()),true);
	glLineWidth(1.0);

	const auto& x = mSoftWorld->GetPositions();
	
	for(auto& c : mSoftWorld->GetConstraints())
	{
		DrawConstraint(c,x);
	}
	if(mIsDrag)
		DrawConstraint(mDragConstraint,x);
	
	for(auto& muscle :mMusculoSkeletalSystem->GetMuscles())
		DrawMuscle(muscle,x);

	DrawSkeleton(mMusculoSkeletalSystem->GetSkeleton());
	if(mIsPlay)
	{
		auto cur_pos = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
		mMusculoSkeletalSystem->GetSkeleton()->setPositions(mController->GetTargetPositions());
		mMusculoSkeletalSystem->GetSkeleton()->computeForwardKinematics(true,false,false);
		DrawSkeleton(mMusculoSkeletalSystem->GetSkeleton(),Eigen::Vector3d(0.8,0.3,0.8));
		mMusculoSkeletalSystem->GetSkeleton()->setPositions(cur_pos);
		mMusculoSkeletalSystem->GetSkeleton()->computeForwardKinematics(true,false,false);
	}
	int ball_index = 0;
	for(auto& ball : mBalls)
	{
		if(ball_index==0)
		DrawSkeleton(ball,Eigen::Vector3d(0.8,0.4,0.4));
		else if(ball_index==1)
		DrawSkeleton(ball,Eigen::Vector3d(0.4,0.8,0.4));
		else if(ball_index==2)
		DrawSkeleton(ball,Eigen::Vector3d(0.4,0.4,0.8));
		else if(ball_index==3)
		DrawSkeleton(ball,Eigen::Vector3d(0.4,0.8,0.8));
		else if(ball_index==4)
		DrawSkeleton(ball,Eigen::Vector3d(0.8,0.4,0.8));
		ball_index++;
	}

	glEnable(GL_DEPTH_TEST);

	glutSwapBuffers();
}
void
SimulationWindow2D::
Keyboard(unsigned char key,int x,int y)
{
	auto& act = mMusculoSkeletalSystem->GetActivationLevel();
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
	// 	break;
	// 	case 'a' : mSimulator->PrintTimer();break;
		case ' ' : mIsPlay = !mIsPlay; break;
		case 'p' : 
			if(!mIsReplay)
			{
				mIsReplay = true; mIsPlay = false; mIsPaused = true; mSimTime = mRigidWorld->getTime();
			}
			else if(mIsPaused)
				mIsPaused = false;
			else
				mIsPaused = true;
			break;
		case '[' : mRecordFrame--; break;
		case ']' : mRecordFrame++; break;
	// 	case 'p' : mSimulator->SetPlay(false);
	// 			   mIsReplay = true;
	// 			    break;
		case 'r' : act.setZero();break;
		case '1' : mDragAnchorPoint.first = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL");std::cout<<"Left Hand"<<std::endl;break;
		case '2' : mDragAnchorPoint.first = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");std::cout<<"Right Hand"<<std::endl;break;
		case 'g' : mMusculoSkeletalSystem->GetSkeleton()->getBodyNode(3)->addExtForce(Eigen::Vector3d::Random()*100,Eigen::Vector3d::Random());break;
		case 27: exit(0);break;
		default : break;
	}

	for(int i=0;i<act.rows();i++)
		if(act[i]>1.0)
			act[i] =1.0;
	mMusculoSkeletalSystem->SetActivationLevel(act);
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
			// int closest_node = mSoftWorld->GetClosestNode(mouse_world);
			// if(closest_node>=0)
			// {
			// 	mDragConstraint->GetI0() = closest_node;
			// 	mDragConstraint->GetP() = mouse_world;
			// 	mSoftWorld->AddConstraint(mDragConstraint);
			// 	std::cout<<closest_node<<std::endl;
			// }
			auto& skel = mMusculoSkeletalSystem->GetSkeleton();
			double min_distance = 1E6;
			BodyNode* target_bn = nullptr;
			for(auto& bn : skel->getBodyNodes())
			{
				const auto& T =bn->getTransform();
				double distance = ( Eigen::Vector3d(mouse_world[0],mouse_world[1],0)-T.translation()).norm();
				if(min_distance>distance)
				{
					target_bn = bn;
					min_distance = distance;
				}
			}
			
		}
	}
	else
	{
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
		// mDragConstraint->GetP() = mouse_world;

		// mController->SetTargetPositions(mController->SolveIK(Eigen::Vector3d(mouse_world[0],mouse_world[1],0),mDragAnchorPoint));
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
	if(mIsPlay){
		while(!TimeStepping()){
		}
		
	}
	else if(mIsReplay){

		if(mIsPaused)
		{
			if(mRecordFrame < 0)
				mRecordFrame = 0;
			if(mRecordFrame>=mRecords.size())
				mRecordFrame=mRecords.size()-1;
			

			SetRecord(mRecords[mRecordFrame]);

			mDisplayTimeout = 1;
		}
		else
		{
			SetRecord(mRecords[mRecordFrame]);
			mRecordFrame+=33;
			if(mRecordFrame>=mRecords.size())
				mRecordFrame=0;
			mDisplayTimeout = 1000/30;	
		}
		
		
	}
	glutPostRedisplay();
	glutTimerFunc(mDisplayTimeout, TimerEvent,1);
}

