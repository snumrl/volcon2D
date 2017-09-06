#include "SimulationWindow2D.h"
#include "MusculoSkeletalSystem.h"
#include "GUI/Camera2D.h"
#include "GUI/GL_function.h"
#include "FEM2D_Interface.h"
#include "DART_Interface.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;

SimulationWindow2D::
SimulationWindow2D()
	:mMouseMode(MOUSE_MODE::CAMERA_CONTROL),
	mDragConstraint(new FEM::AttachmentConstraint(50000.0,0,Eigen::Vector2d(0,0))),
	mMusculoSkeletalSystem(new MusculoSkeletalSystem()),
	mIsPlay(false),mIsReplay(false),mTime(0.0),mSimTime(0.0),mRecordFrame(0)
{
	Initialize();
	mDisplayTimeout = mSoftWorld->GetTimeStep()*1000;
	mTimeStep = mRigidWorld->getTimeStep();
}

void
SimulationWindow2D::
Initialize()
{
	mRigidWorld = std::make_shared<World>();
	mSoftWorld = new FEM::World(
		// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/1000.0,										//time_step
		50, 											//max_iteration	
		0.999											//damping_coeff
		);
	MakeSkeleton(mMusculoSkeletalSystem);
	MakeMuscles("../vmcon2D/export/muscle_parameter.xml",mMusculoSkeletalSystem);

	mRigidWorld->addSkeleton(mMusculoSkeletalSystem->GetSkeleton());
	mMusculoSkeletalSystem->Initialize(mSoftWorld);
	mSoftWorld->Initialize();

	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	// mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 1);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 30);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	mRestPose = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	mPreviousPose = mRestPose;
}
void
SimulationWindow2D::
TimeStepping()
{
	//Compute Optimal Activation Level
	// Compute desired Torque : PD control
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

		Eigen::VectorXd pos_m = mRestPose;
		Eigen::VectorXd vel_m(pos_m.rows());
		vel_m.setZero();

		Eigen::VectorXd pos = skel->getPositions();
		Eigen::VectorXd vel = skel->getVelocities();
		std::cout<<"Pose diff : "<<skel->getPositionDifferences(pos_m,pos).transpose()<<std::endl;
		Eigen::VectorXd kp(pos_m.rows()),kv(pos_m.rows());
		for(int i =0;i<pos_m.rows();i++)
		{
			kp[i] = 500.0;
			kv[i] = 2*sqrt(kp[i]);
		}
		Eigen::VectorXd qdd_desired = 
					skel->getPositionDifferences(pos_m,pos).cwiseProduct(kp) +
					(vel_m - vel).cwiseProduct(kv);


		static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

		if(mTime == 0.0){
			mMuscleOptimizationSolver->Initialize();
			mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	
		}
		else
			mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	
		
		std::cout<<"Desired Torque :\n"<<qdd_desired.transpose()<<std::endl;

	Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
	Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
	Eigen::VectorXd activation = solution.tail(mMusculoSkeletalSystem->GetNumMuscles()); 
	std::cout<<"qdd :\n"<<qdd.transpose()<<std::endl;

	mMusculoSkeletalSystem->SetActivationLevel(activation);	
	mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);
	// skel->setForces(skel->getMassMatrix()*qdd+skel->getCoriolisAndGravityForces());

	mPreviousPose = skel->getPositions();

	if(mSoftWorld->GetTime()<mTime)
		mSoftWorld->TimeStepping();
	mRigidWorld->step();

	mMusculoSkeletalSystem->TransformAttachmentPoints();

	mRecords.push_back(new Record());
	auto rec = mRecords.back();
	rec->time = mTime;
	rec->rigid_body_positions = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	rec->soft_body_positions = mSoftWorld->GetPositions();
	rec->activation_levels = activation;
	for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
		rec->muscle_forces.push_back(std::make_pair(muscle->force_origin,muscle->force_insertion));
	mTime+=mTimeStep;
	std::cout<<std::endl<<std::endl;
}
void
SimulationWindow2D::
SetRecord(Record* rec)
{
	mTime = rec->time;
	mMusculoSkeletalSystem->GetSkeleton()->setPositions(rec->rigid_body_positions);
	mSoftWorld->SetPositions(rec->soft_body_positions);
	int count = 0;
    for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
    {
        muscle->activationLevel = rec->activation_levels[count];
        for(auto& mc : muscle->muscleConstraints)
            mc->SetActivationLevel(muscle->activationLevel);
        muscle->force_origin = rec->muscle_forces[count].first;
        muscle->force_insertion = rec->muscle_forces[count].second;
        muscle->origin->GetP() = GetPoint(muscle->originWayPoints[0]);
        muscle->insertion->GetP() = GetPoint(muscle->insertionWayPoints[0]);
        count++;
    }
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
    DrawStringOnScreen(0.8,0.2,std::to_string(mTime),true);
	glLineWidth(1.0);

	const auto& x = mSoftWorld->GetPositions();
	
	for(auto& c : mSoftWorld->GetConstraints())
	{
		DrawConstraint(c,x);
	}
	if(mIsDrag)
		DrawConstraint(mDragConstraint,x);
	DrawSkeleton(mMusculoSkeletalSystem->GetSkeleton());
	for(auto& muscle :mMusculoSkeletalSystem->GetMuscles())
		DrawMuscle(muscle,x);
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
		case 'p' : mIsReplay = true; mIsPlay = false; mSimTime = mTime;break;
	// 	case 'p' : mSimulator->SetPlay(false);
	// 			   mIsReplay = true;
	// 			    break;
		case 'r' : act.setZero();break;
		case '1' : act[1] += 0.1;break;
		case '2' : act[2] += 0.1;break;
		case '3' : act[3] += 0.1;break;
		case '4' : act[4] += 0.1;break;
		case '5' : act[5] += 0.1;break;
		case '6' : act[6] += 0.1;break;
		case '7' : act[7] += 0.1;break;
		case '8' : act[8] += 0.1;break;
		case '9' : act[9] += 0.1;break;
		case '0' : act[0] += 0.1;break;
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
			int closest_node = mSoftWorld->GetClosestNode(mouse_world);
			if(closest_node>=0)
			{
				mDragConstraint->GetI0() = closest_node;
				mDragConstraint->GetP() = mouse_world;
				mSoftWorld->AddConstraint(mDragConstraint);
				std::cout<<closest_node<<std::endl;
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
	if(mIsPlay){
		TimeStepping();
		
	}
	else if(mIsReplay){
		SetRecord(mRecords[mRecordFrame]);
		mRecordFrame+=33;
		if(mRecordFrame>=mRecords.size())
			mRecordFrame=0;
		mDisplayTimeout = 1000/30;
		
	}
	glutPostRedisplay();
	glutTimerFunc(mDisplayTimeout, TimerEvent,1);
}


