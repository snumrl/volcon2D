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
	mIsPlay(false),mIsReplay(false),mIsPaused(false),mTime(0.0),mSimTime(0.0),mRecordFrame(0)
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
	mRigidWorld->setGravity(Eigen::Vector3d(0,-9.81,0));
	mSoftWorld = new FEM::World(
		// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
		1.0/100.0,										//time_step
		50, 											//max_iteration	
		0.999											//damping_coeff
		);
	MakeSkeleton(mMusculoSkeletalSystem);
	MakeMuscles("../vmcon2D/export/muscle_parameter.xml",mMusculoSkeletalSystem);

	mRigidWorld->addSkeleton(mMusculoSkeletalSystem->GetSkeleton());
	for(int i =0;i<3;i++)
	{
		mBalls.push_back(Skeleton::create("Ball_"+std::to_string(i)));
		MakeBall(mBalls.back(),0.03,0.6);	
		auto pos = mBalls.back()->getPositions();
		// std::cout<<pos.transpose()<<std::endl;
		pos.tail(3) = Eigen::Vector3d(i*0.1,-0.3,0);
		mBalls.back()->setPositions(pos);
	}
	for(auto& ball : mBalls)
		mRigidWorld->addSkeleton(ball);
	mMusculoSkeletalSystem->Initialize(mSoftWorld);
	mSoftWorld->Initialize();

	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	mRestPose = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	mPreviousPose = mRestPose;
	mTargetPositions = mRestPose;

	double kp = 500.0;
	double kv = 2*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);
}
bool
SimulationWindow2D::
TimeStepping()
{
	//Compute Optimal Activation Level
	// Compute desired Torque : PD control
	bool is_fem_updated = false;
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

	if(mSoftWorld->GetTime()<=mTime)
	{
		is_fem_updated =true;

		Eigen::VectorXd qdd_desired = ComputePDForces();

		static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

		if(mTime == 0.0){
			mMuscleOptimizationSolver->Initialize();
			mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);	
		}
		else
			mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	
		
		// std::cout<<"Desired Torque :\n"<<qdd_desired.transpose()<<std::endl;

		Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
		Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
		Eigen::VectorXd activation = solution.tail(mMusculoSkeletalSystem->GetNumMuscles()); 
		// std::cout<<"qdd :\n"<<qdd.transpose()<<std::endl;

		mMusculoSkeletalSystem->SetActivationLevel(activation);	
	}

	mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);

	mPreviousPose = skel->getPositions();
	
	if(is_fem_updated)
	{
		mMusculoSkeletalSystem->TransformAttachmentPoints();
		mSoftWorld->TimeStepping();
	}
	mRigidWorld->step();

		

	mRecords.push_back(new Record());
	auto rec = mRecords.back();
	rec->time = mTime;
	for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
		rec->rigid_body_positions.push_back(mRigidWorld->getSkeleton(i)->getPositions());
	rec->soft_body_positions = mSoftWorld->GetPositions();
	rec->activation_levels = mMusculoSkeletalSystem->GetActivationLevel();
	for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
		rec->muscle_forces.push_back(std::make_pair(muscle->force_origin,muscle->force_insertion));
	mTime+=mTimeStep;
	// std::cout<<std::endl<<std::endl;
	return is_fem_updated;
}
void
SimulationWindow2D::
SetRecord(Record* rec)
{
	mTime = rec->time;
	for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
		mRigidWorld->getSkeleton(i)->setPositions(rec->rigid_body_positions[i]);

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

Eigen::VectorXd
SimulationWindow2D::
ComputePDForces()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();

	Eigen::VectorXd pos_m = mTargetPositions;
	Eigen::VectorXd vel_m(pos_m.rows());
	vel_m.setZero();

	Eigen::VectorXd pos = skel->getPositions();
	Eigen::VectorXd vel = skel->getVelocities();

	Eigen::VectorXd pos_diff =skel->getPositionDifferences(pos_m,pos);
	for(int i = 0;i<pos_diff.rows();i++)
		pos_diff[i] = dart::math::wrapToPi(pos_diff[i]);
	std::cout<<"Pose diff : "<<pos_diff.transpose()<<std::endl;
	Eigen::VectorXd qdd_desired = 
				pos_diff.cwiseProduct(mKp) +
				(vel_m - vel).cwiseProduct(mKv);

	return qdd_desired;
}

Eigen::VectorXd
SimulationWindow2D::
SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap)
{
	auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	Eigen::VectorXd curr_pos = skel->getPositions();

	Eigen::VectorXd result = skel->getPositions();
	for(std::size_t i =0;i<1000;i++)
	{
		dart::math::LinearJacobian J = skel->getLinearJacobian(ap.first,ap.second);
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d inv_singular_value;
		
		inv_singular_value.setZero();
		for(int k=0;k<3;k++)
		{
			if(svd.singularValues()[k]==0)
				inv_singular_value(k,k) = 0.0;
			else
				inv_singular_value(k,k) = 1.0/svd.singularValues()[k];
		}


		Eigen::MatrixXd J_inv= svd.matrixV()*inv_singular_value*svd.matrixU().transpose();

		Eigen::Vector3d dir;
		
		Eigen::Vector2d T_curr = GetPoint(ap);
		
		dir = target_position-Eigen::Vector3d(T_curr[0],T_curr[1],0);
		
		double step_size = 2.0;
		double prev_energy = (dir.norm());
		//Line Minimization
		for(int k=0;k<24;k++)
		{
			step_size*=0.5;
			Eigen::Vector3d reduced_dir = step_size*dir;
			Eigen::VectorXd delta = J_inv*reduced_dir;

			skel->setPositions(result+delta);
			skel->computeForwardKinematics(true,false,false);

			Eigen::Vector3d cur_dir;

			Eigen::Vector2d point = GetPoint(ap);
			cur_dir = (target_position-Eigen::Vector3d(point[0],point[1],0));

			if(prev_energy>cur_dir.norm())
				break;
		}

		result = skel->getPositions();

		if(dir.norm()<1E-4)//||step_size<1E-6)
		{
			//std::cout<<dir.norm()<<std::endl;
			break;
		}
	}

	skel->setPositions(curr_pos);

	return result;
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
	for(auto& ball : mBalls)
		DrawSkeleton(ball);
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
		case 'p' : 
			if(!mIsReplay)
			{
				mIsReplay = true; mIsPlay = false; mIsPaused = true; mSimTime = mTime;
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
			mDragAnchorPoint = std::make_pair(target_bn,Eigen::Vector3d(0,0,0));
			mTargetPositions = SolveIK(Eigen::Vector3d(mouse_world[0],mouse_world[1],0),mDragAnchorPoint);
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
		mTargetPositions = SolveIK(Eigen::Vector3d(mouse_world[0],mouse_world[1],0),mDragAnchorPoint);
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

