#include "Controller.h"
#include "MusculoSkeletalSystem.h"
#include "IKOptimization.h"
#include "FSM_Interface.h"
#include "MuscleOptimization.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;

#include "FSM.h"
#include "FSM_Interface.h"
using namespace dart::simulation;
using namespace dart::dynamics;
Controller::
Controller()
	:mFSM(nullptr)
{
	
}

void
Controller::
Initialize(FEM::World* soft_world,const WorldPtr& rigid_world,MusculoSkeletalSystem* musculo_skeletal_system)
{
	mSoftWorld = soft_world;
	mRigidWorld = rigid_world;
	mMusculoSkeletalSystem = musculo_skeletal_system;

	mMuscleOptimization = new MuscleOptimization(mSoftWorld,mRigidWorld,mMusculoSkeletalSystem);
	mMuscleOptimizationSolver = new IpoptApplication();
	mIKOptimization = new IKOptimization(mMusculoSkeletalSystem->GetSkeleton());
	mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "no");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-4);

	mIKSolver = new IpoptApplication();

	mIKSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mIKSolver->Options()->SetStringValue("jac_c_constant", "yes");
	mIKSolver->Options()->SetStringValue("hessian_constant", "yes");
	mIKSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mIKSolver->Options()->SetIntegerValue("print_level", 2);
	mIKSolver->Options()->SetIntegerValue("max_iter", 1000);
	mIKSolver->Options()->SetNumericValue("tol", 1e-4);

	mRestPose = mMusculoSkeletalSystem->GetSkeleton()->getPositions();
	mPreviousPose = mRestPose;
	mTargetPositions = mRestPose;

	double kp = 200.0;
	double kv = 2*sqrt(kp);
	int n = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(n,kp);
	mKv = Eigen::VectorXd::Constant(n,kv);

	mFSM = new Machine();
	MakeMachine("../vmcon2D/export/juggling.xml",mFSM);
}
const Eigen::VectorXd&
Controller::
GetTargetPositions()
{
	return mTargetPositions;
}
void
Controller::
SetTargetPositions(const Eigen::VectorXd& tp)
{
	mTargetPositions = tp;
}
Eigen::VectorXd
Controller::
ComputeActivationLevels()
{
	auto& skel =mMusculoSkeletalSystem->GetSkeleton();
	Eigen::VectorXd qdd_desired = ComputePDForces();

	static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

	if(mSoftWorld->GetTime() == 0.0){
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

	return activation;
}
Eigen::VectorXd
Controller::
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
	Eigen::VectorXd qdd_desired = 
				pos_diff.cwiseProduct(mKp) +
				(vel_m - vel).cwiseProduct(mKv);

	return qdd_desired;
}

Eigen::VectorXd
Controller::
SolveIK(const Eigen::Vector3d& target_position,AnchorPoint ap)
{
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	ik->AddTargetPositions(ap,target_position);
	mIKSolver->Initialize();
	mIKSolver->OptimizeTNLP(mIKOptimization);

	Eigen::VectorXd solution = ik->GetSolution();
	return solution;
	// auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	// Eigen::VectorXd curr_pos = skel->getPositions();

	// Eigen::VectorXd result = skel->getPositions();
	// for(std::size_t i =0;i<1000;i++)
	// {
	// 	dart::math::LinearJacobian J = skel->getLinearJacobian(ap.first,ap.second);
	// 	Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
	// 	Eigen::Matrix3d inv_singular_value;
		
	// 	inv_singular_value.setZero();
	// 	for(int k=0;k<3;k++)
	// 	{
	// 		if(svd.singularValues()[k]==0)
	// 			inv_singular_value(k,k) = 0.0;
	// 		else
	// 			inv_singular_value(k,k) = 1.0/svd.singularValues()[k];
	// 	}


	// 	// Eigen::MatrixXd J_inv= svd.matrixV()*inv_singular_value*svd.matrixU().transpose();
	// 	Eigen::MatrixXd J_inv = J.transpose();

	// 	Eigen::Vector3d dir;
		
	// 	Eigen::Vector2d T_curr = GetPoint(ap);
		
	// 	dir = target_position-Eigen::Vector3d(T_curr[0],T_curr[1],0);
		
	// 	double step_size = 2.0;
	// 	double prev_energy = (dir.norm());
	// 	//Line Minimization
	// 	for(int k=0;k<24;k++)
	// 	{
	// 		step_size*=0.5;
	// 		Eigen::Vector3d reduced_dir = step_size*dir;
	// 		Eigen::VectorXd delta = J_inv*reduced_dir;

	// 		skel->setPositions(result+delta);
	// 		skel->computeForwardKinematics(true,false,false);

	// 		Eigen::Vector3d cur_dir;

	// 		Eigen::Vector2d point = GetPoint(ap);
	// 		cur_dir = (target_position-Eigen::Vector3d(point[0],point[1],0));

	// 		if(prev_energy>cur_dir.norm())
	// 			break;
	// 	}

	// 	result = skel->getPositions();

	// 	if(dir.norm()<1E-4)//||step_size<1E-6)
	// 	{
	// 		//std::cout<<dir.norm()<<std::endl;
	// 		break;
	// 	}
	// }

	// skel->setPositions(curr_pos);

	// return result;
}


Machine*
Controller::
GetMachine()
{
	return mFSM;
}