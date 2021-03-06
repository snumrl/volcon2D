#include "MuscleOptimization.h"
#include "MusculoSkeletalSystem.h"

using namespace Ipopt;
MuscleOptimization::
MuscleOptimization(FEM::World* soft_world,const dart::simulation::WorldPtr&	rigid_world, MusculoSkeletalSystem* ms)
	:mSoftWorld(soft_world),mRigidWorld(rigid_world),mMusculoSkeletalSystem(ms),mWeightTracking(1.0),mWeightEffort(1.0),mWeightReinforceEndEffector(0.0),mSparseUpdateCount(0)
{
	int num_muscles =  mMusculoSkeletalSystem->GetNumMuscles();
	int dofs 		=  mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();  
	mSolution.resize(dofs +num_muscles);
	mQddDesired.resize(dofs);

	mJt.resize(dofs,(num_muscles*2)*2);			//first 2 : # of forces(insertion, origin) , last 2 : dimension(since we are in 2D)
	mA.resize(num_muscles*2*2,num_muscles);
	mP.resize(num_muscles*2*2);

	mM_minus_JtA.resize(dofs,dofs+num_muscles);
	mJtp_minus_c.resize(dofs);	
	mSolution.setZero();
	mQddDesired.setZero();

	mJt.setZero();
	mA.setZero();
	mP.setZero();

	mM_minus_JtA.setZero();
	mJtp_minus_c.setZero();

	// mEndEffector.push_back(std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL"),Eigen::Vector3d(0,0,0)));
	// mEndEffector.push_back(std::make_pair(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"),Eigen::Vector3d(0,0,0)));

	mYddDesired.resize(2*mEndEffector.size());
	mdydd_dqdd.resize(2*mEndEffector.size(),dofs);

	mYddDesired.setZero();
	mdydd_dqdd.setZero();

}
MuscleOptimization::
~MuscleOptimization()
{

}
void
MuscleOptimization::
Update(const Eigen::VectorXd& qdd_desired)
{
	mQddDesired = qdd_desired;

	//Update Jt
	int index = 0;
	auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	int dofs = skel->getNumDofs();
	for(auto& muscle : mMusculoSkeletalSystem->GetMuscles())
	{
		mJt.block(0,index*4,dofs,2) = (skel->getLinearJacobian(muscle->originWayPoints.back().first,muscle->originWayPoints.back().second).block(0,0,2,dofs)).transpose();
		mJt.block(0,index*4 + 2,dofs,2) = (skel->getLinearJacobian(muscle->insertionWayPoints.back().first,muscle->insertionWayPoints.back().second).block(0,0,2,dofs)).transpose();
		index++;
	}
	int num_muscles =  mMusculoSkeletalSystem->GetNumMuscles();
	mA = mMusculoSkeletalSystem->ComputeForceDerivative(mSoftWorld);
	mP = mMusculoSkeletalSystem->ComputeForce(mSoftWorld);

	mP = mP-mA*mMusculoSkeletalSystem->GetActivationLevel();

	mM_minus_JtA.setZero();
	mM_minus_JtA.block(0,0,dofs,dofs)= mMusculoSkeletalSystem->GetSkeleton()->getMassMatrix();
	mM_minus_JtA.block(0,dofs,dofs,num_muscles)= -mJt*mA;
	mJtp_minus_c = mJt*mP - mMusculoSkeletalSystem->GetSkeleton()->getCoriolisAndGravityForces();

	//Compute Endeffector Things

	// ComputeYdd(mQddDesired,mYddDesired);
	// ComputedYdd();
}

const Eigen::VectorXd&
MuscleOptimization::
GetSolution()
{
	return mSolution;
}

void
MuscleOptimization::
UpdateConstraints(const Eigen::VectorXd& act)
{
	auto& prev_act = mMusculoSkeletalSystem->GetActivationLevel();
	int num_muscles =  mMusculoSkeletalSystem->GetNumMuscles();
	int dofs 		=  mMusculoSkeletalSystem->GetSkeleton()->getNumDofs(); 
	mSparseUpdateCount++;
	if(mSparseUpdateCount==3)
	if( (act-prev_act).norm()>1E-5)
	{
		mMusculoSkeletalSystem->SetActivationLevel(act);
		mSoftWorld->TimeStepping(false);
		//Update A,p
		mA = mMusculoSkeletalSystem->ComputeForceDerivative(mSoftWorld);
		mP = mMusculoSkeletalSystem->ComputeForce(mSoftWorld);

		mP = mP-mA*act;
		// std::cout<<"A : "<<mA.transpose()<<std::endl;
		// std::cout<<"p : "<<mP.transpose()<<std::endl<<std::endl;
		//Update Cache

		mM_minus_JtA.block(0,0,dofs,dofs)= mMusculoSkeletalSystem->GetSkeleton()->getMassMatrix();
		mM_minus_JtA.block(0,dofs,dofs,num_muscles)= -mJt*mA;
		mJtp_minus_c = mJt*mP - mMusculoSkeletalSystem->GetSkeleton()->getCoriolisAndGravityForces();
		mSparseUpdateCount = 0;
	}
}
bool
MuscleOptimization::
get_nlp_info(	Index& n, Index& m, Index& nnz_jac_g,Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	m = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	n = m + mMusculoSkeletalSystem->GetNumMuscles();

	nnz_jac_g = n*m;			//g : full matrix
	nnz_h_lag = m*m+(n-m);				//H : full matrix

	index_style = TNLP::C_STYLE;
}
bool
MuscleOptimization::
get_bounds_info(Index n, Number* x_l, Number* x_u,Index m, Number* g_l, Number* g_u)
{
	for(int i=0;i<m;i++)
	{
		x_l[i] =-1E5;
		x_u[i] =1E5;
	}
	for(int i=m;i<n;i++)
	{
		x_l[i] = 0.0;
		x_u[i] = 1.0;
	}
	for(int i =0;i<m;i++)
		g_l[i]=g_u[i] =0.0;

	return true;
}
bool
MuscleOptimization::
get_starting_point(	Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
	for(int i =0;i<n;i++)
		x[i] = mSolution[i];

	return true;
}
bool
MuscleOptimization::
eval_f(	Index n, const Number* x, bool new_x, Number& obj_value)
{
	double track = 0.0;
	double effort = 0.0;
	double ee = 0.0;

	int m = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	Eigen::VectorXd qdd(m),ydd(2*mEndEffector.size());

	for(int i=0;i<m;i++)
		qdd[i] = x[i];

	track = mWeightTracking*((qdd-mQddDesired).squaredNorm());

	for(int i=m;i<n;i++)
		effort += x[i]*x[i];
	effort *= mWeightEffort;

	ComputeYdd(qdd,ydd);
	ee = mWeightReinforceEndEffector*((ydd-mYddDesired).squaredNorm());

	// std::cout<<"Track : "<<track<<std::endl;
	// std::cout<<"effort : "<<effort<<std::endl;
	// std::cout<<"ee : "<<ee<<std::endl;
	// std::cout<<std::endl;
	obj_value = track + effort + ee;

	return true;
}
bool
MuscleOptimization::
eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	int m = mMusculoSkeletalSystem->GetSkeleton()->getNumDofs();
	Eigen::VectorXd qdd(m),ydd(2*mEndEffector.size());
	for(int i=0;i<m;i++)
		qdd[i] = x[i];

	for(int i=0;i<m;i++)
		grad_f[i] = 2.0*mWeightTracking*(x[i]-mQddDesired[i]);
	for(int i=m;i<n;i++)
		grad_f[i] = 2.0*mWeightEffort*(x[i]);

	ComputeYdd(qdd,ydd);

	auto g_ee = 2.0*mWeightReinforceEndEffector*(mdydd_dqdd.transpose())*(ydd-mYddDesired);

	for(int i =0;i<m;i++)
		grad_f[i] += g_ee[i];

	return true;
}
bool
MuscleOptimization::
eval_g(	Index n, const Number* x, bool new_x, Index m, Number* g)
{
	Eigen::VectorXd eigen_x(n), activation(n-m),eigen_g(m);
	for(int i=0;i<n;i++)
		eigen_x[i] = x[i];

	activation = eigen_x.tail(n-m);
	UpdateConstraints(activation);
	// std::cout<<"act : "<<activation<<std::endl;
	eigen_g = mM_minus_JtA*eigen_x - mJtp_minus_c;
	for(int i = 0;i<m;i++)
		g[i] = eigen_g[i];

	return true;
}
bool
MuscleOptimization::
eval_jac_g( Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
	int nnz = 0;

	if(values == NULL)
	{
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<n;j++)
			{
				iRow[nnz] = i;
				jCol[nnz++] = j;
			}
		}
	}
	else
	{
		Eigen::VectorXd eigen_x(n), activation(n-m),eigen_g(m);
		for(int i=0;i<n;i++)
			eigen_x[i] = x[i];

		activation = eigen_x.tail(n-m);
		UpdateConstraints(activation);
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<n;j++)
			{
				values[nnz++] = mM_minus_JtA(i,j);
			}
		}
	}

	return true;

}
bool
MuscleOptimization::
eval_h( Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	int nnz = 0;

	if(values == NULL)
	{
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<m;j++)
			{
				iRow[nnz] = i;
				jCol[nnz++] = j;
			}
		}

		for(int i =m;i<n;i++)
		{
			iRow[nnz] = i;
			jCol[nnz++] = i;
		}
	}
	else
	{
		Eigen::MatrixXd dyddTdydd = 2.0*mWeightReinforceEndEffector*(mdydd_dqdd.transpose() * mdydd_dqdd);
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<m;j++)
			{
				if(i==j)
					values[nnz++] = dyddTdydd(i,j) + 2.0*obj_factor*mWeightTracking;
				else
					values[nnz++] = dyddTdydd(i,j);
			}
		}
		for(int i=m;i<n;i++)
			values[nnz++] = 2.0*obj_factor*mWeightEffort;
	}

	return true;
}
void
MuscleOptimization::
finalize_solution(	SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda,Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq)
{
	for(int i=0;i<n;i++)
		mSolution[i] = x[i];
}





void
MuscleOptimization::
ComputeYdd(const Eigen::VectorXd& qdd,Eigen::VectorXd& ydd)
{
	ydd.setZero();
	auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	int index = 0;
	for(auto& ee : mEndEffector)
	{
		auto* bn = ee.first;
		auto  xb = ee.second;
		auto  J  = skel->getJacobian(bn);
		auto  Jqdd = J * qdd;
		Eigen::Affine3d T;
		T.linear() 		= dart::math::makeSkewSymmetric(Jqdd.head(3));
		T.translation()	= Jqdd.tail(3);

		ydd.block<2,1>(index*2,0) = (T*xb).block<2,1>(0,0);
		index++;
	}
}
void
MuscleOptimization::
ComputedYdd()
{
	auto& skel = mMusculoSkeletalSystem->GetSkeleton();
	mdydd_dqdd.setZero();
	int index = 0;
	for(auto& ee : mEndEffector)
	{
		auto* bn = ee.first;
		auto  xb = ee.second;
		auto  J  = skel->getJacobian(bn);
		for(int i =0;i<J.cols();i++)
		{
			Eigen::Affine3d T;
			T.linear() 		= dart::math::makeSkewSymmetric(J.block<3,1>(0,i));
			T.translation()	= J.block<3,1>(3,i);

			mdydd_dqdd.block<2,1>(index*2,i) = (T*xb).block<2,1>(0,0);
		}
		index++;
	}
}
