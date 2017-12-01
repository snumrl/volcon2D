#include "VelocityControlDDP.h"
VelocityControlDDP::
VelocityControlDDP(const dart::simulation::WorldPtr& rigid_world,
		FEM::World* soft_world,
		MusculoSkeletalSystem* mMusculoSkeletalSystem,int n,int max_iteration)
		:MusculoSkeletalDDP(rigid_world,soft_world,mMusculoSkeletalSystem,n,max_iteration),
		wa(1E3),wb(1E4),wc(1E3)
{
	Eigen::VectorXd x0(mDofs*2);
	GetState(x0);
	std::vector<Eigen::VectorXd> u0(mN-1,Eigen::VectorXd::Zero(mSu));
	
	dart::math::seedRand();	
	for(int i =0;i<mN-1;i++)
	{
		u0[i] = dart::math::randomVectorXd(mSu,0,0.001);
	}
	Init(x0,u0,Eigen::VectorXd::Constant(mSu,0.0),Eigen::VectorXd::Constant(mSu,1.0));

	x_desired.setZero();

	x_desired[0] = 0.5;
	x_desired[1] = 0.3;
	x_ee_local = Eigen::Vector3d(0,0,0);
}

void
VelocityControlDDP::
EvalCf(const Eigen::VectorXd& x,double& cf)
{
	SetState(x);

	auto x_ee = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getTransform()*x_ee_local;
	cf = 0.5*wb*(x_ee-x_desired).squaredNorm();
	cf += 0.5*wc*(mMusculoSkeletalSystem->GetSkeleton()->getVelocities()).squaredNorm();
}
void
VelocityControlDDP::
EvalCfx(const Eigen::VectorXd& x,Eigen::VectorXd& cfx)
{
	SetState(x);
	
	auto bn_ee = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR");
	auto x_ee = bn_ee->getTransform()*x_ee_local;
	auto J = mMusculoSkeletalSystem->GetSkeleton()->getLinearJacobian(bn_ee,x_ee_local);

	cfx.block(0,0,mDofs,1) = wb*J.transpose()*(x_ee-x_desired);
	cfx.block(mDofs,0,mDofs,1) = wc*(mMusculoSkeletalSystem->GetSkeleton()->getVelocities());
}
void
VelocityControlDDP::
EvalCfxx(const Eigen::VectorXd& x,Eigen::MatrixXd& cfxx)
{
	cfxx.block(0,0,mDofs,mDofs) = wb*Eigen::MatrixXd::Identity(mDofs,mDofs);
	cfxx.block(mDofs,mDofs,mDofs,mDofs) = wc*Eigen::MatrixXd::Identity(mDofs,mDofs);
}
void
VelocityControlDDP::
EvalC( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,double& c)
{
	c = 0.5*wa*u.squaredNorm();
}
void
VelocityControlDDP::
EvalCx( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cx)
{
	cx.setZero();
}
void
VelocityControlDDP::
EvalCu( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cu)
{
	cu = wa*u;
}
void
VelocityControlDDP::
EvalCxx(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxx)
{
	cxx.setZero();
}
void
VelocityControlDDP::
EvalCxu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxu)
{
	cxu.setZero();
}
void
VelocityControlDDP::
EvalCuu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cuu)
{
	cuu = wa*Eigen::MatrixXd::Identity(u.rows(),u.rows());
}
void
VelocityControlDDP::
Evalf(  const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& f)
{
	SetState(x);
	SetControl(u);
	Step();
	GetState(f);
}
void
VelocityControlDDP::
Evalfx( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& fx)
{
	Eigen::VectorXd fx_i_plus(x.rows());
	Eigen::VectorXd fx_i_minus(x.rows());
	Eigen::VectorXd x_i;
	double delta = 0.01;
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;

		x_i[i] = x[i] + delta;
		SetState(x_i);
		// SetControl(u);
		Step(false);
		GetState(fx_i_plus);

		x_i[i] = x[i] - delta;
		SetState(x_i);
		// SetControl(u);
		Step(false);
		GetState(fx_i_minus);

		fx.col(i) = (fx_i_plus - fx_i_minus)/(2*delta);
	}
	
}


void
VelocityControlDDP::
Evalfu( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& fu)
{
	Eigen::VectorXd fu_i_plus(x.rows());
	Eigen::VectorXd fu_i_minus(x.rows());
	Eigen::VectorXd u_i;
	double delta = 0.01;
	for(int i = 0;i<mSu;i++)
	{
		u_i = u;
		double u_plus,u_minus;
		u_plus = u[i] + delta;
		u_minus = u[i] - delta;

		if(u_plus>1.0)
			u_plus = 1.0;
		if(u_minus<0)
			u_minus = 0.0;

		u_i[i] = u_plus;
		SetState(x);
		SetControl(u_i);
		Step();
		GetState(fu_i_plus);

		u_i[i] = u_minus;
		SetState(x);
		SetControl(u_i);
		Step();
		GetState(fu_i_minus);

		fu.col(i) = (fu_i_plus - fu_i_minus)/(u_plus - u_minus);
	}
}