#include "VelocityControlDDP.h"
#include "../BallInfo.h"
#include <fstream>
static std::ifstream param("../vmcon2D/export/param.txt");
VelocityControlDDP::
VelocityControlDDP(const dart::simulation::WorldPtr& rigid_world,
		FEM::World* soft_world,
		MusculoSkeletalSystem* mMusculoSkeletalSystem,
		BallInfo* ball_info,
		const std::vector<Eigen::VectorXd>& u0,
		int n,int max_iteration)
		:MusculoSkeletalDDP(rigid_world,soft_world,mMusculoSkeletalSystem,n,max_iteration),mBallInfo(ball_info),
		w_smooth_pos(0),w_smooth_vel(0),w_smooth2(0),w_smooth_ballpos(0),w_tracking(0),w_dtracking(0)
{
	param>>w_smooth_pos>>w_smooth_vel>>w_smooth2>>w_smooth_ballpos>>w_tracking>>w_dtracking;


	param.close();
	Eigen::VectorXd x0(mDofs*2);
	GetState(x0);
	mSoftX0 = mSoftWorld->GetPositions();

	initial_targets = u0;
	Eigen::VectorXd lower_bound(mSu);
	Eigen::VectorXd upper_bound(mSu);
	for(int i =0;i<mDofs;i++){
		// lower_bound[0] = -0.1;
		// upper_bound[0] = 1.0;

		lower_bound[0] = -0.05;
		upper_bound[0] = 0.8;

		// lower_bound[2] = 0.0;
		// upper_bound[2] = 0.0;
	}
	// for(int i =0;i<mDofs;i++){
	// 	lower_bound[i] = mMusculoSkeletalSystem->GetSkeleton()->getDof(i)->getPositionLowerLimit();
	// 	upper_bound[i] = mMusculoSkeletalSystem->GetSkeleton()->getDof(i)->getPositionUpperLimit();
	// 	lower_bound[i+mDofs] = -10.0;
	// 	upper_bound[i+mDofs] = 10.0;

	// 	lower_bound[i+2*mDofs] = 0;
	// 	upper_bound[i+2*mDofs] = 2000.0;
	// }

	Init(x0,u0,lower_bound,upper_bound);

	target_velocity = Eigen::Vector3d(0,4.0,0);
	ball_offset = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM();
}

double normP(const Eigen::VectorXd& v,double p)
{
	double ret = 0;
	for(int i=0;i<v.rows();i++)
		ret+=pow(v[i],p);
	return ret;
}
void
VelocityControlDDP::
EvalCf(const Eigen::VectorXd& x,double& cf)
{
	SetState(x);

	
	Eigen::Vector3d ball_p,ball_v;

	ball_p = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM();
	ball_v = mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOMLinearVelocity();
	if(cf ==123)
	std::cout<<"Velocity : "<<ball_v.transpose()<<std::endl;
	cf = 0.5*w_tracking*(ball_p - ball_offset).squaredNorm();
	// cf += 0.5*w_dtracking*(ball_v-Eigen::Vector3d(0,5.0,0)).squaredNorm();
	// cf += 0.5*w_dtracking*(ball_v-Eigen::Vector3d(0,5.0,0)).cwiseAbs().array().cube().sum();
	cf += 0.5*w_dtracking*normP(ball_v-target_velocity,2);
	// cf += 0.5*w_dtracking*(ball_v[1]-5.0)*(ball_v[1]-5.0);

}
void
VelocityControlDDP::
EvalCfx(const Eigen::VectorXd& x,Eigen::VectorXd& cfx)
{
	cfx.setZero();

	Eigen::VectorXd x_i;
	double delta = 0.001;
	double cf_minus,cf_plus;
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;
		cf_minus =0;
		cf_plus =0;
		x_i[i] = x[i] + delta;
		EvalCf(x_i,cf_plus);
		x_i[i] = x[i] - delta;
		EvalCf(x_i,cf_minus);

		cfx[i] = (cf_plus - cf_minus)/(2*delta);
	}
}
void
VelocityControlDDP::
EvalCfxx(const Eigen::VectorXd& x,Eigen::MatrixXd& cfxx)
{
	cfxx.resize(mSx,mSx);
	cfxx.setZero();
	Eigen::VectorXd x_i;
	double delta = 0.001;
	Eigen::VectorXd cfx_minus(mSx),cfx_plus(mSx);
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;
		x_i[i] = x[i] + delta;
		EvalCfx(x_i,cfx_plus);
		x_i[i] = x[i] - delta;
		EvalCfx(x_i,cfx_minus);

		cfxx.col(i) = (cfx_plus - cfx_minus)/(2*delta);
	}
}
void
VelocityControlDDP::
EvalC( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,double& c)
{
	// if(t>1)
		// c = 0.5*w_smooth_pos*(u-mu[t-1]).squaredNorm();
	// else
	c = 0.5*w_smooth_pos*(u-initial_targets[t]).squaredNorm();

	// if(t>1)
	// 	c = 0.5*w_smooth_pos*(u.head(mDofs)-2*mu[t-1].head(mDofs)+mu[t-2].head(mDofs)).squaredNorm();
	// else if(t==1)
	// 	c = 0.5*w_smooth_pos*(u.head(mDofs)-mu[t-1].head(mDofs)).squaredNorm();
	// else
	// 	c = 0.5*w_smooth_pos*(u.head(mDofs)-initial_targets[t].head(mDofs)).squaredNorm();

	// if(t>1)
	// 	c += 0.5*w_smooth_vel*(u.block(mDofs,0,mDofs,1)-2*mu[t-1].block(mDofs,0,mDofs,1)+mu[t-2].block(mDofs,0,mDofs,1)).squaredNorm();
	// else if(t==1)
	// 	c += 0.5*w_smooth_vel*(u.block(mDofs,0,mDofs,1)-mu[t-1].block(mDofs,0,mDofs,1)).squaredNorm();
	// else
	// 	c += 0.5*w_smooth_vel*(u.block(mDofs,0,mDofs,1)-initial_targets[t].block(mDofs,0,mDofs,1)).squaredNorm();

	// if(t>1)
	// 	c += 0.5*w_smooth2*(u.tail(mDofs)-2*mu[t-1].tail(mDofs)+mu[t-2].tail(mDofs)).squaredNorm();
	// else if(t==1)
	// 	c += 0.5*w_smooth2*(u.tail(mDofs)-mu[t-1].tail(mDofs)).squaredNorm();
	// else
	// 	c += 0.5*w_smooth2*(u.tail(mDofs)-initial_targets[t].tail(mDofs)).squaredNorm();

	// SetState(x);
	// c += 0.5*w_smooth_ballpos*normP(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOMLinearVelocity()-target_velocity,2);

		// (mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM()[0]-ball_offset[0])*
		// (mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR")->getCOM()[0]-ball_offset[0]);
	// c += 0.5*w_smooth2*u.tail(mDofs).squaredNorm();
}
void
VelocityControlDDP::
EvalCx( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cx)
{
	cx.setZero();

	Eigen::VectorXd x_i;
	double delta = 0.01;
	double c_minus,c_plus;
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;

		x_i[i] = x[i] + delta;
		EvalC(x_i,u,t,c_plus);
		x_i[i] = x[i] - delta;
		EvalC(x_i,u,t,c_minus);

		cx[i] = (c_plus - c_minus)/(2*delta);
	}
}
void
VelocityControlDDP::
EvalCu( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cu)
{
	cu.setZero();
	Eigen::VectorXd u_i;
	double delta = 0.01;
	double c_minus,c_plus;
	for(int i = 0;i<mSu;i++)
	{
		u_i = u;

		u_i[i] = u[i] + delta;
		EvalC(x,u_i,t,c_plus);
		u_i[i] = u[i] - delta;
		EvalC(x,u_i,t,c_minus);

		cu[i] = (c_plus - c_minus)/(2*delta);
	}

	// cu = w_effort*(u-initial_targets[t]);

	// if(t>1)
	// 	cu.head(mDofs) = w_smooth_pos*(u.head(mDofs)-2*mu[t-1].head(mDofs)+mu[t-2].head(mDofs));
	// else if(t==1)
	// 	cu.head(mDofs) = w_smooth_pos*(u.head(mDofs)-mu[t-1].head(mDofs));
	// else
	// 	cu.head(mDofs) = w_smooth_pos*(u.head(mDofs)-initial_targets[t].head(mDofs));

	// if(t>1)
	// 	cu.block(mDofs,0,mDofs,1) = w_smooth_vel*(u.block(mDofs,0,mDofs,1)-2*mu[t-1].block(mDofs,0,mDofs,1)+mu[t-2].block(mDofs,0,mDofs,1));
	// else if(t==1)
	// 	cu.block(mDofs,0,mDofs,1) = w_smooth_vel*(u.block(mDofs,0,mDofs,1)-mu[t-1].block(mDofs,0,mDofs,1));
	// else
	// 	cu.block(mDofs,0,mDofs,1) = w_smooth_vel*(u.block(mDofs,0,mDofs,1)-initial_targets[t].block(mDofs,0,mDofs,1));

	// if(t>1)
	// 	cu.tail(mDofs) = w_smooth2*(u.tail(mDofs)-2*mu[t-1].tail(mDofs)+mu[t-2].tail(mDofs));
	// else if(t==1)
	// 	cu.tail(mDofs) = w_smooth2*(u.tail(mDofs)-mu[t-1].tail(mDofs));
	// else
	// 	cu.tail(mDofs) = w_smooth2*(u.tail(mDofs)-initial_targets[t].tail(mDofs));

	// cu.tail(mDofs) = 0.5*w_smooth2*u.tail(mDofs);
}
void
VelocityControlDDP::
EvalCxx(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxx)
{
	cxx.resize(mSx,mSx);
	cxx.setZero();
	Eigen::VectorXd x_i;
	double delta = 0.01;
	Eigen::VectorXd cx_minus(mSx),cx_plus(mSx);
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;
		x_i[i] = x[i] + delta;
		EvalCx(x_i,u,t,cx_plus);
		x_i[i] = x[i] - delta;
		EvalCx(x_i,u,t,cx_minus);

		cxx.col(i) = (cx_plus - cx_minus)/(2*delta);
	}

}
void
VelocityControlDDP::
EvalCxu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxu)
{
	cxu.setZero();
	cxu.resize(mSx,mSu);
	cxu.setZero();
	Eigen::VectorXd u_i;
	double delta = 0.01;
	Eigen::VectorXd cx_minus(mSx),cx_plus(mSx);
	for(int i = 0;i<mSu;i++)
	{
		u_i = u;
		u_i[i] = u[i] + delta;
		EvalCx(x,u_i,t,cx_plus);
		u_i[i] = u[i] - delta;
		EvalCx(x,u_i,t,cx_minus);

		cxu.col(i) = (cx_plus - cx_minus)/(2*delta);
	}
}
void
VelocityControlDDP::
EvalCuu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cuu)
{
	cuu.setZero();
	Eigen::VectorXd u_i;
	double delta = 0.01;
	Eigen::VectorXd cu_minus(mSu),cu_plus(mSu);
	for(int i = 0;i<mSu;i++)
	{
		u_i = u;
		u_i[i] = u[i] + delta;
		EvalCu(x,u_i,t,cu_plus);
		u_i[i] = u[i] - delta;
		EvalCu(x,u_i,t,cu_minus);

		cuu.col(i) = (cu_plus - cu_minus)/(2*delta);
	}

	// cuu = (w_effort)*Eigen::MatrixXd::Identity(u.rows(),u.rows());
	// if(t>1)
	// cuu = w_smooth*Eigen::MatrixXd::Identity(u.rows(),u.rows());
	// cuu.block(0,0,mDofs,mDofs) = w_smooth_pos*Eigen::MatrixXd::Identity(mDofs,mDofs);
	// cuu.block(mDofs,mDofs,mDofs,mDofs) = w_smooth_vel*Eigen::MatrixXd::Identity(mDofs,mDofs);
	// cuu.block(2*mDofs,2*mDofs,mDofs,mDofs) = w_smooth2*Eigen::MatrixXd::Identity(mDofs,mDofs);
	// std::cout<<cuu<<std::endl;
}
void
VelocityControlDDP::
Evalf(  const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& f)
{
	if(t==0)
		mSoftWorld->SetPositions(mSoftX0);

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
	double delta = 0.001;
	for(int i = 0;i<mSx;i++)
	{
		x_i = x;

		x_i[i] = x[i] + delta;
		SetState(x_i);
		SetControl(u);
		Step();
		GetState(fx_i_plus);

		x_i[i] = x[i] - delta;
		SetState(x_i);
		SetControl(u);
		Step();
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
	double delta = 0.001;
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