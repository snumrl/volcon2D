#ifndef __VELOCITY_CONTROL_DDP_H__
#define __VELOCITY_CONTROL_DDP_H__
#include "MusculoSkeletalDDP.h"
class BallInfo;

class VelocityControlDDP : public MusculoSkeletalDDP
{
public:
	VelocityControlDDP(const dart::simulation::WorldPtr& rigid_world,
		FEM::World* soft_world,
		MusculoSkeletalSystem* mMusculoSkeletalSystem,
		BallInfo* ball_info,
		const std::vector<Eigen::VectorXd>& u0,
		int n,int max_iteration);
public:
	void EvalCf(const Eigen::VectorXd& x,double& cf) override;
	void EvalCfx(const Eigen::VectorXd& x,Eigen::VectorXd& cfx) override;
	void EvalCfxx(const Eigen::VectorXd& x,Eigen::MatrixXd& cfxx) override;

	void EvalC( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,double& c) override;
	void EvalCx( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cx) override;
	void EvalCu( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& cu) override;
	void EvalCxx(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxx) override;
	void EvalCxu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cxu) override;
	void EvalCuu(const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& cuu) override;

	void Evalf(  const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::VectorXd& f) override;
	void Evalfx( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& fx) override;
	void Evalfu( const Eigen::VectorXd& x,const Eigen::VectorXd& u,int t,Eigen::MatrixXd& fu) override;

protected:
	BallInfo*	 mBallInfo;

	
	std::vector<Eigen::VectorXd> initial_targets;
	Eigen::VectorXd mSoftX0;
	double w_smooth_pos,w_smooth_vel,w_smooth2,w_smooth_ballpos,w_tracking,w_dtracking;
	Eigen::Vector3d target_velocity;
	Eigen::Vector3d ball_offset;
};

#endif