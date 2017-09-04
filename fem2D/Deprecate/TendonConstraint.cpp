// #include "TendonConstraint.h"
// #include <iostream>
// const double TendonConstraint::f_toe = 0.33;
// const double TendonConstraint::k_toe = 3.0;
// const double TendonConstraint::e_toe = 0.02;
// const double TendonConstraint::e_t0  = 0.033;
// const double TendonConstraint::k_lin = 32.8407;
// const double TendonConstraint::k_0 = 32.8407;
// // const double TendonConstraint::k_0 = 10.0;
// double
// TendonConstraint::
// E_t(const double& e_t)
// {
// 	double ret = 0;
// 	double e_t_abs = fabs(e_t);
// 	if (e_t_abs <= e_toe)
// 		ret = (k_lin-k_0)/(6.0*e_toe)*(e_t_abs*e_t_abs*e_t_abs) + 0.5*k_0*e_t_abs*e_t_abs;
// 	else
// 		ret = 0.5*k_lin*e_t_abs*e_t_abs + 0.5*(k_0-k_lin)*e_toe*e_t_abs + (k_lin-k_0)/6.0*e_toe*e_toe;

// 	ret = 0.5*k_lin*e_t*e_t;
// 	return ret;
// }
// double
// TendonConstraint::
// f_t(const double& e_t)
// {
// 	double ret = 0;

// 	double e_t_abs = fabs(e_t);
// 	if (e_t_abs <= e_toe)
// 		ret = (k_lin - k_0)/(2.0*e_toe)*e_t_abs*e_t_abs+k_0*e_t_abs;
// 	else
// 		ret = k_lin*e_t_abs + 0.5*(k_0-k_lin)*e_toe;

// 	ret = -ret;
// 	ret = -k_lin*e_t;
// 	return ret;
// }
// double
// TendonConstraint::
// df_t(const double& e_t)
// {
// 	double ret = 0;
// 	double e_t_abs = fabs(e_t);
// 	if (e_t_abs <= e_toe)
// 		ret = (k_lin-k_0)/(e_toe)*e_t_abs + k_0;
// 	else
// 		ret = k_lin;

// 	ret = -ret;
// 	ret = -k_lin;
// 	return ret;
// }
// void
// TendonConstraint::
// ComputeCache(const Eigen::VectorXd& x)
// {
// 	mCachel = 0;
// 	for(int i =0;i<mP.size()-1;i++)
// 		mCachel += (mP[i+1]-mP[i]).norm();
// 	mCachel += (x.block<2,1>(mi0*2,0)-mP[0]).norm();

// 	mCachee_t = (mCachel-ml0)/ml0;
	
// }

// TendonConstraint::
// TendonConstraint(const double& stiffness,const std::vector<Eigen::Vector2d>& p,const int& i0,const double& l0)
// 	:Constraint(stiffness),mP(p),mi0(i0),ml0(l0)
// {
// }
// double
// TendonConstraint::
// EvalPotentialEnergy(const Eigen::VectorXd& x)
// {
// 	ComputeCache(x);

// 	return mStiffness*ml0*E_t(mCachee_t);
// }
// void
// TendonConstraint::
// EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
// {
// 	ComputeCache(x);

// 	Eigen::Vector2d p_x = x.block<2,1>(mi0*2,0) - mP[0];
// 	Eigen::Vector2d g = mStiffness*f_t(mCachee_t)*(p_x.normalized());

// 	gradient.block<2,1>(mi0*2,0) -= g;
// }
// void
// TendonConstraint::
// EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
// {
// 	ComputeCache(x);

// 	Eigen::Vector2d p_x = x.block<2,1>(mi0*2,0) - mP[0];
// 	double l = p_x.norm();
	

// 	double d_f1 = f_t(mCachee_t)/l;
// 	double d_f2 = -f_t(mCachee_t)/(l*l*l) + df_t(mCachee_t)/(ml0*l*l);
// 	auto dxi = dx.block<2,1>(mi0*2,0);

// 	dg.block<2,1>(mi0*2,0) -= mStiffness*(d_f1*dxi + d_f2*(p_x.dot(dxi))*p_x);
// }
// Eigen::Vector2d
// TendonConstraint::
// GetForces(const Eigen::VectorXd& x)
// {
// 	ComputeCache(x);
// 	Eigen::Vector2d f;
// 	f.setZero();
// 	// if(mP.size()==1)
// 	// 	f = f_t(mCachee_t)*((mP[0]-x.block<2,1>(mi0*2,0)).normalized());
// 	// else
// 	// 	f = f_t(mCachee_t)*((mP[mP.size()-1]-mP[mP.size()-2]).normalized());
// 	return f;
// }
// void
// TendonConstraint::
// EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
// {
// 	std::cout<<"TendonConstraint not supported."<<std::endl;
// }
// void
// TendonConstraint::
// EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
// {
// 	std::cout<<"TendonConstraint not supported."<<std::endl;
// }
// void
// TendonConstraint::
// EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
// {
// 	std::cout<<"TendonConstraint not supported."<<std::endl;
// }
// int 
// TendonConstraint::
// GetNumHessianTriplets()
// {
// 	return 0;
// }
// int 
// TendonConstraint::
// GetDof()
// {
// 	return 1;
// }