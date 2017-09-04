// #ifndef __TENDON_CONSTRAINT_H__
// #define __TENDON_CONSTRAINT_H__
// #include "Constraint.h"
// #include <vector>
// class TendonConstraint : public Constraint
// {
// protected:
// 	std::vector<Eigen::Vector2d>  mP;
// 	int mi0;
// 	double ml0;

// 	double mCachel;
// 	double mCachee_t;
// 	static double E_t(const double& e_t); 		// E(e_t)
// 	static double f_t(const double& e_t); 		// f(e_t)
// 	static double df_t(const double& e_t);		// df(e_t)

// 	static const double f_toe;
// 	static const double k_toe;
// 	static const double e_toe;
// 	static const double e_t0;
// 	static const double k_lin;
// 	static const double k_0;
// 	void ComputeCache(const Eigen::VectorXd& x);
// public:
// 	TendonConstraint(const double& stiffness,const std::vector<Eigen::Vector2d>& p,const int& i0,const double& l0);
// 	virtual double  EvalPotentialEnergy(const Eigen::VectorXd& x);
// 	virtual void EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient);
// 	virtual void EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg);

// 	virtual void EvaluateDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d);
// 	virtual void EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets);
// 	virtual void EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets);

// 	virtual int GetNumHessianTriplets();
// 	virtual int GetDof();

// 	std::vector<Eigen::Vector2d>& GetTendonPoints() { return mP; };
// 	const int& GetTendonIndex() { return mi0; };
// 	Eigen::Vector2d GetForces(const Eigen::VectorXd& x);

// 	virtual ConstraintType GetType()	   {return ConstraintType::TENDON; };
// 	virtual void AddOffset(const int& offset){
// 		mi0 +=offset;
// 	};
// };


// #endif