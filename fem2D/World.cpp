#include "World.h"
#include "Constraint/ConstraintHeaders.h"
#include "Mesh/MeshHeaders.h"
#include <iostream>
using namespace FEM;
World::
World(
		IntegrationMethod integration_method,
		double time_step,
		int max_iteration,
		double damping_coeff)
	:mIntegrationMethod(integration_method),
	mTimeStep(time_step),
	mMaxIteration(max_iteration),
	mDampingCoefficinent(damping_coeff),
	mFrame(0),
	mNumVertices(0),
	mConstraintDofs(0),
	mIsInitialized(false)
{
}

void
World::
Initialize()
{
	mTime = 0.0;
	mConstraintDofs = 0;

	for(auto c : mConstraints){
		mConstraintDofs += c->GetDof();
	}

	mV.resize(2*mNumVertices);
	mExternalForces.resize(2*mNumVertices);

	mV.setZero();
	mExternalForces.setZero();

	mMassMatrix.resize(2*mNumVertices,2*mNumVertices);
	mInvMassMatrix.resize(2*mNumVertices,2*mNumVertices);
	mIdentityMatrix.resize(2*mNumVertices,2*mNumVertices);

	std::vector<Eigen::Triplet<double>> i_triplets;
	std::vector<Eigen::Triplet<double>> m_triplets;
	std::vector<Eigen::Triplet<double>> inv_m_triplets;

	i_triplets.reserve(2*mNumVertices);
	m_triplets.reserve(2*mNumVertices);
	inv_m_triplets.reserve(2*mNumVertices);
	
	for(int i = 0;i<mNumVertices;i++)
	{
		m_triplets.push_back(Eigen::Triplet<double>((i)*2+0,(i)*2+0,mUnitMass[i]));
		m_triplets.push_back(Eigen::Triplet<double>((i)*2+1,(i)*2+1,mUnitMass[i]));

		inv_m_triplets.push_back(Eigen::Triplet<double>((i)*2+0,(i)*2+0,1.0/mUnitMass[i]));
		inv_m_triplets.push_back(Eigen::Triplet<double>((i)*2+1,(i)*2+1,1.0/mUnitMass[i]));

		i_triplets.push_back(Eigen::Triplet<double>((i)*2+0,(i)*2+0,1.0));
		i_triplets.push_back(Eigen::Triplet<double>((i)*2+1,(i)*2+1,1.0));
	}

	mMassMatrix.setFromTriplets(m_triplets.cbegin(), m_triplets.cend());
	mInvMassMatrix.setFromTriplets(inv_m_triplets.cbegin(), inv_m_triplets.cend());
	mIdentityMatrix.setFromTriplets(i_triplets.cbegin(), i_triplets.cend());

	mQn.resize(2*mNumVertices);
	mQn.setZero();

	if(mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS || mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)
		PreComputeProjectiveDynamics();

	mIsInitialized = true;
	std::cout<<"Total degree of freedom : "<<mX.rows()<<std::endl;
	std::cout<<"Total constraints : "<<mConstraints.size()<<std::endl;
}

void
World::
TimeStepping(bool isIntegrated)
{
	if(!mIsInitialized){
		std::cout<<"Engine not initialized."<<std::endl;
		return;
	}
	Eigen::VectorXd x_n1(mNumVertices*2);

	ComputeExternalForces();
	mQn = mX + mTimeStep*mV + (mTimeStep*mTimeStep)*(mInvMassMatrix*mExternalForces);

	switch(mIntegrationMethod)
	{
	case SYMPLECTIC_EXPLICIT:
		x_n1 = IntegrateSymplecticExplicit();
	break;
	case SEMI_IMPLICIT:
		x_n1 = IntegrateSemiImplicit();
	break;
	case IMPLICIT_NEWTON_METHOD:
		x_n1 = IntegrateNewtonMethod();
	break;
	case QUASI_STATIC:
		x_n1 = IntegrateQuasiStatic();
	break;
	case PROJECTIVE_DYNAMICS:
		x_n1 = IntegrateProjectiveDynamics();
	break;
	case PROJECTIVE_QUASI_STATIC:
		x_n1 = IntegrateProjectiveQuasiStatic();
	break;
	default:
	break;
	}
	for(auto& c : mConstraints)
	{
		ConstraintType type = c->GetType();
		if(type == ConstraintType::HILL_TYPE_MUSCLE)
		{
			auto* cc = static_cast<HillTypeMuscleConstraint*>(c);
			cc->SetPreviousL(mX);
		}
	}
	UpdatePositionsAndVelocities(x_n1);
	mV *= mDampingCoefficinent;
	if(isIntegrated)
	{	
		mTime += mTimeStep;
		mFrame++;
	}
}
Eigen::VectorXd
World::
ComputeJacobian(std::vector<MuscleConstraint*>& mc
	,std::vector<AttachmentConstraint*>& attachment_vector)
{
	
	//Analytic Jacobian	
	Eigen::VectorXd dg_da(mX.rows());
	dg_da.setZero();
	for(auto& c: mc)
		c->Evaldg_da(mX,dg_da);
	Eigen::VectorXd dx_da;
	if(mIntegrationMethod == QUASI_STATIC || mIntegrationMethod == IMPLICIT_NEWTON_METHOD)
		dx_da= -ConjugateGradient(dg_da,mX);
	else
		dx_da = -mQuasiStaticSolver.solve(dg_da);
	Eigen::VectorXd Xplusdx = mX + dx_da;

	Eigen::VectorXd Ji(mX.rows());
	Eigen::VectorXd ret(attachment_vector.size()*2);

	Ji.setZero();

	for(int i =0;i<attachment_vector.size();i++)	
		attachment_vector[i]->EvalGradient(mX,Ji);
		
	Ji = -Ji;
	
	for(int i =0;i<attachment_vector.size();i++)	
		attachment_vector[i]->EvalGradient(Xplusdx,Ji);

	Ji += dg_da;

	for(int i =0;i<attachment_vector.size();i++)
		ret.block<2,1>(i*2,0) = -Ji.block<2,1>(attachment_vector[i]->GetI0()*2,0);
	
	//Numerical Jacobian
	// for(double delta =0.001;delta<1.0;delta*=2)
	// {
		// double delta = 0.01;
		// double save_activation = mc[0]->GetActivationLevel();
		// Eigen::VectorXd save_X = mX;
		// Eigen::VectorXd g_0(mX.rows()),g_1(mX.rows());
		// g_0.setZero();
		// g_1.setZero();
		
		// double a_minus_dx = save_activation-delta;
		// double a_plus_dx = save_activation+delta;
		// if(a_plus_dx>1.0)
		// 	a_plus_dx =1.0;
		// if(a_minus_dx<0.0)
		// 	a_minus_dx =0.0;


		// for(auto& c : mc)
		// 	c->SetActivationLevel(a_minus_dx);

		// if(mIntegrationMethod == QUASI_STATIC)
		// 	mX = IntegrateQuasiStatic();
		// else
		// 	mX = IntegrateProjectiveQuasiStatic();

		// for(int i =0;i<attachment_vector.size();i++)	
		// 	attachment_vector[i]->EvalGradient(mX,g_0);

		// for(auto& c : mc)
		// 	c->SetActivationLevel(a_plus_dx);

		// if(mIntegrationMethod == QUASI_STATIC)
		// 	mX = IntegrateQuasiStatic();
		// else
		// 	mX = IntegrateProjectiveQuasiStatic();

		// for(int i =0;i<attachment_vector.size();i++)	
		// 	attachment_vector[i]->EvalGradient(mX,g_1);

		// for(auto& c : mc)
		// 	c->SetActivationLevel(save_activation);

		// mX = save_X;
		// Eigen::VectorXd Ji = (g_1 - g_0)*(1.0/(a_plus_dx-a_minus_dx));
		// Eigen::VectorXd ret(attachment_vector.size()*2);
		// for(int i =0;i<attachment_vector.size();i++)
		// 	ret.block<2,1>(i*2,0) = Ji.block<2,1>(attachment_vector[i]->GetI0()*2,0);
	
		// std::cout<<"Numerical (da : "<<a_plus_dx-a_minus_dx<<") :"<<ret.transpose()<<std::endl;
	// }
	// exit(0);
	return ret;
}
void
World::
AddBody(const Eigen::VectorXd& x0,const std::vector<Constraint*>& c,const double& mass)
{
	int nv = ((x0.rows())/2);
	mNumVertices += nv;

	auto temp_X(mX);
	mX.resize(mNumVertices*2);
	mX.head(temp_X.rows()) = temp_X;
	mX.tail(x0.rows()) = x0;
	
	mConstraints.insert(mConstraints.end(), c.begin(), c.end());
	double unit_mass = mass/((double)nv);
	for(int i=0;i<nv;i++)
		mUnitMass.push_back(unit_mass);

	if(mIsInitialized)
		Initialize();
}
void
World::
AddConstraint(Constraint* c)
{
	mConstraints.push_back(c);
	if((mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS||mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)&& mIsInitialized)
	{
		mConstraintDofs = 0;
		for(auto c : mConstraints){
			mConstraintDofs += c->GetDof();
		}
		PreComputeProjectiveDynamics();
	}
}

void
World::
RemoveConstraint(Constraint* c)
{
	bool isRemoved = false;
	for(int i = 0;i<mConstraints.size();i++)
	{
		if(mConstraints[i]==c)
		{
			mConstraints.erase(mConstraints.begin() + i);
			isRemoved = true;
			break;
		}
	}
	if(isRemoved)
		if((mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS||mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)&& mIsInitialized)
		{
			mConstraintDofs = 0;
			for(auto c : mConstraints){
				mConstraintDofs += c->GetDof();
			}
			PreComputeProjectiveDynamics();
		}
	}

int
World::
GetClosestNode(const Eigen::Vector2d& x)
{
	double min_distance = 1E6;
	int ret = -1;
	for(int i=0;i<mNumVertices;i++){
		double distance = (mX.block<2,1>(i*2,0)-x).norm();
		if(min_distance>distance){
			ret = i;
			min_distance=distance;
		}
	}


	if(min_distance>1E-1)
		return -1;
	else
		return ret;
}

double 								
World::
GetTimeStep()
{
	return mTimeStep;
}

double 								
World::
GetTime()
{
	return mTime;
}
void
World::
SetTime(double t)
{
	mTime = t;
}
int    								
World::
GetNumVertices()
{
	return mNumVertices;
}

const Eigen::VectorXd& 				
World::
GetPositions()
{
	return mX;
}
void
World::
SetPositions(const Eigen::VectorXd& p)
{
	mX = p;
}
std::vector<Constraint*>&			
World::
GetConstraints()
{
	return mConstraints;
}










Eigen::VectorXd
World::
IntegrateSymplecticExplicit()
{
	Eigen::VectorXd x_n1(2*mNumVertices);
	Eigen::VectorXd gc(2*mNumVertices);

	EvalConstraintsGradient(mX,gc);

	x_n1 = mX + mTimeStep*mV + (mTimeStep*mTimeStep)*mInvMassMatrix*(-gc + mExternalForces);

	return x_n1;
}
Eigen::VectorXd
World::
IntegrateSemiImplicit()
{
	Eigen::VectorXd v_n1(2*mNumVertices);
	Eigen::VectorXd gc(2*mNumVertices);

	EvalConstraintsGradient(mX,gc);

	Eigen::VectorXd b(2*mNumVertices);

	b = (1.0/(mTimeStep*mTimeStep))*mMassMatrix*mV + (1.0/mTimeStep)*(-gc + mExternalForces);

	v_n1 = ConjugateGradient(b,mQn);

	return mX + mTimeStep*v_n1;
}

Eigen::VectorXd
World::
IntegrateNewtonMethod()
{
	Eigen::VectorXd x_n1(2*mNumVertices);

	x_n1 = mQn;
	int i=0;
	Eigen::VectorXd g(2*mNumVertices);
	g.setZero();
	for(;i<mMaxIteration;i++)
	{
		InversionFree(x_n1);
		EvalGradient(x_n1,g);
		if(g.squaredNorm() < 1E-2){
			break;
		}

		Eigen::VectorXd dir  = -ConjugateGradient(g,x_n1);
		// dir.normalize();
		// double alpha = ComputeStepSize(x_n1,g,dir);
		x_n1 = x_n1 + 0.5*dir;
		// std::cout<<"Iteration ("<<i<<")"<<std::endl;
		// std::cout<<"Force norm : "<<g.norm()<<std::endl;
		// std::cout<<"Alpha : "<<alpha<<std::endl<<std::endl;
	}
	// std::cout<<"force norm : "<<g.norm()<<std::endl<<std::endl;
	return x_n1;
}
Eigen::VectorXd
World::
IntegrateQuasiStatic()
{
	Eigen::VectorXd x_n1(2*mNumVertices);
	x_n1 = mX;
	double alpha=0;
	Eigen::VectorXd g;
	int i=0;
	for(;i<mMaxIteration;i++)
	{
		InversionFree(x_n1);
		EvalConstraintsGradient(x_n1,g);

		if(g.squaredNorm() < 1E-2){
			// std::cout<<" g: "<<g.norm()<<std::endl;
			break;
		}


		Eigen::VectorXd dir  = -ConjugateGradient(g,x_n1);
		// alpha = ComputeStepSize(x_n1,g,dir);
		x_n1 = x_n1 + 0.5*dir;
		// std::cout<<"Iteration ("<<i<<")"<<std::endl;
		// std::cout<<"Force norm : "<<g.squaredNorm()<<std::endl;
		// std::cout<<std::endl;

	}
	

	// std::cout<<"Iteration ("<<i<<")"<<std::endl;
	// std::cout<<"Force norm : "<<g.squaredNorm()<<std::endl;
	return x_n1;	
}


Eigen::VectorXd
World::
IntegrateProjectiveDynamics()
{
	Eigen::VectorXd x_n1(2*mNumVertices);
	Eigen::VectorXd b(2*mNumVertices);
	Eigen::VectorXd d(2*mConstraintDofs);
	b = (1.0/(mTimeStep*mTimeStep))*mMassMatrix*mQn;

	x_n1 = mQn;
	int i=0;
	for(;i<mMaxIteration;i++)
	{
		InversionFree(x_n1);
		EvaluateDVector(x_n1,d);
		x_n1 = mDynamicSolver.solve(b+mJ*d);
	}
	// std::cout<<"Projective : "<<i<<std::endl<<std::endl;
	return x_n1;
}
Eigen::VectorXd
World::
IntegrateProjectiveQuasiStatic()
{
	Eigen::VectorXd x_n1(2*mNumVertices);
	Eigen::VectorXd d(2*mConstraintDofs);

	x_n1 = mX;
	int i=0;
	for(;i<mMaxIteration;i++)
	{
		InversionFree(x_n1);
		EvaluateDVector(x_n1,d);
		x_n1 = mQuasiStaticSolver.solve(mJ*d);
	}
	return x_n1;
}

void
World::
FactorizeLLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& lltSolver)
{
	Eigen::SparseMatrix<double> A_prime = A;
	lltSolver.analyzePattern(A_prime);
	lltSolver.factorize(A_prime);
	double damping = 1E-6;
	bool success = true;
	while (lltSolver.info() != Eigen::Success)
	{
	    damping *= 10;
	    A_prime = A + damping*mIdentityMatrix;
	    lltSolver.factorize(A_prime);
	    success = false;
	}
	if (!success)
	    std::cout << "factorize failure (damping : " << damping<<" )"<<std::endl;
}

void
World::
FactorizeLDLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldltSolver)
{
	Eigen::SparseMatrix<double> A_prime = A;
	ldltSolver.analyzePattern(A_prime);
	ldltSolver.factorize(A_prime);
	double damping = 1E-6;
	bool success = true;
	while (ldltSolver.info() != Eigen::Success)
	{
	    damping *= 10;
	    A_prime = A + damping*mIdentityMatrix;
	    ldltSolver.factorize(A_prime);
	    success = false;
	}
	if (!success)
	    std::cout << "factorize failure (damping : " << damping<<" )"<<std::endl;
}
Eigen::VectorXd
World::
ConjugateGradient(const Eigen::VectorXd& b,const Eigen::VectorXd& x0)
{
	//Solve Ax = b
	Eigen::VectorXd r(2*mNumVertices);
	Eigen::VectorXd p(2*mNumVertices);
	Eigen::VectorXd x(2*mNumVertices);
	Eigen::VectorXd Ax(2*mNumVertices);

	
	if(mIntegrationMethod == QUASI_STATIC)
		EvalConstraintsHessian(x0,x0,Ax);
	else
		EvalHessian(x0,x0,Ax);

	x = x0;
	r = b - Ax;

	p = r;
	double rs_old = r.dot(r);
	int i=0;
	for(;i<10000;i++)
	{
		if(mIntegrationMethod == QUASI_STATIC)
			EvalConstraintsHessian(x0,p,Ax);
		else
			EvalHessian(x0,p,Ax);

		double alpha = rs_old /(p.dot(Ax));
		x = x + alpha * p;
		r = r - alpha * Ax;
		double rs_new = r.dot(r);

		if(rs_new<1E-6){
			rs_old = rs_new;
			break;
		}
		if(rs_old<rs_new)
			break;
		p = r + (rs_new/rs_old)* p;
		rs_old = rs_new;
	}
	// std::cout<<"CG Residual : "<<rs_old<<" iteration : "<<i<<std::endl;
	return x;

}

double
World::
ComputeStepSize(const Eigen::VectorXd& x,const Eigen::VectorXd& gradient,const Eigen::VectorXd& descent)
{
	double alpha = 1.0;
	double c1 = 0.03;

	double current_val;
	double next_step_val;
	Eigen::VectorXd x_next;
	if(mIntegrationMethod == QUASI_STATIC)
		current_val = EvalConstraintsObjectives(x);
	else
		current_val = EvalObjectives(x);

	double g_dot_d = gradient.dot(descent);
	
	for(int i=0;i<32;i++)
	{
		x_next = x + alpha*descent;
		
		if(mIntegrationMethod == QUASI_STATIC)
			next_step_val = EvalConstraintsObjectives(x_next);
		else
			next_step_val = EvalObjectives(x_next);
		if((current_val + alpha*c1*g_dot_d>= next_step_val)){
			break;
		}
		alpha *= 0.5;
	}

	// if(alpha<1E-5)
		// alpha = 0.0;

	// std::cout << "Linesearch Stepsize = " << alpha << std::endl;
	// std::cout << "current energy  = " << next_step_val << std::endl;
	// std::cout << "previous energy = " << current_val << std::endl;
	// std::cout << "rhs (previous energy + alpha * t * gradient.dot(descet_dir)) = " << current_val + alpha*c1*g_dot_d << std::endl;
	return alpha;
}

void
World::
PreComputeProjectiveDynamics()
{
	EvaluateJMatrix(mJ);
	EvaluateLMatrix(mL);
	auto H2ML = (1.0/(mTimeStep*mTimeStep))*mMassMatrix+mL;
	FactorizeLLT(H2ML,mDynamicSolver);
	FactorizeLLT(mL,mQuasiStaticSolver);
}

void
World::
InversionFree(Eigen::VectorXd& x)
{	
	for(auto& c: mConstraints)
	{
		if(c->GetType()==ConstraintType::COROTATE)
		{
			auto* cc = static_cast<CorotateFEMConstraint*>(c);
			cc->AddInversionFreePosition(x);
		}
	}
	
}
double
World::
EvalObjectives(const Eigen::VectorXd& x)
{
	auto x_q = x - mQn;
	double val = EvalConstraintsObjectives(x);
	val += 0.5*(1.0/(mTimeStep*mTimeStep))*(x_q.dot(mMassMatrix*x_q));

	return val;
}

void
World::
EvalGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g)
{
	EvalConstraintsGradient(x,g);
	g += (1.0/(mTimeStep*mTimeStep))*mMassMatrix*(x - mQn);	
}

void
World::
EvalHessian(const Eigen::VectorXd& x,const Eigen::VectorXd& dx,Eigen::VectorXd& dg)
{
	EvalConstraintsHessian(x,dx,dg);
	dg += (1.0/(mTimeStep*mTimeStep))*mMassMatrix*dx;	
}
#pragma omp declare reduction (+: Eigen::VectorXd: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Eigen::VectorXd::Zero(omp_orig.size()))

double 
World::
EvalConstraintsObjectives(const Eigen::VectorXd& x)
{
	double val = 0;

#pragma omp parallel for reduction(+:val)
	for(int i =0;i<mConstraints.size();i++)
		val += mConstraints[i]->EvalPotentialEnergy(x);

	return val;
}

void
World::
EvalConstraintsGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g)
{
	g.resize(mNumVertices*2);
	g.setZero();
#pragma omp parallel for reduction(+:g)
	for(int i =0;i<mConstraints.size();i++){
		mConstraints[i]->EvalGradient(x,g);

	}
}

void
World::
EvalConstraintsHessian(const Eigen::VectorXd& x,const Eigen::VectorXd& dx,Eigen::VectorXd& dg)
{
	dg.resize(2*mNumVertices);
	dg.setZero();
#pragma omp parallel for reduction(+:dg)
	for(int i =0;i<mConstraints.size();i++)
		mConstraints[i]->EvalHessian(x,dx,dg);
}
void
World::
EvaluateDVector(const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.resize(2*mConstraintDofs);

	int index = 0;
#pragma omp parallel for
	for(int i=0;i<mConstraints.size();i++)
	{
		mConstraints[i]->EvaluateDVector(x);
	}

	for(int i =0;i<mConstraints.size();i++)
	{
		mConstraints[i]->GetDVector(index,x,d);
		index+=mConstraints[i]->GetDof();
	}
}
void
World::
EvaluateJMatrix(Eigen::SparseMatrix<double>& J)
{
	J.resize(2*mNumVertices,2*mConstraintDofs);
	std::vector<Eigen::Triplet<double>> J_triplets;

	int index = 0;
	for(int i =0;i<mConstraints.size();i++)
	{
		mConstraints[i]->EvaluateJMatrix(index,J_triplets);
		index+=mConstraints[i]->GetDof();

	}

	J.setFromTriplets(J_triplets.cbegin(), J_triplets.cend());
}
void
World::
EvaluateLMatrix(Eigen::SparseMatrix<double>& L)
{
	L.resize(2*mNumVertices,2*mNumVertices);
	std::vector<Eigen::Triplet<double>> l_triplets;

	for(auto c : mConstraints)
		c->EvaluateLMatrix(l_triplets);

	L.setFromTriplets(l_triplets.cbegin(), l_triplets.cend());
}
void
World::
ComputeExternalForces()
{
	mExternalForces.setZero();

	//Add Gravity forces
	for(int i=0;i<mNumVertices;i++)
		// mExternalForces[2*i+1] = -9.81;
		mExternalForces[2*i+1] = 0.0;

	mExternalForces = mMassMatrix * mExternalForces;

	for(int i=0;i<mNumVertices;i++)
		if(mX[i*2+1]<-8)
		{
			mExternalForces[2*i+1] += -100*(mX[i*2+1]+8);
		}
}
void
World::
UpdatePositionsAndVelocities(const Eigen::VectorXd& x_n1)
{
	mV = (x_n1 - mX)*(1.0/mTimeStep);
	mX = x_n1;
}
