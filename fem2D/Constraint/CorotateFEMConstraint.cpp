#include "CorotateFEMConstraint.h"
#include <Eigen/SVD>
#include <Eigen/LU>
#include <iostream>
using namespace FEM;

void SVD2x2(const Eigen::Matrix2d& F,Eigen::Matrix2d& U, Eigen::Matrix2d& D, Eigen::Matrix2d& V,Eigen::Matrix2d& R)
{
	double EE = (F(0,0)+F(1,1))*0.5;
	double FF = (F(0,0)-F(1,1))*0.5;
	double GG = (F(1,0)+F(0,1))*0.5;
	double HH = (F(1,0)-F(0,1))*0.5;
	double QQ = sqrt(EE*EE+HH*HH);
	double RR = sqrt(FF*FF+GG*GG);
	D(0,0) = QQ+RR;
	D(1,1) = QQ-RR;
	if(QQ-RR<0)
		exit(0);
	double alpha = atan2(GG,FF);
	double beta = atan2(HH,EE);
	double phi = (alpha+beta)*(-0.5);
	double theta = (alpha-beta)*(-0.5);

	double cphi = cos(phi);
	double sphi = sin(phi);
	double cthe = cos(theta);
	double sthe = sin(theta);

	U(0,0) = cphi;
	U(0,1) = sphi;
	U(1,0) = -sphi;
	U(1,1) = cphi;

	V(0,0) = cthe;
	V(0,1) = sthe;
	V(1,0) = -sthe;
	V(1,1) = cthe;

	R(0,0) = cphi*cthe+sphi*sthe;
	R(0,1) = sphi*cthe-cphi*sthe;
	R(1,1) = R(0,0);
	R(1,0) = -R(0,1);
}
bool
CorotateFEMConstraint::
ComputeDeformationGradient(const Eigen::VectorXd& x)
{
	FEMConstraint::ComputeDeformationGradient(x);
	ComputeSVD(mCacheF);
	return true;
}
void
CorotateFEMConstraint::
ComputeSVD(const Eigen::Matrix2d& F)
{
	Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	auto D = svd.singularValues();
	
	mCacheD(0,0) = D[0];
	mCacheD(1,1) = D[1];

	mCacheU = svd.matrixU();
	mCacheV = svd.matrixV();
	mCacheR = mCacheU*mCacheV.transpose();


	// SVD2x2(F,mCacheU,mCacheD,mCacheV,mCacheR);
	mCacheF = F;
}

void
CorotateFEMConstraint::
ComputedP(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dP)
{
	Eigen::Matrix2d dR;
	ComputedR(dF,dR);

	dP =	mMu*(dF - dR) +
			mLambda*(
				(dR.transpose()*mCacheF+mCacheR.transpose()*dF).trace()*mCacheR +
				(mCacheR.transpose()*mCacheF-Eigen::Matrix2d::Identity()).trace()*dR
					);
}
void
CorotateFEMConstraint::
ComputedR(const Eigen::Matrix2d& dF,Eigen::Matrix2d& dR)
{
	double d1 = mCacheD(0,0);
	double d2 = mCacheD(1,1);
	if(fabs(d1-d2)<1E-6){
		Eigen::Matrix2d off_diag_M,M;
		M = mCacheU.transpose()* dF * mCacheV;
		off_diag_M = M;
		
		off_diag_M(0,0) = 0.0;
		off_diag_M(1,1) = 0.0;

		off_diag_M *= (1.0/d1);
		dR = mCacheU*off_diag_M*(mCacheV.transpose());
	}
	else
	{
		Eigen::Matrix2d Ainv;
		
		Eigen::Vector2d uv_tilda,m1221;
		Eigen::Matrix2d U_tilda,V_tilda;

		Ainv << d2,d1,
			    d1,d2;
		Ainv *= (1.0/(d2*d2-d1*d1));

		U_tilda.setZero();
		V_tilda.setZero();

		Eigen::Matrix2d M = mCacheU.transpose()*dF*mCacheV;

		m1221[0] = M(0,1);
		m1221[1] = M(1,0);

		uv_tilda = Ainv*m1221;

		U_tilda(0,1) = uv_tilda[0];
		U_tilda(1,0) = -uv_tilda[0];

		V_tilda(0,1) = uv_tilda[1];
		V_tilda(1,0) = -uv_tilda[1];

		U_tilda = mCacheU*U_tilda;  //dU
		V_tilda = mCacheV*V_tilda;  //dV
		
		dR = (U_tilda*(mCacheV.transpose()) + mCacheU*(V_tilda.transpose()));
	}


}
CorotateFEMConstraint::
CorotateFEMConstraint(const double& stiffness,const double& poisson_ratio,int i0,int i1,int i2,double vol,const Eigen::Matrix2d& invDm)
	:FEMConstraint(stiffness,poisson_ratio,i0,i1,i2,vol,invDm),
	mCacheU(Eigen::Matrix2d::Zero()),
	mCacheV(Eigen::Matrix2d::Zero()),
	mCacheD(Eigen::Matrix2d::Zero()),
	mCacheR(Eigen::Matrix2d::Zero())
{
}
double
CorotateFEMConstraint::
EvalPotentialEnergy(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);

	double vol_preserve_sqrt = (mCacheD-Eigen::Matrix2d::Identity()).trace();
	return mVol*(0.5*mMu*((mCacheF - mCacheR).norm())+0.5*mLambda*vol_preserve_sqrt*vol_preserve_sqrt);
}
void
CorotateFEMConstraint::
EvalGradient(const Eigen::VectorXd& x, Eigen::VectorXd& gradient)
{
	auto p0 = x.block<2,1>(mi0*2,0);
	auto p1 = x.block<2,1>(mi1*2,0);
	auto p2 = x.block<2,1>(mi2*2,0);

	auto d01 = p1 - p0;
	auto d02 = p2 - p0;
	// std::cout<<d01[0]*d02[1] - d01[1]*d02[0]<<std::endl;
	// std::cout<<d01[0]*d02[1] - d01[1]*d02[0]<<std::endl;
	// std::cout<<mCacheF.determinant()<<std::endl;
	ComputeDeformationGradient(x);

	Eigen::Matrix2d P = mMu*(mCacheF - mCacheR) + mLambda*((mCacheR.transpose()*mCacheF-Eigen::Matrix2d::Identity()).trace())*mCacheR;

	P = mVol*P*(mInvDm.transpose());
	gradient.block<2,1>(mi0*2,0) += -(P.block<2,1>(0,0) + P.block<2,1>(0,1));
	gradient.block<2,1>(mi1*2,0) += P.block<2,1>(0,0);
	gradient.block<2,1>(mi2*2,0) += P.block<2,1>(0,1);
}
void
CorotateFEMConstraint::
EvalHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& dx, Eigen::VectorXd& dg)
{
	ComputeDeformationGradient(x);
	Eigen::Matrix2d dDs,dF,dP;
	Eigen::Vector2d dx0(dx.block<2,1>(mi0*2,0));
	dDs.block<2,1>(0,0) = dx.block<2,1>(mi1*2,0)-dx0;
	dDs.block<2,1>(0,1) = dx.block<2,1>(mi2*2,0)-dx0;
	
	dF = dDs*(mInvDm);
	ComputedP(dF,dP);

	dP = mVol * dP * (mInvDm.transpose());

	dg.block<2,1>(mi0*2,0) += -(dP.block<2,1>(0,0) + dP.block<2,1>(0,1));
	dg.block<2,1>(mi1*2,0) += dP.block<2,1>(0,0);
	dg.block<2,1>(mi2*2,0) += dP.block<2,1>(0,1);
}
void
CorotateFEMConstraint::
GetDVector(int index, const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.block<4,1>(2*index,0) = md;
}
void
CorotateFEMConstraint::
EvaluateDVector(const Eigen::VectorXd& x)
{
	ComputeDeformationGradient(x);
	double a1,a2,d1,d2;
	Eigen::Matrix2d D;
	a1 = mCacheD(0,0);
	a2 = mCacheD(1,1);

	D.setZero();

	for(int i=0;i<10;i++)
	{
		d1 = D(0,0);
		d2 = D(1,1);

		double S = (d1*d2 - a1*a2 + 1)/(a1*a1 + a2*a2 + d1*d1 + d2*d2 + 2*a1*d1 + 2*a2*d2);

		D(0,0) = S*(a2 + d2);
		D(1,1) = S*(a1 + d1);

	}

	auto R_star = (0.1*mCacheR + 0.9*mCacheU*(D+mCacheD)*(mCacheV.transpose()));
	md.block<2,1>(0,0) = R_star.block<2,1>(0,0);
	md.block<2,1>(2,0) = R_star.block<2,1>(0,1);
}
void
CorotateFEMConstraint::
EvaluateJMatrix(int index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	Eigen::MatrixXd Ai(2*2,2*3);
	double d11 = mInvDm(0,0);
	double d12 = mInvDm(0,1);
	double d21 = mInvDm(1,0);
	double d22 = mInvDm(1,1);

	Ai<<
		-d11-d21,0,d11,0,d21,0,
		0,-d11-d21,0,d11,0,d21,
		-d12-d22,0,d12,0,d22,0,
		0,-d12-d22,0,d12,0,d22;

	auto MuAiT = 10.0*mMu*mVol*Ai.transpose();

	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+0, MuAiT(2*0+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+1, MuAiT(2*0+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+2, MuAiT(2*0+0,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+0, 2*index+3, MuAiT(2*0+0,2*0+3)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+0, MuAiT(2*0+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+1, MuAiT(2*0+1,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+2, MuAiT(2*0+1,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi0+1, 2*index+3, MuAiT(2*0+1,2*0+3)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+0, MuAiT(2*1+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+1, MuAiT(2*1+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+2, MuAiT(2*1+0,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+0, 2*index+3, MuAiT(2*1+0,2*0+3)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+0, MuAiT(2*1+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+1, MuAiT(2*1+1,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+2, MuAiT(2*1+1,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi1+1, 2*index+3, MuAiT(2*1+1,2*0+3)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+0, MuAiT(2*2+0,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+1, MuAiT(2*2+0,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+2, MuAiT(2*2+0,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+0, 2*index+3, MuAiT(2*2+0,2*0+3)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+0, MuAiT(2*2+1,2*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+1, MuAiT(2*2+1,2*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+2, MuAiT(2*2+1,2*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(2*mi2+1, 2*index+3, MuAiT(2*2+1,2*0+3)));
}

void
CorotateFEMConstraint::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	Eigen::MatrixXd Ai(2*2,2*3);
	double d11 = mInvDm(0,0);
	double d12 = mInvDm(0,1);
	double d21 = mInvDm(1,0);
	double d22 = mInvDm(1,1);

	Ai<<
		-d11-d21,0,d11,0,d21,0,
		0,-d11-d21,0,d11,0,d21,
		-d12-d22,0,d12,0,d22,0,
		0,-d12-d22,0,d12,0,d22;

	auto MuAiTAi = 10.0*mMu*mVol*((Ai.transpose())*Ai);

	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi0+0,MuAiTAi(2*0+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi0+1,MuAiTAi(2*0+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi0+0,MuAiTAi(2*0+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi0+1,MuAiTAi(2*0+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi1+0,MuAiTAi(2*0+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi1+1,MuAiTAi(2*0+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi1+0,MuAiTAi(2*0+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi1+1,MuAiTAi(2*0+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi2+0,MuAiTAi(2*0+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+0,2*mi2+1,MuAiTAi(2*0+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi2+0,MuAiTAi(2*0+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi0+1,2*mi2+1,MuAiTAi(2*0+1, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi0+0,MuAiTAi(2*1+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi0+1,MuAiTAi(2*1+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi0+0,MuAiTAi(2*1+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi0+1,MuAiTAi(2*1+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi1+0,MuAiTAi(2*1+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi1+1,MuAiTAi(2*1+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi1+0,MuAiTAi(2*1+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi1+1,MuAiTAi(2*1+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi2+0,MuAiTAi(2*1+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+0,2*mi2+1,MuAiTAi(2*1+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi2+0,MuAiTAi(2*1+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi1+1,2*mi2+1,MuAiTAi(2*1+1, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi0+0,MuAiTAi(2*2+0, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi0+1,MuAiTAi(2*2+0, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi0+0,MuAiTAi(2*2+1, 2*0+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi0+1,MuAiTAi(2*2+1, 2*0+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi1+0,MuAiTAi(2*2+0, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi1+1,MuAiTAi(2*2+0, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi1+0,MuAiTAi(2*2+1, 2*1+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi1+1,MuAiTAi(2*2+1, 2*1+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi2+0,MuAiTAi(2*2+0, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+0,2*mi2+1,MuAiTAi(2*2+0, 2*2+1)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi2+0,MuAiTAi(2*2+1, 2*2+0)));
	L_triplets.push_back(Eigen::Triplet<double>(2*mi2+1,2*mi2+1,MuAiTAi(2*2+1, 2*2+1)));




}


ConstraintType
CorotateFEMConstraint::
GetType()
{
	return ConstraintType::COROTATE;
}
