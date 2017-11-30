#include "FEM2D_Interface.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "GUI/GL_function.h"
#include "MusculoSkeletalSystem.h"
using namespace FEM;
void
DrawConstraint(Constraint* c,const Eigen::VectorXd& X)
{
    ConstraintType type = c->GetType();
	if(type == ConstraintType::ATTACHMENT)
	{
		auto* cc = static_cast<AttachmentConstraint*>(c);
		DrawLines(cc->GetP(),X.block<2,1>(cc->GetI0()*2,0),Eigen::Vector3d(1,0,0));
	}
	// else if(type == ConstraintType::SPRING)
	// {
	// 	auto* cc = static_cast<SpringConstraint*>(c);
	// 	DrawLines(X.block<2,1>(cc->GetI0()*2,0),X.block<2,1>(cc->GetI1()*2,0),Eigen::Vector3d(0,0,0));
	// }
	else if(type == ConstraintType::COROTATE)
	{
		auto* cc = static_cast<CorotateFEMConstraint*>(c);

		DrawTriangle(X.block<2,1>(cc->GetI0()*2,0),X.block<2,1>(cc->GetI1()*2,0),X.block<2,1>(cc->GetI2()*2,0),Eigen::Vector3d(0.8,0.8,0.8));
		
		DrawLines(X.block<2,1>(cc->GetI0()*2,0),X.block<2,1>(cc->GetI1()*2,0),Eigen::Vector3d(0,0,0));
		DrawLines(X.block<2,1>(cc->GetI1()*2,0),X.block<2,1>(cc->GetI2()*2,0),Eigen::Vector3d(0,0,0));
		DrawLines(X.block<2,1>(cc->GetI2()*2,0),X.block<2,1>(cc->GetI0()*2,0),Eigen::Vector3d(0,0,0));
	}
	// else if(type == ConstraintType::NEO_HOOKEAN)
	// {

	// }
	// else if(type == ConstraintType::VENANT_KIRCHHOFF)
	// {

	// }
	else if(type == ConstraintType::HILL_TYPE_MUSCLE ||type == ConstraintType::LINEAR_MUSCLE)
	{
		auto* cc = static_cast<MuscleConstraint*>(c);
        const auto& act = cc->GetActivationLevel();
        
        DrawTriangle(X.block<2,1>(cc->GetI0()*2,0),X.block<2,1>(cc->GetI1()*2,0),X.block<2,1>(cc->GetI2()*2,0),Eigen::Vector3d(0.8,0.8-0.3*act,0.8-0.3*act));
        
        DrawLines(X.block<2,1>(cc->GetI0()*2,0),X.block<2,1>(cc->GetI1()*2,0),Eigen::Vector3d(0,0,0));
        DrawLines(X.block<2,1>(cc->GetI1()*2,0),X.block<2,1>(cc->GetI2()*2,0),Eigen::Vector3d(0,0,0));
        DrawLines(X.block<2,1>(cc->GetI2()*2,0),X.block<2,1>(cc->GetI0()*2,0),Eigen::Vector3d(0,0,0));



  //       const auto& act = cc->GetActivationLevel();
  //       const auto& dir = cc->GetFiberDirection();
		// const auto& p0 = X.block<2,1>(cc->GetI0()*2,0);
  //       const auto& p1 = X.block<2,1>(cc->GetI1()*2,0);
		// const auto& p2 = X.block<2,1>(cc->GetI2()*2,0);

		// Eigen::Vector2d cen = (1.0/3.0)*(p0+p1+p2);

		// auto proj_0 = (((p2-p1).dot(p0-p1))/((p2-p1).squaredNorm()))*(p2-p1) + p1;
		// auto proj_1 = (((p2-p0).dot(p1-p0))/((p2-p0).squaredNorm()))*(p2-p0) + p0;
		// auto proj_2 = (((p1-p0).dot(p2-p0))/((p1-p0).squaredNorm()))*(p1-p0) + p0;
	
		// double h0 = (p0-proj_0).norm();        
		// double h1 = (p1-proj_1).norm();        
		// double h2 = (p2-proj_2).norm();
		// double min_h = (h0<h1?h0:h1);
		// min_h = (min_h<h2?min_h:h2);
  //       // glLineWidth(3.0);
  //       Eigen::Vector3d(act,0.0,1.0-act);
  //       DrawArrow(cen,-0.2*min_h*dir);
  //       DrawArrow(cen,0.2*min_h*dir);
	}
	// else if(type == ConstraintType::TENDON)
	// {

	// }
}
void
DrawConstraints(std::vector<Constraint*>& cs,const Eigen::VectorXd& X)
{
	for(auto& c: cs)
		DrawConstraint(c,X);
}
void
DrawMuscle(Muscle* muscle,const Eigen::VectorXd& x)
{
    glDisable(GL_DEPTH_TEST);
    for(auto& c: muscle->constraints)
        DrawConstraint(c,x);
    for(auto& c: muscle->muscleConstraints)
        DrawConstraint(c,x);
    // DrawConstraint(muscle->origin,x);
    // DrawConstraint(muscle->insertion,x);


    // glPointSize(5.0);
    // Eigen::Vector3d(1,0,0);
    // DrawPoint(x.block<2,1>(muscle->origin->GetAttachedNode()*2,0));
    // Eigen::Vector3d(0,0,1);
    // DrawPoint(x.block<2,1>(muscle->insertion->GetAttachedNode()*2,0));
    // glPointSize(1.0);

    {
    ;
        std::vector<Eigen::Vector2d> p;
        for(auto& ap : muscle->originWayPoints)
            p.push_back(GetPoint(ap));

        for(int i = 0;i<p.size()-1;i++)
            DrawLines(p[i],p[i+1],Eigen::Vector3d(0,0,1));
        DrawLines(x.block<2,1>(muscle->origin->GetI0()*2,0),p[0],Eigen::Vector3d(0,0,1));

        if(muscle->force_origin.norm()>1E0)
            DrawArrow(p.back(),0.00001*muscle->force_origin,Eigen::Vector3d(0.7,0,0.5));
        glPointSize(2.0);
        for(auto& v : p)
            DrawPoint(v,Eigen::Vector3d(0,0,0));
        glPointSize(1.0);
    }

    {
        std::vector<Eigen::Vector2d> p;
        for(auto& ap : muscle->insertionWayPoints)
            p.push_back(GetPoint(ap));

        for(int i = 0;i<p.size()-1;i++)
            DrawLines(p[i],p[i+1],Eigen::Vector3d(0,0,1));
        DrawLines(x.block<2,1>(muscle->insertion->GetI0()*2,0),p[0],Eigen::Vector3d(0,0,1));
        
        if(muscle->force_insertion.norm()>1E0)
            DrawArrow(p.back(),0.00001*muscle->force_insertion,Eigen::Vector3d(0.7,0,0.5));

        glPointSize(2.0);
        for(auto& v : p)
            DrawPoint(v, Eigen::Vector3d(0,0,0));
        glPointSize(1.0);
    }

    
    
    

    glDisable(GL_DEPTH_TEST);
}