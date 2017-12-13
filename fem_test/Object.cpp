#include "Object.h"
#include <vector>
#include <Eigen/Core>
#include <tinyxml.h>

using namespace FEM;

Object::
Object(const Eigen::Vector2d& trans)
	:mMuscleStiffness(1E6),mYoungsModulus(1E6),mPoissonRatio(0.3),mMesh()
{
	Eigen::Affine2d T;
	T.setIdentity();
	T.translation() = trans;
	// mMesh = new RectangueMesh(0.2,0.1,8,4,T);
	mMesh = new DiamondMesh(0.2,0.2,8,8,T);
}

void
Object::
AddMuscle(
	const std::vector<Eigen::Vector3i>& indexList,
	const Eigen::Vector2d& fiber_direction)
{
	mMuscles.push_back(new Muscle());
	auto muscle = mMuscles.back();
	const auto& vertices = mMesh->GetVertices();
	const auto& triangles = mMesh->GetTriangles();
	std::vector<FEM::Constraint*> constraints;
	Eigen::VectorXd v(vertices.size()*2);

	for(const auto& idx : indexList)
	{
		int i0,i1,i2;
		Eigen::Vector2d p0,p1,p2;
		
		i0 = idx[0];
		i1 = idx[2];
		i2 = idx[1];
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		Eigen::Matrix2d Dm;

		Dm.block<2,1>(0,0) = p1 - p0;
		Dm.block<2,1>(0,1) = p2 - p0;
		if(Dm.determinant()<0)
		{
			i0 = idx[0];
			i1 = idx[1];
			i2 = idx[2];

			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			Dm.block<2,1>(0,0) = p1 - p0;
			Dm.block<2,1>(0,1) = p2 - p0;

			muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
			// muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		}
	}

	// Add Attachment 
	

	
	muscle->activationLevel = 0.0;
}
void
Object::
Initialize(FEM::World* world)
{
	double offset = world->GetNumVertices();
	mMuscles.push_back(new Muscle());
	auto muscle = mMuscles.back();
	const auto& vertices = mMesh->GetVertices();
	const auto& triangles = mMesh->GetTriangles();

	for(const auto& tri : triangles)
	{
		int i0,i1,i2;
		Eigen::Vector2d p0,p1,p2;
		
		i0 = tri[0];
		i1 = tri[2];
		i2 = tri[1];
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		Eigen::Matrix2d Dm;

		Dm.block<2,1>(0,0) = p1 - p0;
		Dm.block<2,1>(0,1) = p2 - p0;
		if(Dm.determinant()<0)
		{
			i0 = tri[0];
			i1 = tri[1];
			i2 = tri[2];
			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			Dm.block<2,1>(0,0) = p1 - p0;
			Dm.block<2,1>(0,1) = p2 - p0;

		}
		if(offset ==0)
		muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,Eigen::Vector2d::UnitX(),i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		else
		muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,Eigen::Vector2d::UnitX(),i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		mConstraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));

	}
	Eigen::VectorXd v(vertices.size()*2);
	for(int i =0;i<vertices.size();i++)
		v.block<2,1>(i*2,0) = vertices[i];

	world->AddBody(v,mConstraints,1.0);
	for(int i =0;i<vertices.size();i++)
	{
		if(vertices[i][0]<-0.49)
		{
		}
		else if(vertices[i][0]>0.49)
		{

		}
		else
			continue;
		Eigen::Vector2d fixed; 
		fixed = vertices[i];
		mAttachementConstraintVector.push_back(new AttachmentConstraint(1E6,i+offset,fixed));			
	}

	for(auto& c: mAttachementConstraintVector)
			world->AddConstraint(c);
	
	for(int i=0;i<mMuscles.size();i++)
	{
		Muscle* muscle = mMuscles[i];

		for(auto& c: muscle->muscleConstraints)
			world->AddConstraint(c);		
	}

	mActivationLevel.resize(mMuscles.size());
	mActivationLevel.setZero();

}

void
Object::
SetActivationLevel(const Eigen::VectorXd& a)
{
	mActivationLevel = a;
	for(int i=0;i<mMuscles.size();i++)
		mMuscles[i]->SetActivationLevel(a[i]);
}
double
Object::
ComputeForce(FEM::World* world)
{
	const Eigen::VectorXd& X = world->GetPositions();
	Eigen::VectorXd grad(X.rows());
	grad.setZero();
	mAttachementConstraintVector[0]->EvalGradient(X,grad);

	return grad[mAttachementConstraintVector[0]->GetI0()*2];

}
// void
// MakeMuscles(const std::string& path,Object* obj)
// {
// 	const auto& vertices = mMesh->GetVertices();
// 	const auto& triangles = mMesh->GetTriangles();
	// TiXmlDocument doc;
	// if(!doc.LoadFile(path))
 //    {
 //        std::cout<<"Cant open XML file : "<<path<<std::endl;
 //        return;
 //    }  

 //    TiXmlElement* muscles = doc.FirstChildElement("Muscles");
	// for(TiXmlElement* fiber = muscles->FirstChildElement("fiber");fiber!=nullptr;fiber = fiber->NextSiblingElement("fiber"))
	// {
	// 	std::vector<Eigen::Vector3i> muscleIndex;
	// 	for(TiXmlElement* element = fiber->FirstChildElement("element");element!=nullptr;element = element->NextSiblingElement("element"))
	// 	{
	// 		int i0,i1,i2;
	// 		i0 = std::stoi(element->Attribute("i0"));
	// 		i1 = std::stoi(element->Attribute("i1"));
	// 		i2 = std::stoi(element->Attribute("i2"));
	// 		muscleIndex.push_back(Eigen::Vector3i(i0,i1,i2));
	// 	}
   	    
 //   	    obj->AddMuscle(muscleIndex,Eigen::Vector2d::UnitX());
	// }     	    
// }

