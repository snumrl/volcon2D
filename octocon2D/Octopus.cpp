#include "Octopus.h"
#include <vector>
#include <Eigen/Core>
#include <tinyxml.h>

using namespace FEM;

Octopus::
Octopus()
	:mMuscleStiffness(2E6),mYoungsModulus(1E7),mPoissonRatio(0.3),mMesh(),mTarget()
{
	// mMesh = new RectangleMesh(0.2,1,2,20);
	mMesh = new OBJLoader("../octocon2D/export/octo_ver4.obj");
}

void 
Octopus::
AddTarget(const Target& target) {
	mTarget.push_back(target);
	// std::cout << mTarget.size() << std::endl;
	// std::cout << target.idx << std::endl;
}

std::vector<Target>
Octopus::
GetTarget() {
	return mTarget;
}

void 
Octopus::
SolveSoftIK(FEM::World* world) {
	std::cout << "SolveSoftIK" << std::endl;
	Eigen::VectorXd X = world->GetPositions();

	Eigen::VectorXd X_prime(world->GetNumVertices()*2);
	X_prime.setZero();
	for(int i=0; i<mTarget.size(); i++) {
		X_prime.block<2,1>(mTarget[i].idx*2,0) = X.block<2,1>(mTarget[i].idx*2,0) - mTarget[i].coord;
	}

	
}

void
Octopus::
AddMuscle(
	const std::vector<Eigen::Vector3i> indexList,
	const Eigen::Vector2d& fiber_direction
	)
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

			// muscle->constraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));
			// muscle->constraints.push_back(new SpringConstraint(10000.0,i0,i1,(p0-p1).norm()));
			// muscle->constraints.push_back(new SpringConstraint(10000.0,i1,i2,(p1-p2).norm()));
			// muscle->constraints.push_back(new SpringConstraint(10000.0,i2,i0,(p2-p0).norm()));
			muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
			// muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		}
	}

	// Add Attachment 
	// AttachmentConstraint* ground = new AttachmentConstraint(1E8,0,Eigen::Vector2d(0,-0.5));
	// mAttachementConstraintVector.push_back(ground);
	std::vector<AttachmentConstraint*> ground;
	std::vector<int> groundIdx;

	for(int i=981; i<=995; i++) {
		groundIdx.push_back(i);
	}

	for(int i=0; i<groundIdx.size(); i++) {		
		Eigen::Vector2d fixed; 
		fixed[0] = i / 20.0 - 0.35;
		fixed[1] = 1.15;
		mAttachementConstraintVector.push_back(new AttachmentConstraint(1E8,groundIdx[i],fixed));
	}


	muscle->activationLevel = 0.0;
}
void
Octopus::
Initialize(FEM::World* world)
{
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

			mConstraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));
		}
	}

	Eigen::VectorXd v(vertices.size()*2);
	for(int i =0;i<vertices.size();i++)
		v.block<2,1>(i*2,0) = vertices[i];

	world->AddBody(v,mConstraints,1.0);
	
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
Octopus::
SetActivationLevel(const Eigen::VectorXd& a)
{
	// std::cout << "****" << mMuscles.size() << std::endl;
	mActivationLevel = a;
	for(int i=0;i<mMuscles.size();i++)
		for(auto& mc : mMuscles[i]->muscleConstraints)
			mc->SetActivationLevel(a[i]);
}

void
MakeMuscles(const std::string& path,Octopus* ms)
{
	TiXmlDocument doc;
	if(!doc.LoadFile(path))
    {
        std::cout<<"Cant open XML file : "<<path<<std::endl;
        return;
    }  

    TiXmlElement* muscles = doc.FirstChildElement("Muscles");

    for(TiXmlElement* leg = muscles->FirstChildElement("leg");leg!=nullptr;leg = leg->NextSiblingElement("leg")) {
    	for(TiXmlElement* fiber = leg->FirstChildElement("fiber");fiber!=nullptr;fiber = fiber->NextSiblingElement("fiber")) {
    	
	    	int start_num, end_num;
	    	int start_idx1, end_idx1;
	    	int start_idx2, end_idx2;
	    	int start_idx3, end_idx3;

	    	for(TiXmlElement* start = fiber->FirstChildElement("start");start!=nullptr;start = start->NextSiblingElement("start")) {
	        	start_num = std::stod(start->Attribute("num"));
	        	start_idx1 = std::stod(start->Attribute("idx1"));
				start_idx2 = std::stod(start->Attribute("idx2"));
				start_idx3 = std::stod(start->Attribute("idx3"));			
	   	    }

	   	    for(TiXmlElement* end = fiber->FirstChildElement("end");end!=nullptr;end = end->NextSiblingElement("end")) {
	        	end_num = std::stod(end->Attribute("num"));
	        	end_idx1 = std::stod(end->Attribute("idx1"));
				end_idx2 = std::stod(end->Attribute("idx2"));
				end_idx3 = std::stod(end->Attribute("idx3"));
	   	    }   

	   	    int cellNum = end_num - start_num + 1;
	   	    int dist = start_idx3 - start_idx1;

	   	    std::vector<Eigen::Vector3i> muscleIndex;

	   	    for(int i=0; i<cellNum/2; i++) {
	   	    	int idx1 = start_idx1 + dist*i;
	   	    	int idx2 = start_idx2 + dist*i;
	   	    	int idx3 = start_idx3 + dist*i;

	   	    	muscleIndex.push_back(Eigen::Vector3i(idx1,idx2,idx3));
	   	    }

			for(int i=0; i<cellNum/2; i++) {
				int idx1 = end_idx1 - dist*i;
	   	    	int idx2 = end_idx2 - dist*i;
	   	    	int idx3 = end_idx3 - dist*i;

	   	    	muscleIndex.push_back(Eigen::Vector3i(idx1,idx2,idx3));    
	   	    }   	

	   	    ms->AddMuscle(muscleIndex,Eigen::Vector2d::UnitY());
   	    }     	    
    }

    // std::cout << "# (Muscle Fiber) : " << muscleIndex.size() << std::endl;
}

