#include "Object.h"
#include <vector>
#include <Eigen/Core>
#include <tinyxml.h>
#include "muscleDirection.h"

using namespace FEM;

#define USE_COMPLEX_MESH 1

Object::
Object(const Eigen::Vector2d& trans)
	:mMuscleStiffness(1E6),mYoungsModulus(1E6),mPoissonRatio(0.3),mMesh()
{
	Eigen::Affine2d T;
	T.setIdentity();
	T.translation() = trans;
	// mMesh = new RectanguleMesh(0.2,0.1,8,4,T);
#if USE_COMPLEX_MESH
	mMesh = new OBJLoader("../fem_test/export/simple.obj",T);
#else
	mMesh = new DiamondMesh(0.2,0.2,8,8,T);
#endif
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

	std::vector<Eigen::Vector3d> v_vector;
	for(int i=0;i<v_vector.size();i++)
	{
		Eigen::Vector3d v = v_vector[i];
	}

	for(auto v: v_vector)
	{

	}

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

		//yul
		//In hear!!!
		//Instead  Eigen::Vector2d::UnitX(), use direction from poisson's equation!
		 


       //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        // change the musccle fiber direction using Poisson's equation
        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //manually given
        //Eigen::Vector2d dir(1, 2);



#if USE_COMPLEX_MESH

        int nx = 40;
        int ny = 30;
        int nn = 15;
        Eigen::MatrixXd v;
        
        get_muscle_direction(v, nx, ny, nn);

#else

        // int nx = 12;
        // int ny = 10;
        // int nn = 5;
        int nx = 22;
        int ny = 8;
        int nn = 4;

        Eigen::MatrixXd v;
        
        get_muscle_direction(v, nx, ny, nn);
#endif
		//std::cout<<"POISSON!!! Potential : \n "<<v<<std::endl;

		Eigen::MatrixXd delta_mat_x;
		delta_mat_x.resize(ny, nx);
		delta_mat_x.setZero();

		Eigen::MatrixXd delta_mat_y;
		delta_mat_y.resize(ny, nx);
		delta_mat_y.setZero();

		Eigen::Vector2d temp_delta;
		Eigen::Matrix2d rot_90;
		rot_90(0,0) = 0;	rot_90(0,1) = -1;	rot_90(1,0) = 1;	rot_90(1,1) = 0; 
		//rot_90(0,0) = 0;	rot_90(0,1) = 1;	rot_90(1,0) = -1;	rot_90(1,1) = 0; 

		// std::cout<<"v col 15 : \n "<<v.col(15).transpose()<<std::endl;
		//std::cout<<"v(15,1) : \n "<<v(15,1)<<std::endl;
		// std::cout<<"v(0,15) : \n "<<v(0,15)<<std::endl;
		// std::cout<<"v(2,15) : \n "<<v(0,15)<<std::endl;
		// std::cout<<"v(7,2) : \n "<<v(7,2)<<std::endl;
		// std::cout<<"v(8,3) : \n "<<v(8,3)<<std::endl;
		//std::cout<<"rot_90 : \n "<<rot_90<<std::endl;

		for(int yidx = 1; yidx < ny-1 ; yidx++)
			for (int xidx = 1 ; xidx < nx-1 ; xidx++)
			{
				temp_delta(0)= v(yidx, xidx+1) - v(yidx, xidx-1);
				temp_delta(1)= v(yidx+1, xidx)-v(yidx-1, xidx);

				temp_delta = rot_90 * temp_delta;

				delta_mat_x(yidx, xidx) = temp_delta[0];
				delta_mat_y(yidx, xidx) = temp_delta[1];
			}


			//&&&&&&&&&&&$$$$$$$$$$$$
			//mapping!!
			//$$$$$$$$$$$$$$$$$$$$$$$

			Eigen::Vector2d cen = (1.0/3.0)*(p0+p1+p2);
			//Eigen::Vector2d approximate_cen (((cen[0] + 0.5) * 22), ((cen[1] + 0.1) * 8));
			
			// Eigen::Vector2d approximate_cen (abs((cen[0] + 0.5) * 22), abs((cen[1] + 0.1) * 8));
			Eigen::Vector2d approximate_cen (std::round(std::abs((cen[1] + 0.1)-0.2) * 16)+1, std::abs((cen[0] + 0.5) * 22)+1);		
			// std::cout<<"approximate_cen :  "<<approximate_cen.transpose()<<std::endl;

			//std::cout<<"delta_mat_x\n"<<delta_mat_x<<std::endl;


			Eigen::Vector2d dir (delta_mat_x(approximate_cen[0], approximate_cen[1]), delta_mat_y(approximate_cen[0], approximate_cen[1]));
			
			dir.normalize();
			//std::cout<<"dir(0) : "<<dir(0)<<std::endl;
			//if (dir(0) != dir(1))
			if (std::isnan(dir(0)))
			{
				//std::cout<<"isnan"<<std::endl;
				dir = Eigen::Vector2d::UnitX();
			}

			if (std::isinf(dir(0)))
			{
				//std::cout<<"isinf"<<std::endl;
				dir = Eigen::Vector2d::UnitX();
			}
			// std::cout<<"direction :  "<<dir.transpose()<<std::endl;

			//dir = Eigen::Vector2d::UnitX();
		if(offset !=0)
			muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,dir,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		else
			muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,dir,i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		mConstraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));


		// Eigen::MatrixXd v_indicator;
		// v_indicator.resize(ny, nx);
  //       v_indicator.setZero();



		// for(int yidx = 1; yidx < ny-1 ; yidx++)
		// 	for (int xidx = 1 ; xidx < nx-1 ; xidx++)
		// 		// int yidx = 3;
		// 		// int xidx = 2;
		// 		if ( v(yidx, xidx)< 10 && v(yidx, xidx) > -10)
		// 		{

		// 		// std::cout<<"v(yidx, xidx):"<<v(yidx, xidx)<<std::endl;
		// 		// std::cout<<"UP   v:"<<v(yidx-1, xidx)<<std::endl;
		// 		// std::cout<<"DOWN v:"<<v(yidx+1, xidx)<<std::endl;
		// 		// std::cout<<"LEFT v:"<<v(yidx, xidx-1)<<std::endl;
		// 		// std::cout<<"RIGHTv:"<<v(yidx, xidx+1)<<std::endl;
		// 		// std::cout<<"UPPER LEFT CORNER  v:"<<v(yidx-1, xidx-1)<<std::endl;
		// 		// std::cout<<"UPPER RIRHT CORNER v:"<<v(yidx-1, xidx+1)<<std::endl;
		// 		// std::cout<<"LOWER LEFT CORNER  v:"<<v(yidx+1, xidx-1)<<std::endl;
		// 		// std::cout<<"LOWER RIRHT CORNER v:"<<v(yidx+1, xidx+1)<<std::endl;


		// 			// std::cout<<"abs v:"<<std::abs(v(yidx-1,xidx-1)- v(yidx,xidx))<<std::endl;
  //        			// std::cout<<"value:"<<(v(yidx-1,xidx-1)- v(yidx,xidx))<<std::endl;

		// 			int local_min = 100;
		// 			// if (local_min >= std::abs(v(yidx-1,xidx-1)- v(yidx,xidx)))	
		// 			// {
		// 			// 	local_min = std::abs(v(yidx-1,xidx-1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 1;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx-1,xidx)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx-1,xidx)- v(yidx,xidx));
		// 			// 	// v_indicator(yidx,xidx) = 2;
		// 			// 	v_indicator(yidx,xidx) = 0;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx-1,xidx+1)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx-1,xidx+1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 3;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx,xidx+1)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx,xidx+1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 4;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx+1,xidx+1)- v(yidx,xidx)))	
		// 			// {
		// 			// 	local_min = std::abs(v(yidx+1,xidx+1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 5;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx+1,xidx)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx+1,xidx)- v(yidx,xidx));
		// 			// 	//v_indicator(yidx,xidx) = 6;
		// 			// 	v_indicator(yidx,xidx) = 0;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx+1,xidx-1)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx+1,xidx-1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 7;
		// 			// }
		// 			// if (local_min >= std::abs(v(yidx,xidx-1)- v(yidx,xidx)))
		// 			// {
		// 			// 	local_min = std::abs(v(yidx,xidx-1)- v(yidx,xidx));
		// 			// 	v_indicator(yidx,xidx) = 8;
		// 			// }


		// 			//		//
		// 			// <--> //
		// 			//		//
		// 			//plain
		// 			if (local_min >= std::abs(v(yidx,xidx-1)- v(yidx,xidx)) || local_min >= std::abs(v(yidx,xidx+1)- v(yidx,xidx)))
		// 			{
		// 				if (std::abs(v(yidx,xidx-1)- v(yidx,xidx)) < std::abs(v(yidx,xidx+1)- v(yidx,xidx)))
		// 					local_min = std::abs(v(yidx,xidx-1)- v(yidx,xidx));
		// 				else
		// 					local_min = std::abs(v(yidx,xidx+1)- v(yidx,xidx));

		// 				v_indicator(yidx,xidx) = 1;
		// 			}
					
		// 			// \	//
		// 			//  \   //
		// 			//	 \	//
		// 			// upper - left TO lower-Right : 1,5
		// 			if (local_min >= std::abs(v(yidx-1,xidx-1)- v(yidx,xidx)) || local_min >= std::abs(v(yidx+1,xidx+1)- v(yidx,xidx)))	
		// 			{
		// 				if (std::abs(v(yidx-1,xidx-1)- v(yidx,xidx)) > std::abs(v(yidx+1,xidx+1)- v(yidx,xidx)))
		// 					local_min = std::abs(v(yidx+1,xidx+1)- v(yidx,xidx));
		// 				else
		// 					local_min = std::abs(v(yidx-1,xidx-1)- v(yidx,xidx));
		// 				v_indicator(yidx,xidx) = 2;
		// 			}

		// 			// 	 /  //
		// 			//  /   //
		// 			// /	//
		// 			// upper-right TO lower-left : 3,7
		// 			if (local_min >= std::abs(v(yidx-1,xidx+1)- v(yidx,xidx)) || local_min >= std::abs(v(yidx+1,xidx-1)- v(yidx,xidx)))
		// 			{
		// 				if (std::abs(v(yidx-1,xidx+1)- v(yidx,xidx)) > std::abs(v(yidx+1,xidx-1)- v(yidx,xidx)))
		// 					local_min = std::abs(v(yidx+1,xidx-1)- v(yidx,xidx));
		// 				else
		// 					local_min = std::abs(v(yidx-1,xidx+1)- v(yidx,xidx));
		// 				v_indicator(yidx,xidx) = 3;
		// 			}

		// 		}


		// //std::cout<<"v_indicator : \n "<<v_indicator<<std::endl;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		// if(offset ==0)
		// muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,Eigen::Vector2d::UnitX(),i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		// else
		// muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,Eigen::Vector2d::UnitX(),i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		// mConstraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0+offset,i1+offset,i2+offset,1.0/6.0*(Dm.determinant()),Dm.inverse()));

	}
	Eigen::VectorXd v(vertices.size()*2);
	for(int i =0;i<vertices.size();i++)
		v.block<2,1>(i*2,0) = vertices[i];

	world->AddBody(v,mConstraints,1.0);

#if USE_COMPLEX_MESH
	Eigen::Vector2d fixed; 
	fixed = vertices[5];
	mAttachementConstraintVector.push_back(new AttachmentConstraint(1E6,5+offset,fixed));
	fixed = vertices[395];
	mAttachementConstraintVector.push_back(new AttachmentConstraint(1E6,395+offset,fixed));
#else

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
#endif
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

