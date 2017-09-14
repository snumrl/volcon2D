#include "Octopus.h"


void
Octopus::
AddMuscle(const std::vector<int>& element_index_list,const Eigen::Vector3d& fiber_direction)
{
	mMuscles.push_back(new Muscle());
	auto& muscle = mMuscles.back();
	const auto& vertices = mMesh->GetVertices();
	const auto& triangles = mMesh->GetTriangles();
	for(auto index : element_index_list)
	{
		int i0,i1,i2;
		Eigen::Vector2d p0,p1,p2;

		i0 = triangles[index][0];
		i1 = triangles[index][1];
		i2 = triangles[index][2];

		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		Eigen::Matrix2d Dm;

		Dm.block<2,1>(0,0) = p1 - p0;
		Dm.block<2,1>(0,1) = p2 - p0;
		if(Dm.determinant()<0)
		{
			i0 = tri[0];
			i1 = tri[2];
			i2 = tri[1];
			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			Dm.block<2,1>(0,0) = p1 - p0;
			Dm.block<2,1>(0,1) = p2 - p0;
		}
		muscle->constraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
	}
}

void
Octopus::
Initialize(FEM::World* world)
{
	const auto& vertices = mMesh->GetVertices();
	const auto& triangles = mMesh->GetTriangles();
	std::vector<FEM::Constraint*> constraints;
	Eigen::VectorXd v(vertices.size()*2);
	for(int i =0;i<vertices.size();i++)
		v.block<2,1>(i*2,0) = vertices[i];
	for(const auto& tri : triangles)
	{
		int i0,i1,i2;
		Eigen::Vector2d p0,p1,p2;
		
		i0 = tri[0];
		i1 = tri[1];
		i2 = tri[2];
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		Eigen::Matrix2d Dm;

		Dm.block<2,1>(0,0) = p1 - p0;
		Dm.block<2,1>(0,1) = p2 - p0;
		if(Dm.determinant()<0)
		{
			i0 = tri[0];
			i1 = tri[2];
			i2 = tri[1];
			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			Dm.block<2,1>(0,0) = p1 - p0;
			Dm.block<2,1>(0,1) = p2 - p0;
		}

		constraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));
	}

}