#include "Mesh.h"
#include "RectangleMesh.h"
#include "DiamondMesh.h"
#include <cassert>
#include <algorithm>
#include <iostream>
using namespace FEM;
bool is_big(Eigen::Vector2d& p)
{
	if(p.norm()>1E6)
		return true;
	else
		return false;
}
DiamondMesh::
DiamondMesh(double _w,double _h,int _nw,int _nh,const Eigen::Affine2d& T)
	:RectangleMesh(_w,_h,_nw,_nh,T)
{
	assert(mNx>=mNy);

	for(auto& v : mVertices)
		v = T.inverse()*v;
	std::vector<int> delete_index;

	for(int i =0;i<mVertices.size();i++)
		if(!CheckInside(mVertices[i]))
			delete_index.push_back(i);

	std::vector<int> accumulate_index;
	accumulate_index.resize(mVertices.size(),0);

	for(int i =0;i<delete_index.size();i++)
	{
		for(int j=delete_index[i];j<mVertices.size();j++)
			accumulate_index[j]+=1;
	}

	for(int i=0;i<delete_index.size();i++){
		mVertices[delete_index[i]] = Eigen::Vector2d(1E8,1E8);
	}
	mVertices.erase(std::remove_if (mVertices.begin(), mVertices.end(), is_big), mVertices.end());
	for(int i=0;i<mTriangles.size();i++)
	{
		Eigen::Vector3i tri = mTriangles[i];
		bool is_outside = false;
		for(int j=0;j<3;j++)
			if(find(delete_index.begin(), delete_index.end(), tri[j]) != delete_index.end())
				is_outside = true;
		if(is_outside){
			mTriangles.erase(mTriangles.begin()+i);
			i--;
		}
		else
			mTriangles[i] = Eigen::Vector3i(
									tri[0]-accumulate_index[tri[0]],
									tri[1]-accumulate_index[tri[1]],
									tri[2]-accumulate_index[tri[2]]);

	}

	double min_x = 2;
	mStartPointIndex = -1;
	double max_x = -2;
	mEndPointIndex = -1;
	for(int i =0;i<mVertices.size();i++)
	{
		if(mVertices[i][0]<min_x)
		{
			mStartPointIndex = i;
			min_x = mVertices[i][0];
		}
		if(mVertices[i][0]>max_x)
		{
			mEndPointIndex = i;
			max_x = mVertices[i][0];
		}
	}
	double length = max_x-min_x;
	for(auto& v : mVertices){
		v[0] *= 1.0/length;
	}
	for(auto& v : mVertices)
		v = T*v;
}

bool 
DiamondMesh::
CheckInside(const Eigen::Vector2d& p)
{
	bool ret = true;
	if(fabs(p[0] - p[1])>0.5*mX+1E-5)
		ret = false;
	if(fabs(p[0] + p[1])>0.5*mX+1E-5)
		ret = false;
	
	return ret;
}
