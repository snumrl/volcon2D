#include "Mesh.h"
#include "RectangleMesh.h"
#include <cassert>
#include <algorithm>
#include <iostream>
using namespace FEM;
RectangleMesh::
RectangleMesh(double _w,double _h,int _nw,int _nh,const Eigen::Affine2d& T)
	:Mesh(),mX(_w),mY(_h),mNx(_nw),mNy(_nh)
{
	double dx = mX /(double)mNx;
	double dy = mY /(double)mNy;

	for(int i=0;i<mNx;i++)
	{
		for(int j=0;j<mNy;j++)
		{
			mVertices.push_back(Eigen::Vector2d(i*dx-0.5*mX,j*dy-0.5*mY));
			mVertices.push_back(Eigen::Vector2d(i*dx-0.5*mX+0.5*dx,j*dy-0.5*mY+0.5*dy));
		}
		mVertices.push_back(Eigen::Vector2d(i*dx-0.5*mX,mNy*dy-0.5*mY));
	}
	for(int j=0;j<mNy+1;j++)
	{
		mVertices.push_back(Eigen::Vector2d(mNx*dx-0.5*mX,j*dy-0.5*mY));
	}


	for(auto& v : mVertices)
		v = T*v;

	for(int i=0 ;i<mNx-1;i++)
		for(int j=0;j<mNy;j++)
		{
			int index = (2*mNy+1)*i+2*j;
			int n21 = 2*mNy+1;
			mTriangles.push_back(Eigen::Vector3i(index+1,index,index+2));
			mTriangles.push_back(Eigen::Vector3i(index+1,index+n21,index));
			mTriangles.push_back(Eigen::Vector3i(index+1,index+2,index+2+n21));
			mTriangles.push_back(Eigen::Vector3i(index+1,index+2+n21,index+n21));
		}

	for(int j=0;j<mNy;j++)
	{
		int index = (2*mNy+1)*(mNx-1)+2*j;
		int n21 = 2*mNy+1-j;
		mTriangles.push_back(Eigen::Vector3i(index+1,index,index+2));
		mTriangles.push_back(Eigen::Vector3i(index+1,index+n21,index));
		mTriangles.push_back(Eigen::Vector3i(index+1,index+2,index+1+n21));
		mTriangles.push_back(Eigen::Vector3i(index+1,index+1+n21,index+n21));
	}
}

