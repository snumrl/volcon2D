#ifndef __DIAMOND_MESH_H__
#define __DIAMOND_MESH_H__
#include <Eigen/Core>

namespace FEM
{
class RectangleMesh;
class DiamondMesh : public RectangleMesh
{
protected:
	int mStartPointIndex,mEndPointIndex;
public:
	DiamondMesh(double _w=1.0,double _h=1.0,int _nw = 5,int _nh = 5,const Eigen::Affine2d& T = Eigen::Affine2d::Identity());
	bool CheckInside(const Eigen::Vector2d& p);
	int GetStartingPointIndex(){return mStartPointIndex;};
	int GetEndingPointIndex(){return mEndPointIndex;};
};
};
#endif
