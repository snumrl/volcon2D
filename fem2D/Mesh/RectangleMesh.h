#ifndef __RECTANGLE_MESH_H__
#define __RECTANGLE_MESH_H__
#include <Eigen/Core>

namespace FEM
{
class Mesh;
class RectangleMesh : public Mesh
{
protected:
	int mNx,mNy;
	double mX,mY;
public:
	RectangleMesh(double _w=1.0,double _h=1.0,int _nw = 5,int _nh = 5,const Eigen::Affine2d& T = Eigen::Affine2d::Identity());
	
};
};
#endif
