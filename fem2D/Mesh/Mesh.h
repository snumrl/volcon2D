#ifndef __MESH_H__
#define __MESH_H__
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace FEM
{
class Mesh
{
protected:
	std::vector<Eigen::Vector2d> mVertices;
	std::vector<Eigen::Vector3i> mTriangles;
public:
	Mesh(){};
	virtual const std::vector<Eigen::Vector3i>& GetTriangles(){return mTriangles;};
	virtual const std::vector<Eigen::Vector2d>& GetVertices(){return mVertices;};
	virtual void Clear() {mVertices.clear(); mTriangles.clear();};
};
};
#endif
