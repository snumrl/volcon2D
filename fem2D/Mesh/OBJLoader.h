#ifndef __OBJLOADER_H__
#define __OBJLOADER_H__
#include <Eigen/Core>

namespace FEM
{
class Mesh;

class OBJLoader : public Mesh
{
public:
	OBJLoader();
	OBJLoader(const std::string& obj_file,const Eigen::Affine2d& T = Eigen::Affine2d::Identity());

	void Load(const std::string& obj_file,const Eigen::Affine2d& T = Eigen::Affine2d::Identity());
};

};


#endif
