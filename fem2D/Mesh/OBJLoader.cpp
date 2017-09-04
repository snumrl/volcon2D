#include "Mesh.h"
#include "OBJLoader.h"
#include <fstream>
#include <iostream>

using namespace FEM;
void
MakeSlashToSpace(std::string& str)
{
	for(int i=0;i<str.size();i++)
		if(str[i]=='/')
			str[i] = ' ';
}


OBJLoader::
OBJLoader()
	:Mesh()
{

}
OBJLoader::
OBJLoader(const std::string& obj_file,const Eigen::Affine2d& T)
	:Mesh()
{
	Load(obj_file,T);
}

void
OBJLoader::
Load(const std::string& obj_file,const Eigen::Affine2d& T)
{
	std::ifstream ifs(obj_file);
	if(!(ifs.is_open()))
	{
		std::cout<<"Can't read file "<<obj_file<<std::endl;
		return;
	}
	std::string str;
	std::string index;
	std::stringstream ss;

	while(!ifs.eof())
	{
		str.clear();
		index.clear();
		ss.clear();

		std::getline(ifs,str);
		ss.str(str);
		ss>>index;

		if(!index.compare("v"))
		{
			double x,y;
			ss>>x>>y;
			mVertices.push_back(Eigen::Vector2d(x,y));
		}
		else if(!index.compare("f"))
		{
			std::string token;
			int index[3];
			for(int i=0;i<3;i++)
			{
				int vi = 0,vti = 0,vni = 0;
				ss>>token;
				MakeSlashToSpace(token);
				std::stringstream ss_token(token);
				ss_token>>vi;
				index[i] = vi-1;
			}
			mTriangles.push_back(Eigen::Vector3i(index[0],index[1],index[2]));
		}
	}
	ifs.close();

	for(auto& v : mVertices)
		v = T*v;
}
