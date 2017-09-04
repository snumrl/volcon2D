#ifndef __CAMERA_2D__
#define __CAMERA_2D__
#include <Eigen/Core>


class Camera2D
{
private:
	double mCx,mCy;
	double mdx;

public:
	Camera2D();

	void Apply();
	void Pan(int x,int y,int prev_x,int prev_y);
	void Translate(int x,int y,int prev_x,int prev_y);
	Eigen::Vector2d GetWorldPosition(int x,int y);
};


#endif