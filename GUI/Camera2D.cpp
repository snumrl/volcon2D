#include "Camera2D.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <GL/glut.h>
Camera2D::
Camera2D()
	:mCx(0.0),mCy(0.5),mdx(2.0)
{

}

void
Camera2D::
Apply()
{
	GLint w = glutGet(GLUT_WINDOW_WIDTH);
	GLint h = glutGet(GLUT_WINDOW_HEIGHT);
	double aspect_ratio_inv = (double)h/(double)w;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(mCx - mdx*0.5,mCx + mdx*0.5,
			mCy - mdx*aspect_ratio_inv*0.5,mCy + mdx*aspect_ratio_inv*0.5,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
void
Camera2D::
Pan(int x,int y,int prev_x,int prev_y)
{
	float delta = (float)prev_y - (float)y;
	delta*=0.01;

	mdx += delta;
	if(mdx<0.1)
		mdx=0.1;
	else if(mdx>100)
		mdx = 100.0;
	
}
void
Camera2D::
Translate(int x,int y,int prev_x,int prev_y)
{
	Eigen::Vector2d delta = GetWorldPosition(prev_x,prev_y)-GetWorldPosition(x,y);

	mCx += delta[0];
	mCy += delta[1];
}
Eigen::Vector2d
Camera2D::
GetWorldPosition(int x,int y)
{
	GLint w = glutGet(GLUT_WINDOW_WIDTH);
	GLint h = glutGet(GLUT_WINDOW_HEIGHT);

	Eigen::Vector2d screen_pos((double)x/(double)w-0.5,(1.0-(double)y/(double)h)-0.5);
	Eigen::Vector2d pos = mdx*screen_pos+Eigen::Vector2d(mCx,mCy);

	return pos;
}