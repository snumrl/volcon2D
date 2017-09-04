#include "Window2D.h"
#include "Camera2D.h"
#include <iostream>
#include <Eigen/Geometry>
#include "GL/glut.h"
std::vector<Window2D*> Window2D::mWindows;
std::vector<int> Window2D::mWinIDs;
Window2D::
Window2D()
	:mCamera(new Camera2D()),mIsDrag(false),mMouseType(0),mPrevX(0),mPrevY(0),mDisplayTimeout(1.0/30.0)
{
}
Window2D::
~Window2D()
{
}
void
Window2D::
InitWindow(int _w,int _h,const char* _name)
{
	mWindows.push_back(this);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE | GLUT_ACCUM);
	glutInitWindowPosition(150, 100);
	glutInitWindowSize(_w, _h);
	mWinIDs.push_back(glutCreateWindow(_name));
	  glutDisplayFunc(DisplayEvent);
  glutReshapeFunc(ReshapeEvent);
  glutKeyboardFunc(KeyboardEvent);
  glutMouseFunc(MouseEvent);
  glutMotionFunc(MotionEvent);
  glutTimerFunc(mDisplayTimeout, TimerEvent, 0);
}
inline Window2D* Window2D::current() {
  int id = glutGetWindow();
  for (unsigned int i = 0; i < mWinIDs.size(); i++) {
    if (mWinIDs.at(i) == id) {
      return mWindows.at(i);
    }
  }
  std::cout << "An unknown error occurred!" << std::endl;
  exit(0);
}
void
Window2D::
DisplayEvent()
{
	current()->Display();
}
void
Window2D::
KeyboardEvent(unsigned char key,int x,int y)
{
	current()->Keyboard(key,x,y);
}
void
Window2D::
MouseEvent(int button, int state, int x, int y)
{
	current()->Mouse(button,state,x,y);
}
void
Window2D::
MotionEvent(int x, int y)
{
	current()->Motion(x,y);
}
void
Window2D::
ReshapeEvent(int w, int h)
{
	current()->Reshape(w,h);
}
void
Window2D::
TimerEvent(int value)
{
	current()->Timer(value);
}