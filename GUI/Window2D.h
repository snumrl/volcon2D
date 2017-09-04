#ifndef __WINDOW2D_H__
#define __WINDOW2D_H__
#include <vector>
class Camera2D;

class Window2D
{
protected:
	Camera2D* 	mCamera;
	bool 		mIsDrag;
	int 			mMouseType;
	int 			mPrevX,mPrevY;
	int 			mDisplayTimeout;
public:
	Window2D();
	~Window2D();

	virtual void InitWindow(int _w,int _h,const char* _name);
	static void DisplayEvent();
	static void KeyboardEvent(unsigned char key,int x,int y);
	static void MouseEvent(int button, int state, int x, int y);
	static void MotionEvent(int x, int y);
	static void ReshapeEvent(int w, int h);
	static void TimerEvent(int value);

	static Window2D* current();
	static std::vector<Window2D*> mWindows;
	static std::vector<int> mWinIDs;
protected:

	virtual void Display() = 0;
	virtual void Keyboard(unsigned char key,int x,int y) = 0;
	virtual void Mouse(int button, int state, int x, int y) = 0;
	virtual void Motion(int x, int y) = 0;
	virtual void Reshape(int w, int h) = 0;
	virtual void Timer(int value) = 0;
};

#endif
