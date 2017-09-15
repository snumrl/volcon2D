#include "GL_function.h"
#include "GL/glut.h"

void
DrawStringOnScreen(float _x, float _y, const std::string& _s,bool _bigFont,const Eigen::Vector3d& color)
{
    // draws text on the screen
    GLint oldMode;
    glGetIntegerv(GL_MATRIX_MODE, &oldMode);
    glMatrixMode(GL_PROJECTION);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRasterPos2f(_x, _y);
    glColor3f(color[0],color[1],color[2]);
    unsigned int length = _s.length();
    for (unsigned int c = 0; c < length; c++) {
    if (_bigFont)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, _s.at(c) );
    else
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, _s.at(c) );
    }  
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(oldMode);
}
void
DrawPoint(const Eigen::Vector2d& p,const Eigen::Vector3d& color)
{
    glBegin(GL_POINTS);
    glColor3f(color[0],color[1],color[2]);
    glVertex2f(p[0],p[1]);
    glEnd();
}

void
DrawLines(const Eigen::Vector2d& p,const Eigen::Vector2d& q,const Eigen::Vector3d& color)
{
    glBegin(GL_LINES);
    glColor3f(color[0],color[1],color[2]);
    glVertex2f(p[0],p[1]);
    glVertex2f(q[0],q[1]);
    glEnd();
}

void
DrawTriangle(const Eigen::Vector2d& a,const Eigen::Vector2d& b,const Eigen::Vector2d& c,const Eigen::Vector3d& color)
{
    glBegin(GL_TRIANGLES);
    glColor3f(color[0],color[1],color[2]);
    glVertex2f(a[0],a[1]);
    glVertex2f(b[0],b[1]);
    glVertex2f(c[0],c[1]);
    glEnd();
}
void
DrawArrow(const Eigen::Vector2d& p, const Eigen::Vector2d& v,const Eigen::Vector3d& color)
{
    double thickness = 0.004;
    glLineWidth(thickness);
    glBegin(GL_LINES);
    glColor3f(color[0],color[1],color[2]);
    glVertex2f(p[0], p[1]);
    glVertex2f(p[0]+v[0], p[1]+v[1]);
    glEnd();

    // draw arrowhead as a triangle
    double theta = atan2(v[1], v[0]);
    
    glPushMatrix();
    glTranslatef(p[0]+v[0], p[1]+v[1], 0.0);
    glRotatef(theta*180.0/3.141592, 0.0, 0.0, 1.0);
    // glTranslatef(thickness, 0.0, 0.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(0.0, thickness);
    glVertex2f(2*thickness, 0.0);
    glVertex2f(0.0, -thickness);
    glEnd();
    glPopMatrix();
    glLineWidth(1.0);
}
void
DrawSphere(const Eigen::Isometry3d& T,const float& r,const Eigen::Vector3d& color)
{
    glPushMatrix();
    glMultMatrixd(T.matrix().data());
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(color[0],color[1],color[2]);
    glVertex2f(0,0);
    for(float phi = 0;phi<=3.141592*2.0+1E-4;phi+=3.141592/8.0)
    {
        float x = r*cos(phi);
        float y = r*sin(phi);                            
        glVertex2f(x,y);
    }
    glEnd();
    glColor3f(0,0,0);
    glBegin(GL_LINE_LOOP);
    for(double phi = 0;phi<=3.141592*2;phi+=3.141592/8.0)
        glVertex2f(r*cos(phi),r*sin(phi));
    glEnd();

    glPopMatrix();
}
void
DrawCapsule(const Eigen::Isometry3d& T,float w,float h,const Eigen::Vector3d& color)
{
    glPushMatrix();
    glMultMatrixd(T.matrix().data());
    if(w<h)
    {
        double diff = h - w;
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_QUADS);
        glVertex2f(w/2,diff/2);
        glVertex2f(-w/2,diff/2);
        glVertex2f(-w/2,-diff/2);
        glVertex2f(w/2,-diff/2);
        glEnd();

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0,diff/2);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = 0.5*w*cos(phi);
            double y = diff/2 + 0.5*w*sin(phi);                            
            glVertex2f(x,y);
        }
        glEnd();

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0,-diff/2);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = 0.5*w*cos(phi);
            double y = -diff/2 - 0.5*w*sin(phi);                            
            glVertex2f(x,y);
        }
        glEnd();
        glColor3f(0,0,0);
        glBegin(GL_LINES);
        glVertex2f(w/2,diff/2);
        glVertex2f(-w/2,diff/2);

        glVertex2f(-w/2,diff/2);
        glVertex2f(-w/2,-diff/2);

        glVertex2f(-w/2,-diff/2);
        glVertex2f(w/2,-diff/2);

        glVertex2f(w/2,-diff/2);
        glVertex2f(w/2,diff/2);
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex2f(0,diff/2);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = 0.5*w*cos(phi);
            double y = diff/2 + 0.5*w*sin(phi);                            
            glVertex2f(x,y);
        }
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex2f(0,-diff/2);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = 0.5*w*cos(phi);
            double y = -diff/2 - 0.5*w*sin(phi);                            
            glVertex2f(x,y);
        }
        glEnd();
    }
    else
    {
        double diff = w - h;
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_QUADS);
        glVertex2f(diff/2,h/2);
        glVertex2f(-diff/2,h/2);
        glVertex2f(-diff/2,-h/2);
        glVertex2f(diff/2,-h/2);
        glEnd();

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(diff/2,0);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = diff/2 +0.5*h*sin(phi);
            double y = 0.5*h*cos(phi);                            
            glVertex2f(x,y);
        }
        glEnd();

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(-diff/2,0);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = -diff/2 - 0.5*h*sin(phi);
            double y = 0.5*h*cos(phi);                            
            glVertex2f(x,y);
        }
        glEnd();
        glColor3f(0,0,0);
        glBegin(GL_LINES);
        glVertex2f(diff/2,h/2);
        glVertex2f(-diff/2,h/2);
        glVertex2f(-diff/2,-h/2);
        glVertex2f(diff/2,-h/2);
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex2f(diff/2,0);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = diff/2 + 0.5*h*sin(phi);
            double y = 0.5*h*cos(phi);                            
            glVertex2f(x,y);
        }
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex2f(-diff/2,0);
        for(double phi = 0;phi<=3.141592;phi+=3.141592/8.0)
        {
            double x = -diff/2 - 0.5*h*sin(phi);
            double y = 0.5*h*cos(phi);                            
            glVertex2f(x,y);
        }
        glEnd();
    }
    glPopMatrix();
}
