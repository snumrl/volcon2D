#ifndef __GL_FUNCTION_H__
#define __GL_FUNCTION_H__
#include <Eigen/Core>
#include <Eigen/Geometry>

void DrawStringOnScreen(float _x, float _y, const std::string& _s,bool _bigFont);
void DrawPoint(const Eigen::Vector2d& p);
void DrawLines(const Eigen::Vector2d& p,const Eigen::Vector2d& q);
void DrawTriangle(const Eigen::Vector2d& a,const Eigen::Vector2d& b,const Eigen::Vector2d& c);
void DrawArrow(const Eigen::Vector2d& p, const Eigen::Vector2d& v);
void DrawSphere(const Eigen::Isometry3d& T,const float& r);
void DrawCapsule(const Eigen::Isometry3d& T,float w,float h);

#endif
