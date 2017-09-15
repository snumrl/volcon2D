#ifndef __GL_FUNCTION_H__
#define __GL_FUNCTION_H__
#include <Eigen/Core>
#include <Eigen/Geometry>

void DrawStringOnScreen(float _x, float _y, const std::string& _s,bool _bigFont,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawPoint(const Eigen::Vector2d& p,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawLines(const Eigen::Vector2d& p,const Eigen::Vector2d& q,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawTriangle(const Eigen::Vector2d& a,const Eigen::Vector2d& b,const Eigen::Vector2d& c,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawArrow(const Eigen::Vector2d& p, const Eigen::Vector2d& v,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawSphere(const Eigen::Isometry3d& T,const float& r,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));
void DrawCapsule(const Eigen::Isometry3d& T,float w,float h,const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8));

#endif
