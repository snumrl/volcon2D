#ifndef __DART_INTERFACE_H__
#define __DART_INTERFACE_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

void MakeRootBody(const dart::dynamics::SkeletonPtr& skel, const std::string& name,const Eigen::Vector3d& size,const Eigen::Vector3d& c_to_joint,const double& mass);
void MakeBody(const dart::dynamics::SkeletonPtr& skel,const dart::dynamics::BodyNodePtr& parent,const std::string& name,const Eigen::Vector3d& size,const Eigen::Vector3d& p_to_joint,const Eigen::Vector3d& c_to_joint,const double& mass);
void MakeWeldBody(const dart::dynamics::SkeletonPtr& skel,const dart::dynamics::BodyNodePtr& parent,const std::string& name,const double& radius,const Eigen::Vector3d& p_to_joint,const Eigen::Vector3d& c_to_joint,const double& mass);
void MakeBall(const dart::dynamics::SkeletonPtr& skel,const double& radius,const double& mass);
void DrawSkeleton(const dart::dynamics::SkeletonPtr& skel);
#endif
