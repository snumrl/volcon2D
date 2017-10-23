#include "DART_Interface.h"
#include "GUI/GL_function.h"
#include <Eigen/Geometry>
#include <vector>
#include "GL/glut.h"
using namespace dart::dynamics;
using namespace dart::simulation;
/*************************************************************
    MakeRootBody function
This function makes root body for skeleton. 
globalOffset param is center of mass of body.
*************************************************************/
void
MakeRootBody(const SkeletonPtr& skel, const std::string& name,const Eigen::Vector3d& size,const Eigen::Vector3d& c_to_joint,const double& mass)
{
    ShapePtr shape = std::shared_ptr<BoxShape>(new BoxShape(size));
    dart::dynamics::Inertia inertia;
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));

    RevoluteJoint::Properties prop;
    prop.mName = name + "_joint";
    prop.mAxis = Eigen::Vector3d::UnitZ();
    prop.mT_ChildBodyToJoint.translation() = c_to_joint;

    BodyNodePtr bn = skel->createJointAndBodyNodePair<RevoluteJoint>(
      nullptr,prop,BodyNode::AspectProperties(name)).second;
    auto sn = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(shape);
    bn->setInertia(inertia);
}
/*************************************************************
    MakeBody function
This function makes articulated body of skeleton.
Input is straightforward, except globaloffset and jointoffset.
globaloffset parameter is center of mass of body expressed in world coordinate.
jointoffset parameter is position of joint expressed in child(this) body coordinate. 
*************************************************************/
void
MakeBody(const SkeletonPtr& skel,const BodyNodePtr& parent,const std::string& name,
    const Eigen::Vector3d& size,
    const Eigen::Vector3d& p_to_joint,
    const Eigen::Vector3d& c_to_joint,const double& mass)
{
    ShapePtr shape = std::shared_ptr<BoxShape>(new BoxShape(size));
    dart::dynamics::Inertia inertia;
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));

    RevoluteJoint::Properties prop;
    prop.mName = name + "_joint";
    prop.mAxis = Eigen::Vector3d::UnitZ();
    prop.mT_ParentBodyToJoint.translation() = p_to_joint;
    prop.mT_ChildBodyToJoint.translation() = c_to_joint;

    BodyNodePtr bn = skel->createJointAndBodyNodePair<RevoluteJoint>(
      parent,prop,BodyNode::AspectProperties(name)).second;

    auto sn = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(shape);

    bn->setInertia(inertia);
}
void MakeWeldBody(const SkeletonPtr& skel,
    const BodyNodePtr& parent,const std::string& name,const double& radius,
    const Eigen::Vector3d& p_to_joint,const Eigen::Vector3d& c_to_joint,const double& mass)
{
    ShapePtr shape = std::shared_ptr<SphereShape>(new SphereShape(radius));
    dart::dynamics::Inertia inertia;
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));

    WeldJoint::Properties prop;
    prop.mName = name + "_joint";
    prop.mT_ParentBodyToJoint.translation() = p_to_joint;
    prop.mT_ChildBodyToJoint.translation() = c_to_joint;

    BodyNodePtr bn = skel->createJointAndBodyNodePair<WeldJoint>(
      parent,prop,BodyNode::AspectProperties(name)).second;

    auto sn = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(shape);

    bn->setInertia(inertia);
}
void
MakeBall(const SkeletonPtr& skel,const double& radius,const double& mass)
{
    ShapePtr shape = std::shared_ptr<SphereShape>(new SphereShape(radius));
    dart::dynamics::Inertia inertia;
    inertia.setMass(mass);
    inertia.setMoment(shape->computeInertia(mass));

    FreeJoint::Properties prop;
    prop.mT_ParentBodyToJoint.setIdentity();
    prop.mT_ChildBodyToJoint.setIdentity();

    BodyNodePtr bn = skel->createJointAndBodyNodePair<FreeJoint>(
      nullptr,prop,BodyNode::AspectProperties("ball")).second;

    auto sn = bn->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(shape);
    bn->setCollidable(false);
    bn->setInertia(inertia);
}
void
DrawSkeleton(const dart::dynamics::SkeletonPtr& skel,const Eigen::Vector3d& color)
{
    glDisable(GL_DEPTH_TEST);

    for(int j =0;j<skel->getNumBodyNodes();j++)
    {
        auto bn = skel->getBodyNode(j);
        auto shapeNodes = bn->getShapeNodesWith<VisualAspect>();
        auto T = bn->getTransform();
        for (auto& sn : shapeNodes)
        {
            auto shape = sn->getShape().get();

            if (shape->is<BoxShape>())
            {
                const auto* box = static_cast<const BoxShape*>(shape);
                const auto& size = box->getSize();
                DrawCapsule(T,size[0],size[1],color);
            }
            else if (shape->is<SphereShape>())
            {
                const auto* sphere = static_cast<const SphereShape*>(shape);
                const auto& r = sphere->getRadius();
                DrawSphere(T,r,color);

            }
            
        }

    }
    for(int j =0;j<skel->getNumBodyNodes();j++)
    {
        auto bn = skel->getBodyNode(j);
        auto T = bn->getTransform();
        Eigen::Vector3d p = T.translation();
        Eigen::Vector3d v = bn->getCOMLinearVelocity();
        v*= 0.05;
        DrawArrow(p.block<2,1>(0,0),v.block<2,1>(0,0),Eigen::Vector3d(0,0,1));
    }
    glEnable(GL_DEPTH_TEST);
}
