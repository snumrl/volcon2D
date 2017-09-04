#ifndef __FEM2D_INTERFACE_H__
#define __FEM2D_INTERFACE_H__
#include <vector>
#include <Eigen/Core>
#include "fem2D/Constraint/ConstraintHeaders.h"
class Muscle;

void DrawConstraint(FEM::Constraint* c,const Eigen::VectorXd& x);
void DrawConstraints(std::vector<FEM::Constraint*>& cs,const Eigen::VectorXd& x);

#endif