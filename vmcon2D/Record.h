#ifndef __RECORD_H__
#define __RECORD_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include <Eigen/Core>
#include <Eigen/StdVector>
class MusculoSkeletalSystem;
class Controller;
namespace FEM
{
	class World;
};
class Record
{
public:
	
	void Set(const dart::simulation::WorldPtr& rigid_world,
			 FEM::World* soft_world,
			 MusculoSkeletalSystem* musculo_skeletal_system,
			 Controller* controller=nullptr);
	void Get(const dart::simulation::WorldPtr& rigid_world,
			 FEM::World* soft_world,
			 MusculoSkeletalSystem* musculo_skeletal_system,
			 Controller* controller=nullptr);
	double			t;
	std::vector<Eigen::VectorXd> rigid_body_positions;
	std::vector<Eigen::VectorXd> rigid_body_velocities;
	Eigen::VectorXd soft_body_positions;
	Eigen::VectorXd activation_levels;
	std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>> muscle_forces;

	Eigen::VectorXd target_positions;
	Eigen::VectorXd target_positions2;
	Eigen::VectorXd target_velocities;
};

#endif