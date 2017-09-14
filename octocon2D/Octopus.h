#ifndef __OCTOPUS_H__
#define __OCTOPUS_H__
#include <vector>
#include <string>
#include <Eigen/Core>
struct Muscle
{
	std::vector<FEM::MuscleConstraint*>		constraints;
	double									activationLevel;
};

class Octopus
{
private:
	FEM::Mesh*								mMesh;
	std::vector<Muscle*>					mMuscles;
	Eigen::VectorXd							mActivationLevels;
public:
	Octopus(const std::string& mesh_file);

	void AddMuscle(const std::vector<int>& element_index_list,const Eigen::Vector3d& fiber_direction);
	void Initialize(FEM::World* world);
};



#endif