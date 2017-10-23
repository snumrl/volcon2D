#include "fem2D/World.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "MusculoSkeletalSystem.h"
#include "Controller.h"
#include "FSM.h"
#include "Record.h"

using namespace dart::dynamics;
using namespace dart::simulation;
void
Record::
Set(const dart::simulation::WorldPtr& rigid_world,
	FEM::World* soft_world,
	MusculoSkeletalSystem* musculo_skeletal_system,
	Controller* controller)
{
	t = rigid_world->getTime();
	rigid_body_positions.resize(rigid_world->getNumSkeletons());
	rigid_body_velocities.resize(rigid_world->getNumSkeletons());
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_body_positions[i] = rigid_world->getSkeleton(i)->getPositions();

	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_body_velocities[i] = rigid_world->getSkeleton(i)->getVelocities();

	soft_body_positions = soft_world->GetPositions();
	activation_levels = musculo_skeletal_system->GetActivationLevel();
	for(auto& muscle : musculo_skeletal_system->GetMuscles())
		muscle_forces.push_back(std::make_pair(muscle->force_origin,muscle->force_insertion));
	if(controller!=nullptr)
		state = controller->GetMachine()->GetCurrentState();
}

void
Record::
Get(const dart::simulation::WorldPtr& rigid_world,
	FEM::World* soft_world,
	MusculoSkeletalSystem* musculo_skeletal_system,
	Controller* controller)
{
	rigid_world->setTime(t);
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_world->getSkeleton(i)->setPositions(rigid_body_positions[i]);
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_world->getSkeleton(i)->setVelocities(rigid_body_velocities[i]);
	
	soft_world->SetTime(t);
	soft_world->SetPositions(soft_body_positions);
	int count = 0;
    for(auto& muscle : musculo_skeletal_system->GetMuscles())
    {
        muscle->activationLevel = activation_levels[count];
        for(auto& mc : muscle->muscleConstraints)
            mc->SetActivationLevel(muscle->activationLevel);
        muscle->force_origin = muscle_forces[count].first;
        muscle->force_insertion = muscle_forces[count].second;
        muscle->origin->GetP() = GetPoint(muscle->originWayPoints[0]);
        muscle->insertion->GetP() = GetPoint(muscle->insertionWayPoints[0]);
        count++;
    }
    if(controller!=nullptr)
    	controller->GetMachine()->SetCurrentState(state);
}

