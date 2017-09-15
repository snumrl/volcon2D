#ifndef __FSM_INTERFACE_H__
#define __FSM_INTERFACE_H__
#include <string>
#include <map>
#include <Eigen/Core>
class Machine;
class State;
void DrawMachine(Machine* machine,double x,double y,std::map<State*,Eigen::Vector2d> state_positions = std::map<State*,Eigen::Vector2d>());

void MakeMachine(const std::string& file_path,Machine* machine);
#endif