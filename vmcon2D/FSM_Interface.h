#ifndef __FSM_INTERFACE_H__
#define __FSM_INTERFACE_H__
#include <string>
class Machine;

void DrawMachine(Machine* machine,double x,double y);

void MakeMachine(const std::string& file_path,Machine* machine);
#endif