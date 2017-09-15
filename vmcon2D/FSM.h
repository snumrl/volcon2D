#ifndef __FSM_H__
#define __FSM_H__
#include <vector>
#include <string>
#include <map>
class State
{
private:
	std::map<std::string,State*>	mEvents;

public:
	void 	AddEvent(const std::string& name,State* next_state);
	State*	GetNextState(const std::string& event_name);
	std::map<std::string,State*>& GetEvents();
};
class Machine
{
private:
	State*					mCurrentState;

	std::map<std::string,State*>		mStates;
public:
	Machine();
	
	State* 	AddState(const std::string& name);
	void 	AddEvent(State* state_from,State* state_to,const std::string& name);

	State* 	GetCurrentState();
	std::map<std::string,State*>& GetStates();
	void	Trigger(const std::string& event_name);

};



#endif