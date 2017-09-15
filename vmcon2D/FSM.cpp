#include "FSM.h"
#include <iostream>
void
State::
AddEvent(const std::string& name,State* next_state)
{
	mEvents.insert(std::make_pair(name,next_state));	
}

State*
State::
GetNextState(const std::string& event_name)
{
	if(mEvents.find(event_name)==mEvents.end())
	{
		std::cout<<"NO event name : "<<event_name<<std::endl;
		return nullptr;
	}
	return mEvents.at(event_name);
}
std::map<std::string,State*>&
State::
GetEvents()
{
	return mEvents;
}
Machine::
Machine()
	:mCurrentState(nullptr)
{

}

State*
Machine::
AddState(const std::string& name)
{
	mStates.insert(std::make_pair(name,new State()));
	if(mStates.size()==1)
		mCurrentState = mStates.at(name);
		

	return mStates.at(name);
}
void 	
Machine::
AddEvent(State* state_from,State* state_to,const std::string& name)
{
	state_from->AddEvent(name,state_to);
}


State* 	
Machine::
GetCurrentState()
{
	return mCurrentState;
}
std::map<std::string,State*>&
Machine::
GetStates()
{
	return mStates;
}
void	
Machine::
Trigger(const std::string& name)
{
	State* next_state;
	next_state = mCurrentState->GetNextState(name);
	if(next_state!=nullptr)
		mCurrentState = next_state;
}