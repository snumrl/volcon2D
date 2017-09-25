#ifndef __FSM_H__
#define __FSM_H__
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <Eigen/Core>

class State
{
private:
	std::map<std::string,State*>	mEvents;
	double							mTimeElapsed;
	std::vector<std::pair<double,Eigen::VectorXd>>	mMotion;
public:
	State(const std::string& motion_path);
	void 	AddEvent(const std::string& name,State* next_state);
	State*	GetNextState(const std::string& event_name);
	std::map<std::string,State*>& GetEvents();
	std::string GetMotion(Eigen::VectorXd& motion);

	const double& GetTime();
	void SetTime(const double& t);
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

	void 	SetCurrentState(State* s);
	State* 	GetCurrentState();
	std::map<std::string,State*>& GetStates();
	void	Trigger(const std::string& event_name);
	void GetMotion(Eigen::VectorXd& motion);

};



#endif