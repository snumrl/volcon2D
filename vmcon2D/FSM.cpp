#include "FSM.h"
#include <iostream>
#include <fstream>
#include <string>

void ReadMotion(const std::string path,std::vector<std::pair<double,Eigen::VectorXd>>& motion)
{
	std::ifstream ifs(path);
	std::string str;
	std::stringstream ss;
	std::string t;
	std::string token;
	std::vector<std::string> token_vector;
	while(!ifs.eof())
	{
		str.clear();
		ss.clear();
		token_vector.clear();

		std::getline(ifs,str);
		ss.str(str);
		if(ss>>t){
			while(!ss.eof())
			{
				ss>>token;
				token_vector.push_back(token);
			}
			Eigen::VectorXd vec(token_vector.size());
			for(int i =0;i<token_vector.size();i++)
				vec[i] = std::stod(token_vector[i]);

			motion.push_back(std::make_pair(std::stod(t),vec));
		}
	}
	// std::cout<<path<<std::endl;
	// for(int i =0;i<motion.size();i++)
	// {
	// 	std::cout<<"time : "<<motion[i].first<<std::endl<<motion[i].second.transpose()<<std::endl;
	// }
}
State::
State(const std::string& motion_path)
{
	std::string path_to_motion_dir = "../vmcon2D/export/motions/";
	
	ReadMotion(path_to_motion_dir+motion_path,mMotion);
}
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
		// std::cout<<"NO event name : "<<event_name<<std::endl;
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

std::string
State::
GetMotion(Eigen::VectorXd& motion)
{
	int k =0,k1 =0;
	for(int i =0;i<mMotion.size();i++)
	{
		if(mMotion[i].first<mTimeElapsed)
			k=i;
	}

	if(k==mMotion.size()-1)
		return std::string("end");

	k1 = k+1;
	double t = mTimeElapsed-mMotion[k].first;
	double dt = mMotion[k1].first-mMotion[k].first;

	t/= dt;
	motion = (1.0-t)*(mMotion[k].second) + (t)*(mMotion[k1].second);

	return std::string("no_event");
}











const double&
State::
GetTime()
{
	return mTimeElapsed;
}
void
State::
SetTime(const double& t)
{
	mTimeElapsed = t;
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
	mStates.insert(std::make_pair(name,new State(name)));
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


void
Machine::
SetCurrentState(State* s)
{
	mCurrentState = s;
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
	if(next_state!=nullptr){
		mCurrentState = next_state;
		mCurrentState->SetTime(0.0);
	}
	
}

void
Machine::
GetMotion(Eigen::VectorXd& ret)
{
	std::string event = mCurrentState->GetMotion(ret);	
	Trigger(event);
}