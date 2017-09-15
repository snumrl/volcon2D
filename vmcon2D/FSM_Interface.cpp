#include "FSM_Interface.h"
#include "FSM.h"
#include "GUI/GL_function.h"
#include <tinyxml.h>
#include "GL/glut.h"

void
DrawMachine(Machine* machine,double x,double y,std::map<State*,Eigen::Vector2d> state_positions)
{
	   // draws text on the screen
    GLint oldMode;
    glGetIntegerv(GL_MATRIX_MODE, &oldMode);
    glMatrixMode(GL_PROJECTION);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    auto& states = machine->GetStates();

    double radius = 0.03;
    double Radius = 0.1;
    double fontsize = 0.007;
    int slice = states.size();
    double phi = 0.0;
    Eigen::Isometry3d T;
    T.setIdentity();

    bool isStatePostionsNull = false;
    if(state_positions.size()==0)
        isStatePostionsNull = true;
    for(auto& state_pair : states)
    {
    	auto state = state_pair.second;

        if(isStatePostionsNull)
    	   state_positions.insert(std::make_pair(state,Eigen::Vector2d(Radius*cos(phi),Radius*sin(phi))));

        T.translation()[0] = x+state_positions[state][0];
    	T.translation()[1] = y+state_positions[state][1];

        if(machine->GetCurrentState() == state)
        {

            DrawSphere(T,radius,Eigen::Vector3d(0.9,0.5,0.5));    
            DrawStringOnScreen(x+state_positions[state][0]-0.5*fontsize*state_pair.first.size(),y+state_positions[state][1]-0.5*fontsize,state_pair.first,false,Eigen::Vector3d(0,0,0));

        }
        else
        {
            DrawSphere(T,radius,Eigen::Vector3d(0.9,0.9,0.9));    
            DrawStringOnScreen(x+state_positions[state][0]-0.5*fontsize*state_pair.first.size(),y+state_positions[state][1]-0.5*fontsize,state_pair.first,false,Eigen::Vector3d(0,0,0));
        }

    	phi += 2*3.141592/(double)slice;
    }

	for(auto& state_pair : states)
    {
    	auto state = state_pair.second;

        auto& event = state->GetEvents();
        Eigen::Vector2d p0 = state_positions[state];
        for(auto& next_state : event)
        {
            Eigen::Vector2d p1 = state_positions[next_state.second];

            Eigen::Vector2d p,v;
            p[0] = p0[0]+x;
            p[1] = p0[1]+y;
            v[0] = p1[0]-p[0]+x;
            v[1] = p1[1]-p[1]+y;

            p += radius*v.normalized();
            v -= 2*radius*(v.normalized()).eval();

            Eigen::Vector2d m = p+0.5*v;
            DrawArrow(p,v,Eigen::Vector3d(0,0,0));
            // DrawStringOnScreen(m[0]-0.5*fontsize*next_state.first.size(),m[1]-0.5*fontsize,next_state.first,false);


        }
    }
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(oldMode);
}

void
MakeMachine(const std::string& path,Machine* machine)
{
    TiXmlDocument doc;
    if(!doc.LoadFile(path))
    {
        std::cout<<"Cant open XML file : "<<path<<std::endl;
        return;
    }

    TiXmlElement* machine_xml = doc.FirstChildElement("Machine");

    for(TiXmlElement* state_xml = machine_xml->FirstChildElement("State");state_xml!=nullptr;state_xml = state_xml->NextSiblingElement("State"))
        machine->AddState(state_xml->Attribute("name"));
    auto& states = machine->GetStates();
    for(TiXmlElement* event_xml = machine_xml->FirstChildElement("Event");event_xml!=nullptr;event_xml = event_xml->NextSiblingElement("Event"))
        machine->AddEvent(states[event_xml->Attribute("from")],states[event_xml->Attribute("to")],event_xml->Attribute("name"));
}