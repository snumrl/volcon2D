#include "FSM_Interface.h"
#include "FSM.h"
#include "GUI/GL_function.h"
#include <tinyxml.h>
#include "GL/glut.h"

void
DrawMachine(Machine* machine,double x,double y)
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
    ;
    std::map<State*,Eigen::Vector3d> state_position;
    Eigen::Isometry3d T;
    T.setIdentity();


    for(auto& state_pair : states)
    {
    	auto state = state_pair.second;

    	state_position.insert(std::make_pair(state,Eigen::Vector3d(x + Radius*cos(phi),x+Radius*sin(phi),0)));
    	T.translation() = state_position[state];

        if(machine->GetCurrentState() == state)
        {
            DrawSphere(T,radius,Eigen::Vector3d(0.9,0.5,0.5));    
            DrawStringOnScreen(x + Radius*cos(phi)-0.5*fontsize*state_pair.first.size(),x+Radius*sin(phi)-0.5*fontsize,state_pair.first,false,Eigen::Vector3d(0,0,0));

        }
        else
        {
            DrawSphere(T,radius,Eigen::Vector3d(0.9,0.9,0.9));    
            DrawStringOnScreen(x + Radius*cos(phi)-0.5*fontsize*state_pair.first.size(),x+Radius*sin(phi)-0.5*fontsize,state_pair.first,false,Eigen::Vector3d(0,0,0));
        }

    	phi += 2*3.141592/(double)slice;
    }

	for(auto& state_pair : states)
    {
    	auto state = state_pair.second;

        auto& event = state->GetEvents();
        Eigen::Vector3d p0 = state_position[state];
        for(auto& next_state : event)
        {
            Eigen::Vector3d p1 = state_position[next_state.second];

            Eigen::Vector2d p,v;
            p[0] = p0[0];
            p[1] = p0[1];
            v[0] = p1[0]-p[0];
            v[1] = p1[1]-p[1];

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