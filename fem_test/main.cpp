#include "SimulationWindow2D.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "fem2D/World.h"
#include "GL/glut.h"


int main(int argc,char** argv)
{
	SimulationWindow2D simwindow;
	glutInit(&argc, argv);
	simwindow.InitWindow(800,800,"vmcon2D");
	glutMainLoop();
}
