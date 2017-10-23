#include "SimulationWindow2D.h"
#include "fem2D/Mesh/MeshHeaders.h"
#include "fem2D/Constraint/ConstraintHeaders.h"
#include "fem2D/World.h"
//#include <omp.h>
#include "GL/glut.h"


int main(int argc,char** argv)
{
	//omp_set_num_threads(8);
	//Eigen::setNbThreads(8);
	// FEM::World* soft_world = new FEM::World(
	// 	// FEM::IntegrationMethod::IMPLICIT_NEWTON_METHOD,		//Integration Method
	// 	// FEM::IntegrationMethod::QUASI_STATIC,		//Integration Method
	// 	// FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,		//Integration Method
	// 	FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,		//Integration Method
	// 	1.0/120.0,										//time_step
	// 	100, 											//max_iteration	
	// 	0.999											//damping_coeff
	// 	);

	// double youngs_modulus = 1E4;
	// double poisson_ratio = 0.3;
	// double muscle_stiffness =1E5;
	
	// FEM::RectangleMesh *dm = new FEM::RectangleMesh(0.1,1,2,20);
	
	// //Initialize Vertices
	// const auto& X = dm->GetVertices();
	// Eigen::VectorXd v(X.size()*2);

	// for(int i =0;i<X.size();i++)
	// 	v.block<2,1>(i*2,0) = X[i];

	// //Initialize Constraints
	// std::vector<FEM::Constraint*> constraint_vector;

	// for(const auto& tri : dm->GetTriangles())
	// {
	// 	int i0,i1,i2;
	// 	Eigen::Vector2d p0,p1,p2;
		
	// 	i0 = tri[0];
	// 	i1 = tri[1];
	// 	i2 = tri[2];
	// 	p0 = X[i0];
	// 	p1 = X[i1];
	// 	p2 = X[i2];

	// 	Eigen::Matrix2d Dm;

	// 	Dm.block<2,1>(0,0) = p1 - p0;
	// 	Dm.block<2,1>(0,1) = p2 - p0;
	// 	if(Dm.determinant()<0)
	// 	{
	// 		i0 = tri[0];
	// 		i1 = tri[2];
	// 		i2 = tri[1];
	// 		p0 = X[i0];
	// 		p1 = X[i1];
	// 		p2 = X[i2];
	// 		Dm.block<2,1>(0,0) = p1 - p0;
	// 		Dm.block<2,1>(0,1) = p2 - p0;

	// 		constraint_vector.push_back(new FEM::CorotateFEMConstraint(youngs_modulus,poisson_ratio,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));
	// 		constraint_vector.push_back(new FEM::LinearMuscleConstraint(muscle_stiffness,Eigen::Vector2d::UnitY(),i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
	// 		//constraint_vector.push_back(new FEM::LinearMuscleConstraint(muscle_stiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
	// 	}
	// }

	// constraint_vector.push_back(new FEM::AttachmentConstraint(1000000,0,Eigen::Vector2d(-0.05,-0.5)));
	// constraint_vector.push_back(new FEM::AttachmentConstraint(1000000,41,Eigen::Vector2d(0,-0.5)));
	// constraint_vector.push_back(new FEM::AttachmentConstraint(1000000,82,Eigen::Vector2d(0.05,-0.5)));
	
	//Initialize world
	// soft_world->AddBody(v, constraint_vector,0.1);
	// soft_world->Initialize();
	//SimulationWindow2D simwindow(soft_world);
	SimulationWindow2D simwindow;
	glutInit(&argc, argv);
	simwindow.InitWindow(800,800,"vmcon2D");
	glutMainLoop();
}
