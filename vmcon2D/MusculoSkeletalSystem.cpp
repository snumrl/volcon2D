#include "MusculoSkeletalSystem.h"
#include "DART_Interface.h"

#include <tinyxml.h>

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace FEM;
Eigen::Vector2d GetPoint(const AnchorPoint& p)
{
	Eigen::Vector3d point =  p.first->getTransform()*p.second;

	return Eigen::Vector2d(point[0],point[1]);
}
void
Muscle::
TransferForce(Eigen::Vector2d& f_origin,Eigen::Vector2d& f_insertion)
{
	int no = originWayPoints.size();
	int ni = insertionWayPoints.size();

	if(no>1)
	{
		Eigen::Vector2d u = (GetPoint(insertionWayPoints[0])-GetPoint(originWayPoints[0])).normalized();
		Eigen::Vector2d v = (GetPoint(originWayPoints[no-2])-GetPoint(originWayPoints[no-1])).normalized();
		double angle = acos(u.dot(v));
		double sign = u[0]*v[1] - u[1]*v[0];
		if(sign<0)
			angle = -angle;
		Eigen::Rotation2D<double> R(angle);
		f_origin = R*f_origin;
	}

	if(ni>1)
	{
		Eigen::Vector2d u = (GetPoint(originWayPoints[0])-GetPoint(insertionWayPoints[0])).normalized();
		Eigen::Vector2d v = (GetPoint(insertionWayPoints[ni-2])-GetPoint(insertionWayPoints[ni-1])).normalized();
		double angle = acos(u.dot(v));
		double sign = u[0]*v[1] - u[1]*v[0];
		if(sign<0)
			angle = -angle;
		Eigen::Rotation2D<double> R(angle);
		f_insertion = R*f_insertion;
	}
}

MusculoSkeletalSystem::
MusculoSkeletalSystem()
	:mTendonStiffness(1E5),mMuscleStiffness(1E5),mYoungsModulus(5E5),mPoissonRatio(0.3)
{

}
void
MusculoSkeletalSystem::
AddMuscle(
	const std::vector<AnchorPoint>& origin,
	const int&			origin_i,
	const std::vector<AnchorPoint>& insertion,
	const int&			insertion_i,
	const Eigen::Vector2d& fiber_direction,
	Mesh* mesh)
{
	mMuscles.push_back(new Muscle());
	auto muscle = mMuscles.back();
	muscle->mesh = mesh;
	muscle->originWayPoints = origin;
	muscle->insertionWayPoints = insertion;
	std::vector<Eigen::Vector2d> p_origin,p_insertion;
	double l0_origin=0,l0_insertion=0;
	const auto& vertices = muscle->mesh->GetVertices();
	const auto& triangles = muscle->mesh->GetTriangles();
	for(int i=0;i<origin.size();i++)
		p_origin.push_back(GetPoint(origin[i]));

	for(int i=0;i<p_origin.size()-1;i++)
		l0_origin += (p_origin[i] - p_origin[i+1]).norm();
	l0_origin += (p_origin[0]- vertices[origin_i]).norm();
	for(int i=0;i<insertion.size();i++)
		p_insertion.push_back(GetPoint(insertion[i]));
	
	for(int i=0;i<p_insertion.size()-1;i++)
		l0_insertion += (p_insertion[i] - p_insertion[i+1]).norm();
	l0_insertion += (p_insertion[0]- vertices[insertion_i]).norm();


	for(const auto& tri : triangles)
	{
		int i0,i1,i2;
		Eigen::Vector2d p0,p1,p2;
		
		i0 = tri[0];
		i1 = tri[1];
		i2 = tri[2];
		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];

		Eigen::Matrix2d Dm;

		Dm.block<2,1>(0,0) = p1 - p0;
		Dm.block<2,1>(0,1) = p2 - p0;
		if(Dm.determinant()<0)
		{
			i0 = tri[0];
			i1 = tri[2];
			i2 = tri[1];
			p0 = vertices[i0];
			p1 = vertices[i1];
			p2 = vertices[i2];
			Dm.block<2,1>(0,0) = p1 - p0;
			Dm.block<2,1>(0,1) = p2 - p0;

			muscle->constraints.push_back(new CorotateFEMConstraint(mYoungsModulus,mPoissonRatio,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));

			// muscle->constraints.push_back(new SpringConstraint(10000.0,i0,i1,(p0-p1).norm()));
			// muscle->constraints.push_back(new SpringConstraint(10000.0,i1,i2,(p1-p2).norm()));
			// muscle->constraints.push_back(new SpringConstraint(10000.0,i2,i0,(p2-p0).norm()));
			muscle->muscleConstraints.push_back(new LinearMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
			// muscle->muscleConstraints.push_back(new HillTypeMuscleConstraint(mMuscleStiffness,fiber_direction,i0,i1,i2,1.0/6.0*(Dm.determinant()),Dm.inverse()));	
		}
		
	}
	muscle->origin = new AttachmentConstraint(mTendonStiffness,origin_i,p_origin[0]);
	muscle->insertion = new AttachmentConstraint(mTendonStiffness,insertion_i,p_insertion[0]);
	mAttachementConstraintVector.push_back(muscle->origin);
	mAttachementConstraintVector.push_back(muscle->insertion);
	muscle->activationLevel = 0.0;
}

void
MusculoSkeletalSystem::
Initialize(FEM::World* world)
{
	int offset = world->GetNumVertices();
	for(int i=0;i<mMuscles.size();i++)
	{
		Muscle* muscle = mMuscles[i];

		const std::vector<Eigen::Vector2d>& vertices = muscle->mesh->GetVertices();
		muscle->origin->AddOffset(offset);
		muscle->insertion->AddOffset(offset);

		for(auto& c : muscle->constraints)
			c->AddOffset(offset);
		for(auto& c : muscle->muscleConstraints)
			c->AddOffset(offset);
		
		Eigen::VectorXd v(vertices.size()*2);
		for(int i =0;i<vertices.size();i++)
			v.block<2,1>(i*2,0) = vertices[i];

		world->AddBody(v,muscle->constraints,1.0);
		for(auto& c: muscle->muscleConstraints)
			world->AddConstraint(c);

		world->AddConstraint(muscle->origin);
		world->AddConstraint(muscle->insertion);

		offset += vertices.size();
	}

	mActivationLevel.resize(mMuscles.size());
	mActivationLevel.setZero();

}

void
MusculoSkeletalSystem::
TransformAttachmentPoints()
{
	for(auto& muscle : mMuscles)
	{
		auto& origin_way_points = muscle->originWayPoints;
		auto& insertion_way_points = muscle->insertionWayPoints;

		auto po = GetPoint(origin_way_points[0]);
		auto pi = GetPoint(insertion_way_points[0]);
		muscle->origin->GetP() = po;
		muscle->insertion->GetP() = pi;
	}
}
void
MusculoSkeletalSystem::
SetActivationLevel(const Eigen::VectorXd& a)
{
	mActivationLevel = a;
	for(int i=0;i<mMuscles.size();i++)
		for(auto& mc : mMuscles[i]->muscleConstraints)
			mc->SetActivationLevel(a[i]);
}
void
MusculoSkeletalSystem::
ApplyForcesToSkeletons(FEM::World* world)
{
	Eigen::VectorXd X = world->GetPositions();
	Eigen::VectorXd force_origin(X.rows()),force_insertion(X.rows());
	Eigen::Vector2d fo,fi;
	
	for(auto& muscle : mMuscles)
	{
		auto& origin_way_points = muscle->originWayPoints;
		auto& insertion_way_points = muscle->insertionWayPoints;
		int no = origin_way_points.size();
		int ni = insertion_way_points.size();

		force_origin.setZero();
		force_insertion.setZero();

		muscle->origin->EvalGradient(X,force_origin);
		muscle->insertion->EvalGradient(X,force_insertion);
		
		fo = force_origin.block<2,1>(muscle->origin->GetI0()*2,0);
		fi = force_insertion.block<2,1>(muscle->insertion->GetI0()*2,0);
		
		muscle->TransferForce(fo,fi);

		Eigen::Vector3d f_origin_3D(fo[0],fo[1],0.0);
		Eigen::Vector3d f_insertion_3D(fi[0],fi[1],0.0);

		origin_way_points[no-1].first->addExtForce(f_origin_3D,origin_way_points[no-1].second);
		insertion_way_points[ni -1].first->addExtForce(f_insertion_3D,insertion_way_points[ni-1].second);

		muscle->force_origin =fo;
		muscle->force_insertion =fi;
	}	
}



Eigen::MatrixXd
MusculoSkeletalSystem::
ComputeForceDerivative(FEM::World* world)
{
	Eigen::VectorXd X = world->GetPositions();
	Eigen::MatrixXd J(mMuscles.size()*4,mMuscles.size());
	J.setZero();
	for(int i=0;i<mMuscles.size();i++)
	{
		auto& muscle = mMuscles[i];
		J.block(0,i,mMuscles.size()*4,1) = world->ComputeJacobian(muscle->muscleConstraints,mAttachementConstraintVector);

		for(int j=0;j<mMuscles.size();j++){
			Eigen::Vector2d fo,fi;
			fo = J.block<2,1>(j*4+0,i);
			fi = J.block<2,1>(j*4+2,i);
			mMuscles[j]->TransferForce(fo,fi);
			J.block<2,1>(j*4+0,i) = fo;
			J.block<2,1>(j*4+2,i) = fi;
		}
	}
	// for(int i =0;i<J.rows();i++)
	// 	for(int j =0;j<J.cols();j++)
	// 		if(fabs(J(i,j))<1E-4)
	// 			J(i,j) =0.0;				

	return J;
}
Eigen::VectorXd
MusculoSkeletalSystem::
ComputeForce(FEM::World* world)
{
	Eigen::VectorXd X = world->GetPositions();	
	Eigen::VectorXd b(mMuscles.size()*4);

	Eigen::VectorXd force_origin(X.rows()),force_insertion(X.rows());
	Eigen::Vector2d fo,fi;
	for(int i=0;i<mMuscles.size();i++)
	{
		auto& muscle = mMuscles[i];
		auto& origin_way_points = muscle->originWayPoints;
		auto& insertion_way_points = muscle->insertionWayPoints;

		force_origin.setZero();
		force_insertion.setZero();

		muscle->origin->EvalGradient(X,force_origin);
		muscle->insertion->EvalGradient(X,force_insertion);
		
		fo = force_origin.block<2,1>(muscle->origin->GetI0()*2,0);
		fi = force_insertion.block<2,1>(muscle->insertion->GetI0()*2,0);
		
		muscle->TransferForce(fo,fi);

		b.block<2,1>(4*i+0,0) = fo;
		b.block<2,1>(4*i+2,0) = fi;
	}
	// for(int i =0;i<b.rows();i++)
	// 	if(fabs(b[i])<1E-4)
	// 		b[i] = 0.0;

	return b;
}
void
MakeSkeleton(MusculoSkeletalSystem* ms)
{
	ms->GetSkeleton() = Skeleton::create("human");
	auto& skel = ms->GetSkeleton();
	MakeRootBody(skel,"Torso",Eigen::Vector3d(0.05,0.6,0.0),Eigen::Vector3d(0,-0.3,0),5);

	MakeBody(skel,skel->getBodyNode("Torso"),"NeckR",
		Eigen::Vector3d(0.3,0.05,0.0),
		Eigen::Vector3d(0.0,0.3,0),
		Eigen::Vector3d(-0.15,0.0,0),2);

	// MakeBody(skel,skel->getBodyNode("Torso"),"NeckL",
	// 	Eigen::Vector3d(0.3,0.05,0.0),
	// 	Eigen::Vector3d(0.0,0.3,0),
	// 	Eigen::Vector3d(0.15,0.0,0),5);

	MakeBody(skel,skel->getBodyNode("NeckR"),"ShoulderR",
		Eigen::Vector3d(0.3,0.05,0.0),
		Eigen::Vector3d(0.15,0.0,0),
		Eigen::Vector3d(-0.15,0.0,0),1);

	// MakeBody(skel,skel->getBodyNode("NeckL"),"ShoulderL",
	// 	Eigen::Vector3d(0.3,0.05,0.0),
	// 	Eigen::Vector3d(-0.15,0.0,0),
	// 	Eigen::Vector3d(0.15,0.0,0),5);

	MakeBody(skel,skel->getBodyNode("ShoulderR"),"ElbowR",
		Eigen::Vector3d(0.3,0.05,0.0),
		Eigen::Vector3d(0.15,0.0,0),
		Eigen::Vector3d(-0.15,0.0,0),1);

	// MakeBody(skel,skel->getBodyNode("ShoulderL"),"ElbowL",
	// 	Eigen::Vector3d(0.3,0.05,0.0),
	// 	Eigen::Vector3d(-0.15,0.0,0),
	// 	Eigen::Vector3d(0.15,0.0,0),5);

	MakeWeldBody(skel,skel->getBodyNode("Torso"),"Head",
		0.07,
		Eigen::Vector3d(0,0.40,0),
		Eigen::Vector3d(0,0,0),
		1);

	MakeWeldBody(skel,skel->getBodyNode("ElbowR"),"HandR",
		0.02,
		Eigen::Vector3d(0.17,0,0),
		Eigen::Vector3d(0,0,0),
		0.5);

	// MakeWeldBody(skel,skel->getBodyNode("ElbowL"),"HandL",
	// 	0.02,
	// 	Eigen::Vector3d(-0.17,0.0,0),
	// 	Eigen::Vector3d(0,0,0),
	// 	3);



	auto pos = skel->getPositions();
	// pos[0] = 0.0;

	pos[0] = 0.1;

	pos[1] = -1.0;

	pos[2] = -1.0;

	// 	pos[1] = 0.1;
	// pos[2] = -0.1;

	// pos[3] = -1.0;
	// pos[4] = 1.0;

	// pos[5] = -1.0;
	// pos[6] = 1.0;

	
	skel->setPositions(pos);
	skel->computeForwardKinematics(true,false,false);
	// skel->getDof(0)->setPositionLimits(-0.0,0.0);
	// skel->getDof(1)->setPositionLimits(0.0,0.0);
	// skel->getDof(2)->setPositionLimits(-0.0,0.0);
	// skel->getDof(3)->setPositionLimits(0.0,0.0);
	// skel->getDof(4)->setPositionLimits(0.0,0.0);

	// skel->getDof(0)->setPositionLimits(-0.1,0.1);
	// skel->getDof(1)->setPositionLimits(-0.4,0.2);
	// skel->getDof(2)->setPositionLimits(-0.2,0.4);
	// skel->getDof(3)->setPositionLimits(-1.57,0.0);
	// skel->getDof(4)->setPositionLimits(0.0,1.57);

	// skel->getDof(5)->setPositionLimits(-2.0,2.0);
	// skel->getDof(6)->setPositionLimits(-2.0,2.0);
	
	// skel->getDof(0)->setPositionLimits(0.0,0.0);
	// skel->getDof(1)->setPositionLimits(-0.1,0.2);
	// // skel->getDof(2)->setPositionLimits(-0.2,0.0);
	// skel->getDof(2)->setPositionLimits(-1.0,-1.0);
	// // skel->getDof(4)->setPositionLimits(0.0,1.57);
	// skel->getDof(3)->setPositionLimits(-1.0,-1.0);
	// // skel->getDof(6)->setPositionLimits(-2.0,2.0);
	
	for(int i =0;i<skel->getNumDofs();i++)
		skel->getDof(i)->getJoint()->setPositionLimitEnforced(true);
	for(int i=0;i<skel->getNumBodyNodes();i++)
		skel->getBodyNode(i)->setCollidable(false);
}
void
MakeMuscles(const std::string& path,MusculoSkeletalSystem* ms)
{
    auto& skel = ms->GetSkeleton();

    TiXmlDocument doc;
    if(!doc.LoadFile(path))
    {
        std::cout<<"Cant open XML file : "<<path<<std::endl;
        return;
    }

    TiXmlElement* muscles = doc.FirstChildElement("Muscles");

    for(TiXmlElement* unit = muscles->FirstChildElement("unit");unit!=nullptr;unit = unit->NextSiblingElement("unit"))
    {
        TiXmlElement* ori = unit->FirstChildElement("origin");
        std::vector<AnchorPoint> p_ori,p_ins;
       
        for(TiXmlElement* anc = ori->FirstChildElement("anchor");anc!=nullptr;anc = anc->NextSiblingElement("anchor"))   
        {
            std::string body_name = anc->Attribute("body");
            double x = std::stod(anc->Attribute("x"));
            double y = std::stod(anc->Attribute("y"));
            p_ori.push_back(AnchorPoint(skel->getBodyNode(body_name.c_str()),Eigen::Vector3d(x,y,0.0)));
        }
        
        TiXmlElement* ins = unit->FirstChildElement("insertion");
        for(TiXmlElement* anc = ins->FirstChildElement("anchor");anc!=nullptr;anc = anc->NextSiblingElement("anchor"))   
        {
            std::string body_name = anc->Attribute("body");
            double x = std::stod(anc->Attribute("x"));
            double y = std::stod(anc->Attribute("y"));
            p_ins.push_back(AnchorPoint(skel->getBodyNode(body_name.c_str()),Eigen::Vector3d(x,y,0.0)));
        }

        Eigen::Vector2d muscle_start,muscle_end;

        muscle_start = GetPoint(p_ori[0]);
        muscle_end = GetPoint(p_ins[0]);

        // p_ori.erase(p_ori.begin());
        // p_ins.erase(p_ins.begin());

        double len = (muscle_start - muscle_end).norm();
        Eigen::Vector2d unit_dir = (muscle_start - muscle_end).normalized();
        double cosa = unit_dir[0];
        double sina = unit_dir[1];
        double angle = atan2(sina,cosa);

        

        TiXmlElement* mesh_element = unit->FirstChildElement("mesh");
        std::string mesh_type = mesh_element->Attribute("type");

        Eigen::Affine2d T = Eigen::Affine2d::Identity();
        Eigen::Rotation2D<double> R(angle);
        T.translation() = 0.5*(muscle_start + muscle_end);
        T.linear()      = R*Eigen::Scaling(Eigen::Vector2d(len,len*std::stod(mesh_element->Attribute("ratio"))));
        if(!mesh_type.compare("Rectangle"))
        {
            int nx = std::stoi(mesh_element->Attribute("nx"));
            int ny = std::stoi(mesh_element->Attribute("ny"));
            // double ratio = std::stod(mesh_element->Attribute("ratio"));

            DiamondMesh* dm = new DiamondMesh(1.0,(double)ny/(double)nx,nx,ny,T);
            int i_ori = dm->GetEndingPointIndex();
            int i_ins = dm->GetStartingPointIndex();
            // int i_ori = (ny+1)*(nx)+ny/2;

            // int i_ins = ny/2;
            ms->AddMuscle(p_ori,i_ori,p_ins,i_ins,unit_dir,dm);
        }
        else
        {
            int i_ori,i_ins;

            i_ori = std::stoi(unit->FirstChildElement("mesh")->Attribute("origin_index"));
            i_ins = std::stoi(unit->FirstChildElement("mesh")->Attribute("insertion_index"));
            ms->AddMuscle(p_ori,i_ori,p_ins,i_ins,unit_dir,
                new OBJLoader(unit->FirstChildElement("mesh")->Attribute("path"),T));    
        }
    }
}