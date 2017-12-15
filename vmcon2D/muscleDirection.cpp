#include <iostream>
#include <cmath>
#include <cstdlib>
#include "muscleDirection.h"

void get_muscle_direction(Eigen::MatrixXd& v, int nx, int ny, int nn)
{

	//std::cout << "hihihhihihihihhi " <<std::endl;
	v.resize(ny, nx);
	v.setZero();

	Eigen::MatrixXd isBoundary(ny, nx);
	isBoundary.setZero();

	//boundary condition

	for (int xidx = 0 ; xidx < nx ; xidx++)
	{
		v(0,xidx) = 10;
		isBoundary(0,xidx) = 1;
		v(ny-1,xidx) = -10;
		isBoundary(ny-1,xidx) = 1;
	}
	

	for (int i= 0; i < nn ; i++)
	{
		v(nn-1-i,i) = 10.;
		v(nn-1-i,nx-i-1) = 10.;
		isBoundary(nn-i-1,i) = 1;
		isBoundary(nn-i-1,nx-i-1) = 1;
	} 

	for (int i= nn; i < ny-nn ; i++)
		v(i,0) = -20. * i / (ny+1-2*nn) + 10. *(ny-1)/(ny+1-2*nn);

	for(int yidx = 1; yidx < ny-1 ; yidx++)
		isBoundary(yidx,0) = 1;


	for (int i= 0; i < nn ; i++)
	{
		v(ny-nn+i,i) = -10.;
		v(ny-nn+i,nx-i-1) = -10.;
		isBoundary(ny-nn+i,i) = 1;
		isBoundary(ny-nn+i,nx-i-1) = 1;
	} 

	for (int i= nn; i < ny-nn ; i++)
		v(i,nx-1) = -20. * i / (ny+1-2*nn) + 10. *(ny-1)/(ny+1-2*nn);

	for(int yidx = 1; yidx < ny-1 ; yidx++)
		isBoundary(yidx,nx-1) = 1;


	//std::cout<<"isBoundary : \n"<<isBoundary<<std::endl;

	//std::cout<<"V : \n"<<v<<std::endl;

	// solve Poission's equation
	double eps = 2.0e-2;
	double gap = 100.;
	double temp = 0;
	int firstFlag = 1;
	int num_iter = 0;
	//while (gap>eps)
	while(num_iter < 100)
	{
		//std::cout << "1111gap: " << gap <<std::endl;
		for(int yidx = 1; yidx < ny-1 ; yidx++)
			for (int xidx = 1 ; xidx < nx-1 ; xidx++)
				if (isBoundary(yidx, xidx) < 1.)
					v(yidx,xidx) = .25*(v(yidx+1, xidx) + v(yidx-1, xidx) + v(yidx, xidx+1) + v(yidx, xidx-1));
		if (firstFlag == 1)
		{
			gap = 0;
			for(int yidx = 1; yidx < ny-1 ; yidx++)
				for (int xidx = 1 ; xidx < nx-1 ; xidx++)
					gap = gap + v(yidx,xidx);
			gap = gap /((ny-2)*(nx-2));
			firstFlag = 0;
			//std::cout << "first gap: " << gap <<std::endl;
		}
		else
		{
			temp = 0;
			for(int yidx = 1; yidx < ny-1 ; yidx++)
				for (int xidx = 1 ; xidx < nx-1 ; xidx++)
					temp = temp + v(yidx,xidx);
			temp = temp/((ny-2)*(nx-2));
			//std::cout << "!!!!!gap: " << gap <<std::endl;
			//std::cout << "!!!!!temp: " << temp <<std::endl;
			gap = std::abs(gap - temp);
			//std::cout << "result gap: " << gap <<std::endl;
		}
		num_iter++;
	}
	//std::cout << "iter_num: " << num_iter <<std::endl;
}