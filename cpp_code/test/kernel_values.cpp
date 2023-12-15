#include <iostream>
#include <cmath>
#include "Vector.h"
#include <fstream>

double h=1;
double h_1	= 1/h;							// 1/h
double h8	= pow(h,8);						// h^8
double h2	= h*h;							// h^2
double h2_1	= 1/h2;							// 1/h^2
double h3_1	= 1/(h*h*h);						// 1/h^3
double h6	= h2*h2*h2;						// h^6
double Cd   = 8/(M_PI*h2*h);				// Constant for the smoothing kernel
double spiky_constant=6/(M_PI*pow(h,4));			// Constant for the spiky kernel
double smoothing_kernel(vector2D<double> vec)
	{
		double r2=vec.norm2();
		double r=sqrt(r2);
		double value=0;
		double value2=0;
		double h2_r2=h2-r2;

		value=h2_r2*h2_r2*h2_r2;
		if(r<h)
			{return 4*value/(M_PI*h8);}
		else
			{return 0;}
	}
vector2D<double> smoothing_kernel_gradient(vector2D<double> vec)
	{
		double r2=vec.norm2();
		double r=sqrt(r2);
		double value=0;
		double h2_r2=h2-r2;

		value=h2_r2*h2_r2*r;
		if(r<h)
			{return -24*value/(M_PI*h8)*unit(vec);}
		else
			{return 0*unit(vec);}
	}
double spiky_kernel(vector2D<double> vec)
	{
		double value=0;
		double r=vec.norm();
		if(r<h)
			{value=(h-r)*(h-r);}
		else 
			{return 0;}
		

		return value/spiky_constant;
	}
vector2D<double> spiky_gradient(vector2D<double> vec)
	{
		double value=0;
		double r=vec.norm();

		if(r<h)
			{value=(h-r);}
		else 
			{return 0*unit(vec);}

		return -2*value*spiky_constant*unit(vec);
	}

std::ofstream file;
int main()
	{
		file.open("kernel_values.csv");

		for(double ii=0;ii<2;ii+=0.01)
			{
				vector2D<double> vec(ii,0);
				file<<ii<<","<<smoothing_kernel(vec)<<","<<spiky_kernel(vec)<<","
				<<smoothing_kernel_gradient(vec).x<<","<<smoothing_kernel_gradient(vec).y<<","
				<<spiky_gradient(vec).x<<","<<spiky_gradient(vec).y<<std::endl;
				
			}
		file.close();


	}