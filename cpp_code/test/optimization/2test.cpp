#include <iostream>
#include <chrono>
#include <armadillo>
#include "functions.h"
using namespace std;


int main()
	{


		arma::mat a = arma::randu<arma::mat>(1000,1000);

		int N = 1000;
		double duration[N];
		double rand_int=0;
		std::chrono::time_point <std::chrono::high_resolution_clock> t1,t2;
		for(int ii=0; ii<N;ii++)
			{
				t1 = std::chrono::high_resolution_clock::now();
				arma::mat P = Pressure(a,0.1,1);
				t2 = std::chrono::high_resolution_clock::now();
				duration[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
				displayProgressBar(float(ii)/float(N),50);
			}

		std::cout<<"Pressure: "<<sum_array(duration,N)/N<<" ns"<<std::endl;
		//  Pressure: 2.12879e+07 ns
		// 	Pressure: 3.59613e+06 ns -O3
	}

