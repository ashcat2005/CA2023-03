#include <iostream>
#include <chrono>
#include <armadillo>
#include "functions.h"
using namespace std;


int N0 = 100; 		// ! ??

const int N = 1000; 		// Number of particles
double t = 0; 		// current time 
double tMax = 10; 	// maximum time
double dt = 0.04; 	// time step
double dt_2=dt/2; 	// half time step
int nStep = tMax/dt; // number of steps
double M = 2; 		// star mass
double R = 0.75; 		// star radius
double h = 0.1; 	// smoothing length
double k = 0.1;		// equation of state constant
double n = 1; 		// polytropic index
double nu = 2; 		// viscosity constant
const int d = 3; 		// dimension

double lambda=  (2.*k*(1.+n)* pow(M_PI,-3./(2.*n)) / (R*R)) * pow(M*tgamma(5./2.+n) / (R*R*R*tgamma(1.+n)), 1./n);// external force constant

// double lambda=tgamma(5/2+n) ;
double M0 = M/N; 	// single particle mass

int main()
	{
		arma::mat pos = arma::randu<arma::mat>(1000,3);
		arma::mat vel = arma::randu<arma::mat>(1000,3);
		Crandom rand64(1);	// random number generator
		int N = 1000;
		double duration[N];

		std::chrono::time_point <std::chrono::high_resolution_clock> t1,t2;
		for(int ii=0; ii<N;ii++)
			{
				t1 = std::chrono::high_resolution_clock::now();
				arma::cube q=PairwiseDistance(pos,pos);
				t2 = std::chrono::high_resolution_clock::now();
				duration[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
				displayProgressBar(float(ii)/float(N),50);
			}
		
		std::cout<<"Pairwise "<<sum_array(duration,N)/N<<" ns"<<std::endl;
		// 2.99204e+07 ns
		// 1.80338e+07 ns -O3

	}