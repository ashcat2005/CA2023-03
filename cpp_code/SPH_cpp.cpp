#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <armadillo>
#include "functions.h"

using namespace std;
// SPH code for star structure simulation //
// - Simulation Parameters - //

int N0 = 100; 		// ! ??

const int N = 1; 		// Number of particles
double t = 0; 		// current time 
double tMax = 1; 	// maximum time
double dt = 0.04; 	// time step
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


arma::mat pos=arma::ones(N,d); // position matrix
arma::mat vel=arma::zeros(N,d); // velocity matrix. jsdoc




int main()
	{

		Crandom rand64(1);	// random number generator
		arma::mat acc=Acceleration(pos,vel,M0,h,k,n,lambda,nu,rand64);
		// ! No se que sea esto de rlin y rr
		arma::mat rr=arma::zeros(N0,d);
		arma::vec rlin=arma::linspace(0,1,N0); 
		rr.col(0)=rlin;
		arma::vec rho_analytic = lambda*(R*R-arma::square(rlin))/(4*k);

		arma::cube data = arma::zeros(nStep,N,d);
		arma::mat rho_data= arma::zeros(nStep,N);
		arma::mat density_data= arma::zeros(nStep,N0);
		// std::cout<<"x\n";
		for(int ii=0; ii<nStep;ii++)
			{
				acc.print();
				vel += acc*dt/2;
				pos += vel*dt;
				acc=Acceleration(pos,vel,M0,h,k,n,lambda,nu,rand64);
				vel += acc*dt/2;

				t += dt;

				data.row(ii)=pos;			
			}
	}
