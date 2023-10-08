#include <iostream>
#include <chrono>
#include <armadillo>
#include "functions.h"
using namespace std;


int main()
	{
		// write a duration test for the functions in functions.h
		Crandom rand64(1);	// random number generator
		arma::mat a = arma::randu<arma::mat>(1000,1000);
		arma::mat b = arma::randu<arma::mat>(1000,1000);
		int N = 10000;
		double duration[N];
		double rand_int=0;
		std::chrono::time_point <std::chrono::high_resolution_clock> t1,t2;
		for(int ii=0; ii<N;ii++)
			{	rand_int=rand64.int32()%1000;
				// std::cout<<rand_int<<"\n";
				t1 = std::chrono::high_resolution_clock::now();
				arma::mat c = sum_vectors(a.col(rand_int),b.row(rand_int));
				t2 = std::chrono::high_resolution_clock::now();
				duration[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
				displayProgressBar(float(ii)/float(N),50);
			}
		std::cout<<"sum_vectors: "<<sum_array(duration,N)/N<<" ns"<<std::endl;
		//    	sum_vectors: 4.3945e+06 ns
		// 		sum_vectors: 2.62085e+06 ns -O3

	}