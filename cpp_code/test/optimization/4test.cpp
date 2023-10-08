#include <iostream>
// #include <armadillo>
#include "functions.h"
// arma::mat subtract_vectors(arma::colvec a, arma::rowvec b)
// 	{	
// 		arma::mat c = arma::zeros(a.n_elem,b.n_elem);
// 		c.each_col() += a;
// 		c.each_row() -= b;
// 		return c ;
// 	}


arma::mat subtract_vectors_2(arma::colvec a, arma::rowvec b)
	{	
		return arma::repmat(a, 1, b.n_elem) - arma::repmat(b, a.n_elem, 1);
	}

arma::mat subtract_vectors_3(arma::colvec a, arma::rowvec b)
	{	
		arma::mat c = arma::zeros(a.n_elem,b.n_elem);

		for(int i=0;i<a.n_elem;i++)
			{
				for(int j=0;j<b.n_elem;j++)
					{
						c(i,j)=a(i)-b(j);
					}
			}
		return c;
	}

int main()
{
    arma::vec a = arma::randu<arma::vec>(1000);
    arma::vec b = arma::randu<arma::vec>(1000);

	int N = 5000;
	double duration1[N];
	double duration2[N];
	double duration3[N];


	std::chrono::time_point <std::chrono::high_resolution_clock> 
	ta,tb,tc,td,te,tf;
	for(int ii=0; ii<N;ii++)
		{
			ta = std::chrono::high_resolution_clock::now();
			arma::mat c = subtract_vectors(a,b.t());
			tb = std::chrono::high_resolution_clock::now();
			duration1[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( tb - ta ).count();


			tc = std::chrono::high_resolution_clock::now();
			arma::mat d = subtract_vectors_2(a,b.t());
			td = std::chrono::high_resolution_clock::now();
			duration2[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( td - tc ).count();


			te = std::chrono::high_resolution_clock::now();
			arma::mat f = subtract_vectors_3(a,b.t());
			tf = std::chrono::high_resolution_clock::now();
			duration3[ii] = std::chrono::duration_cast<std::chrono::nanoseconds>( tf - te ).count();

			displayProgressBar(float(ii)/float(N),50);
		}
	std::cout<<"\n";

	std::cout	<<"subtract_vectors_1: "<<sum_array(duration1,N)/N<<" ns \n"
				<<"subtract_vectors_2: "<<sum_array(duration2,N)/N<<" ns \n"
				<<"subtract_vectors_3: "<<sum_array(duration3,N)/N<<" ns \n";





    return 0;
}