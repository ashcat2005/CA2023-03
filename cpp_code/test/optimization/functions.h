#include <armadillo>
#include "Random64.h"
using namespace std;

/// @brief Reproduces the behaviour of numpy's when subtracting a row vector from a column vector. Returns the matrix of a-b
/// @param a vector 1
/// @param b vector 2
arma::mat subtract_vectors(arma::colvec a, arma::rowvec b)
	{	
		arma::mat c = arma::zeros(a.n_elem,b.n_elem);
		c.each_col() += a;
		c.each_row() -= b;
		return c ;
	}
/// @brief Reproduces the behaviour of numpy's when adding a row vector from a column vector. Returns the matrix of a+b
/// @param a 
/// @param b 
/// @return 
arma::mat sum_vectors(arma::colvec a, arma::rowvec b)
	{	

		arma::mat c = arma::zeros(a.n_elem,b.n_elem);
		c.each_col() += a;
		c.each_row() += b;
		// std::cout<<"y1\n";
		return c ;
	}
arma::cube gradiant(arma::cube dq,double h)
	{
		arma::mat x=dq.slice(0);
		arma::mat y=dq.slice(1);
		arma::mat z=dq.slice(2);

		arma::mat r= arma::sqrt(arma::square(x)+arma::square(y)+arma::square(z));
		arma::mat dw = (-2* arma::exp(-pow(r/h,2)))/(pow(h,5)*pow(M_PI,3/2));

		arma::mat dwx=dw*x;
		arma::mat dwy=dw*y;
		arma::mat dwz=dw*z;

		arma::cube gradiant=arma::join_slices(dwx,dwy);
		gradiant=arma::join_slices(gradiant,dwz);
		return gradiant;

	}




/// @brief Calculates the distances between 2 sets of points
arma::cube PairwiseDistance(arma::mat position_set1,arma::mat position_set2)
	{	

		int M = position_set1.n_rows; 	// number of particles in set1
		int N = position_set2.n_rows;	// number of particles in set2

		arma::mat dx = subtract_vectors(position_set1.col(0),position_set2.col(0).t());
		arma::mat dy = subtract_vectors(position_set1.col(1),position_set2.col(1).t());
		arma::mat dz = subtract_vectors(position_set1.col(2),position_set2.col(2).t());


		arma::cube dq=arma::join_slices(dx,dy);
		dq=arma::join_slices(dq,dz);


		return dq;
	}

arma::mat W(arma::cube dq,double h) // ! No se bien que es esto
	{	
		arma::mat dx=dq.slice(0);
		arma::mat dy=dq.slice(1);
		arma::mat dz=dq.slice(2);
		arma::mat r= arma::sqrt(arma::square(dx)+arma::square(dy)+arma::square(dz));
		arma::mat W=pow(1.f/(h*sqrt(M_PI)),3)*arma::exp(-pow(r/h,2));


		
		return W;
	}

/// @brief Calculates the density of the particles
/// @param sampleling_positions Positions where the density is calculated
/// @param pos_matrix position matrix of the particles
/// @param M0 Single particle mass
/// @param h smoothing length
/// @return rho, M x 1 vector with the density of the particles // ! en el codigo original es un vector de aceleracion
arma::mat Density(arma::mat sampleling_positions,arma::mat pos_matrix, double M0,double h)
	{
		// arma::vec rho=arma::zeros(sampleling_positions.n_rows);
		arma::cube dq=PairwiseDistance(sampleling_positions,pos_matrix);
		// auto rho=M0*W(dq,h);
		arma::mat temp_matrix=M0*W(dq,h);
		arma::mat rho = arma::sum(temp_matrix,1);
		

		return rho;
	}



/// @brief Calculates preassure vector
/// @param rho Density vector
/// @param k equation of state constant
/// @param n polytropic index
arma::mat Pressure(arma::mat rho,double k, double n)
	{
		arma::mat P=k*arma::pow(rho,1+(1/n));
		return P;
	}


arma::mat gravity_force(arma::mat pos, double lambda)
	{
		return -lambda*pos;
	}
arma::mat viscous_force(arma::mat velocity, double nu)
	{
		return -nu*velocity;
	}


/// @brief Calculates the acceleration of the particles
/// @param pos Position matrix
/// @param vel Velocity matrix
/// @param M0 Mass of a single particle
/// @param h smoothing length
/// @param k equation of state constant
/// @param n polytropic index
/// @param lambda external force constant
/// @param nu viscosity constant
/// @return The x,y,z components of the acceleration of the particles
arma:: mat Acceleration( arma::mat pos,arma::mat vel, double M0, double h,double k, double n, double lambda, double nu,Crandom & rand64)
	{
		int N=pos.n_rows;
		arma::mat rho=Density(pos,pos,M0,h);
		arma::mat P=Pressure(rho,k,n);
		arma::cube dq=PairwiseDistance(pos,pos);
		arma::cube gradW=gradiant(dq,h);

		arma::mat term1=M0*sum_vectors(P/arma::pow(rho,2), P.t()/arma::pow(rho.t(),2));		
		
		arma::vec pressure_acc_x= - arma::reshape(arma::sum(term1*gradW.slice(0),1),N,1);
		arma::vec pressure_acc_y= - arma::reshape(arma::sum(term1*gradW.slice(1),1),N,1);
		arma::vec pressure_acc_z= - arma::reshape(arma::sum(term1*gradW.slice(2),1),N,1);
		
		arma::mat acc=arma::join_rows(pressure_acc_x,pressure_acc_y);
		acc=arma::join_horiz(acc,pressure_acc_z);
		acc+=gravity_force(pos,lambda);
		acc+=viscous_force(vel,nu);

		return acc;
	}

/// @brief fills a vector with a value
	/// @param vec 
	/// @param value 
void fill_arma(arma::vec & vec, double value)
	{
		int q=vec.size();
		for(int ii=0;ii<q;ii++)
			{vec(ii)=value;}
	}



/// @brief fills a matrix with random values (0,1)
/// @param mat 
/// @param rand64 random number generator
void fill_arma_random(arma::mat & mat,Crandom rand64)
	{
		int q=mat.n_rows;
		int p=mat.n_cols;
		for(int ii=0;ii<q;ii++)
			{
				for(int jj=0;jj<p;jj++)
					{
						mat(ii,jj)=rand64.r();
					}
			}
	}
void fill_arma_random(arma::vec & vec,Crandom rand64)
	{
		int q=vec.size();
		for(int ii=0;ii<q;ii++)
			{
				vec(ii)=rand64.r();
			}
	}

void displayProgressBar(float progress, int barWidth = 50, std::string line="=")
	{
		std::cout << "[";
		int pos = static_cast<int>(progress * barWidth);
		for (int i = 0; i < barWidth; ++i) 
			{
				if (i < pos)
					{std::cout << line;}
				else if (i == pos)
					{std::cout << ">";}
				else
					{std::cout << " ";}
			}
		std::cout << "] " << static_cast<int>(progress * 100.0)<< "%\r";
		std::cout.flush();
	}

double sum_array(double arr[], int size) {
	double sum = 0.0;
	for (int i = 0; i < size; i++) {
		sum += arr[i];
	}
	return sum;
}



