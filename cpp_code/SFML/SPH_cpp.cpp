#include <cmath>
#include <iostream>
#include <fstream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>	// To threads in SFML
#include "Random64.h"
#include "Vector.h"
using namespace std;

ofstream salida;

double G=50;

const int N_particles=800;		// Number of particles
double t     = 0;  		// current time of the simulation
double tEnd  = 12; 		// time at which simulation ends
double dt    = 0.04;   	// timestep
double M     = 2;  		// star mass
double R     = 0.75;  	// star radius
double h     = 1;  	// smoothing length
double k     = 0.1;  	// equation of state constant
double n     = 1;  		// polytropic index
double nu    = 1;  		// damping
double lambda=  (2.*k*(1.+n)* pow(M_PI,-3./(2.*n)) / (R*R)) * pow(M*tgamma(5./2.+n) / (R*R*R*tgamma(1.+n)), 1./n);
double h2	= h*h;		// h^2
double W_c	= pow(1.f/(h*sqrt(M_PI)),3);	// Constant for the smoothing kernel
double Cd   = 1.f/(4*M_PI*h2*h);	// Constant for the smoothing kernel

// TODO a velocity map to show the velocity of the particles to capture a range of color




sf::Color HSV_color(float H, float S, float V)
{
	float C = S * V; // Chroma
	float HPrime = std::fmod(H / 60, 6.f); // H'
	float X = C * (1 - std::fabs(std::fmod(HPrime, 2.f) - 1));
	float M = V - C;

	float R = 0.f;
	float G = 0.f;
	float B = 0.f;

	switch (static_cast<int>(HPrime))
	{
	case 0: R = C; G = X;        break; // [0, 1)
	case 1: R = X; G = C;        break; // [1, 2)
	case 2:        G = C; B = X; break; // [2, 3)
	case 3:        G = X; B = C; break; // [3, 4)
	case 4: R = X;        B = C; break; // [4, 5)
	case 5: R = C;        B = X; break; // [5, 6)
	}

	R += M;
	G += M;
	B += M;

	sf::Color color;
	color.r = static_cast<sf::Uint8>(std::round(R * 255));
	color.g = static_cast<sf::Uint8>(std::round(G * 255));
	color.b = static_cast<sf::Uint8>(std::round(B * 255));

	return color;
}


class Particle
{
	public:

		vector3D<double> pos, vel, acc;
		double mass;
		double preasure;
		double density;

		void show_pos(void)
			{
				pos.show();
			}	
		void start(Crandom & rand64, float max_radius=100, bool random_vel=false,bool random_acc=false, float max_vel=1, float max_acc=1)
			{
				pos.x = (rand64.r()-0.5)*max_radius;
				pos.y = (rand64.r()-0.5)*max_radius;
				pos.z = (rand64.r()-0.5)*max_radius;
				// pos.show();
				vel.load(0,0,0);
				acc.load(0,0,0);

				if(random_vel)
					{
						vel.x = rand64.r()*max_vel;
						vel.y = rand64.r()*max_vel;
						vel.z = rand64.r()*max_vel;
					}

				if(random_acc)
					{
						acc.x = rand64.r();
						acc.y = rand64.r();
						acc.z = rand64.r();
					}
				mass=1;
			}


		vector3D<double> force()
			{
				vector3D<double> force(0,0,0);
				if (pos.norm2()<0.1)
					{return force;}
				double r2=pos.norm();
				force=-unit(pos)*G/(1+r2);
				
				return force;
			}

		void move(double dt)
			{	
				vector3D<double> f=force();

				acc=f/mass;
				vel+=acc*dt;
				pos+=vel*dt;
			}
		void show_SFML(sf::RenderWindow & window,double x_screen, double y_screen,float scale, float radius)
			{
				sf::CircleShape shape(radius*scale);
				// vel.show();
				shape.setFillColor(HSV_color(density*50,0.8,0.8));
				// shape.setFillColor(sf::Color::White);
				shape.setPosition((pos.x*scale)+x_screen,pos.y*scale+y_screen);
				window.draw(shape);
			}

		friend class Interact;

};


class Interact
	{
		public:

			double Densities[N_particles];
			double Preasures[N_particles];


			double smoothing_kernel(vector3D<double> r)
				{
					double q=r.norm();
					double val=0;
					double q_2=2-q;
					double q_1=1-q;
					if(0<=q && q<=1)
						{val=q_2*q_2*q_2  -  4*(q_1*q_1*q_1);}
					else if(1<=q && q<=2)
						{val=q_2*q_2*q_2;}

					else if(q>2)
						{val=0;}

					else{
						try{
							r.show();
							std::cout<<r.norm()<<"\n";
							throw std::runtime_error("Error in smothing_kernel2: q out of range. r.norm()<0");
						}
						catch(const std::exception& e){
							std::cerr << e.what() << '\n';
						}
					}
					return Cd*val;
				}

			vector3D<double> smoothing_derivative(vector3D<double> distance)
				{
					double q=distance.norm();
					double val=0;

					if(0<q && q<1)
						{val=9*(q*q)-12*q;}
					else if(1<q && q<2)
						{val=-3*(2-q)*(2-q);}
					
					else
						{	
							val=0;
							return val*distance;
						}
					
					return Cd*val*unit(distance);
				}

			double Density(Particle * particles ,vector3D<double> position)
				{
					double density=0;
					for(int ii=0; ii<N_particles; ii++)
						{density+=particles[ii].mass*smoothing_kernel(position-particles[ii].pos);}

					
					return density;
				}

			double Preasure(double density)
				{
					return k*pow(density,1.f+1.f/n);
				}

			void Calculate_properties(Particle * particles)
				{
					for(int ii=0; ii<N_particles; ii++)
						{
							Densities[ii]=Density(particles,particles[ii].pos);
							particles[ii].density=Densities[ii];
							Preasures[ii]=Preasure(Densities[ii]);
							particles[ii].preasure=Preasures[ii];
						}
				}

			vector3D<double> gravity_force(vector3D<double> pos, int index)   
				{	
					return -lambda*pos;
				}
			vector3D<double> viscous_force(vector3D<double> vel)   {return -nu*vel;}

			vector3D<double> particle_interactions(Particle * particles, int index)
			{
				vector3D<double> acc(0,0,0);
				vector3D<double> gradient(0,0,0);
				
				for(int jj=0; jj<N_particles; jj++)
					{	
						if(jj!=index){
						// std::cout<<"jj="<<jj<<"\t"<<"index="<<index<<"\n";
						// acc.show();
						gradient=smoothing_derivative(particles[index].pos-particles[jj].pos);
						// gradient.show();
						acc-=particles[jj].mass*Preasures[jj]/Densities[jj]*gradient;
						// acc.show();
						acc-=Preasures[index]/Densities[index]*gradient;
						// std::cout<<"Preasures["<<index<<"] = "<<Preasures[index]<<"\t" <<"Densities["<<index<<"] = "<<Densities[index]<<"\n";

						// acc.show();
						// if(Preasures[jj]/Densities[jj])
						// 	{
						// 		// std::cout<<"Particle interactions:\n";
						// 		std::cout<<"Preasures[jj]/Densities[jj]="<<Preasures[jj]/Densities[jj]<<"\n";
						// 	}
						}
						
					}
				
				if (gradient.norm()>1e5)
					{
						std::cout << "Particle interactions:\n";
						std::cout << "gradient.norm()=" << gradient.norm() << "\n";
					}
				if(acc.norm()>1e5)
					{	
						// std::cout<<"Particle interactions:\n";
						std::cout<<"acc.norm()="<<acc.norm()<<"\n";
					}
				return acc;
			}

			void time_step(Particle * particles, double dt)
				{	
					Calculate_properties(particles);
					vector3D<double> pos_prev(0,0,0);
					vector3D<double> vel_prev(0,0,0);
					vector3D<double> acc_prev(0,0,0);
					vector3D<double> grav_prev(0,0,0);
					for(int ii=0;ii<N_particles;ii++)
						{	
							// std::cout<<"ii="<<ii<<"\t";

							// if (ii==33)
							// 	{
							// 		// particles[ii].pos.show_2();std::cout<<"\t";
							// 		// particles[ii].vel.show_2();std::cout<<"\t";
							// 		// particles[ii].acc.show_2();std::cout<<"\n";
							// 		pos_prev=particles[ii].pos;
							// 		vel_prev=particles[ii].vel;
							// 		acc_prev=particles[ii].acc;

							// 	}
							particles[ii].vel+=particles[ii].acc*dt/2;	// kick
							particles[ii].pos+=particles[ii].vel*dt;	// drift
							// particles[ii].acc.load(0,0,0);

							
							particles[ii].acc=gravity_force(particles[ii].pos, ii);
							particles[ii].acc+=viscous_force(particles[ii].vel);
							particles[ii].acc+=particle_interactions(particles,ii);

							if(ii==33)
								{
									std::cout<<particles[ii].preasure<<"\t";
									std::cout<<particles[ii].density<<"\n";

								}

							particles[ii].vel+=particles[ii].acc*dt/2;	// kick

						}
					
				}
			
			
	};


int main()
	{	

		salida.open("data.csv");

		
		double t_max=1000;
		double dt=0.001;
		double N_steps=t_max/dt;
		int	screen_x=1600;
		int	screen_y=1200;
		int t_current=0;
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SFML Window");

		Interact interaction;
		Particle particle[N_particles];
		Crandom rand64(3);

		for(int ii=0; ii<N_particles; ii++)
			{
				particle[ii].start(rand64,50.0F,false,false,100);
			}




		// Run the program as long as the window is open
		while (window.isOpen())
		{
			t_current++;
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();// "close requested" event: we close the window
				
			}
			if (t_current>N_steps)
					{window.close();}

			// Clear the window with a black color
			window.clear(sf::Color::Black);



			interaction.time_step(particle,dt);
			for(int ii=0; ii<N_particles; ii++)
				{
					particle[ii].show_SFML(window,screen_x/2,screen_y/2,50.f,0.1f );
					// salida<<t_current<<","<<particle[0].pos.x<<"\t"<<particle[0].vel.x<<"\t"<<particle[0].acc.x<<"\n";
				}

			// window.draw(...);


			window.display();
			sf::sleep(sf::milliseconds(5));
		}
	salida.close();
	}