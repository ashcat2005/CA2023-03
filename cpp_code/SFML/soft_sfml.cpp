#include <cmath>
#include <iostream>
#include <fstream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>	// To threads in SFML
#include "Random64.h"
#include "Vector.h"
#include <chrono>

ofstream times;

// Define bounds
const double Lx = 500;
const double Ly = 500;
const double Lz = 500;
const float half_Lx = Lx/2;
const float half_Ly = Ly/2;
const float half_Lz = Lz/2;



bool central_force=false;
// Define constants for liquid- this only works for not central force
const float preasure_multiplier=100;
const float target_density=1;
const double G = 50;
float grav_const=10;
// Define constants for drawing
float draw_radius=5;
float scale=1;

// Define constants for simulation
const int N_particles=1400;


double tEnd  = 12; 		// time at which simulation ends
double dt    = 0.04;   	// timestep
double M     = 2;  		// star mass
double R     = 0.75;  	// star radius
double h     = 1;  		// smoothing length
double k     = 0.1;  	// equation of state constant
double n     = 1;  		// polytropic index
double nu    = 1;  		// damping
double lambda=  (2.*k*(1.+n)* pow(M_PI,-3./(2.*n)) / (R*R)) * pow(M*tgamma(5./2.+n) / (R*R*R*tgamma(1.+n)), 1./n);
double h2	= h*h;		// h^2
double W_c	= pow(1.f/(h*sqrt(M_PI)),3);	// Constant for the smoothing kernel
double Cd   = 1.f/(4*M_PI*h2*h);			// Constant for the smoothing kernel


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
// double G=50;

class Particle
{
	public:

		vector3D<double> pos, vel, acc;
		double mass;
		double preasure;
		double density;
		sf::CircleShape shape;

		void show_pos(void){pos.show();}
		
		void start(Crandom & rand64, float max_radius=100, bool random_vel=false, float max_vel=1, float max_acc=1)
			{
				pos.x = (rand64.r()-0.5)*max_radius; 	// max_radius is the spawing radius
				pos.y = (rand64.r()-0.5)*max_radius;
				pos.z = (rand64.r()-0.5)*max_radius;
				pos.z = 0;
				// pos.y = 0;
				// pos.x = 0;

				// shape(radius*scale);
				shape.setRadius(draw_radius);
				shape.setOrigin(draw_radius,draw_radius);
				// pos.show();
				vel.load(0,0,0);
				acc.load(0,0,0);
				// vel.x = 15;
				if(random_vel)
					{
						vel.x = rand64.r()*max_vel;
						vel.y = rand64.r()*max_vel;
						vel.z = rand64.r()*max_vel;
					}

				mass=1;
			}

		void show_SFML(sf::RenderWindow & window,double x_screen, double y_screen,float scale, float radius)
			{
				shape.setFillColor(HSV_color(preasure*1000,0.8,0.8));
				// shape.setFillColor(sf::Color::White);
				shape.setPosition((pos.x*scale)+x_screen,pos.y*scale+y_screen);
				window.draw(shape);
			}

		friend class Interact;

};

class Interact
{
	
	public:

		int X_domain,Y_domain;
		sf::RectangleShape bounding_box;
		double Densities[N_particles];
		double Preasures[N_particles];

		vector3D<double> gravity_force(vector3D<double> pos, int index)   
				{	
					vector3D<double> dir(0,-1,0);
					if(central_force)
						{
							return -lambda*pos;
						}
					else
						{
							return -grav_const*M*dir;
						}	
				}

		Interact(int Lx,int Ly,double x_screen, double y_screen,float scale=1.f)
			{
				X_domain=Lx;
				Y_domain=Ly;
				bounding_box.setSize(sf::Vector2f(X_domain*scale,Y_domain*scale));
				bounding_box.setOutlineThickness(2);
				bounding_box.setOutlineColor(sf::Color::White);
				bounding_box.setFillColor(sf::Color::Transparent);
				bounding_box.setPosition((x_screen-X_domain*scale)*0.5f,(y_screen-Y_domain*scale)*0.5f);

			}

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

		float Density_to_preasure(float density)
			{
					float Density_delta=density-target_density;
					float preasure_delta=Density_delta*preasure_multiplier;
					return preasure_delta;
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

		void bounding_collision(Particle & particle)
			{
				if(particle.pos.x>half_Lx)
					{
						particle.pos.x=half_Lx;
						particle.vel.x*=-1;
					}
				if(particle.pos.x<-half_Lx)
					{
						particle.pos.x=-half_Lx;
						particle.vel.x*=-1;
					}
				if(particle.pos.y>half_Ly)
					{
						particle.pos.y=half_Ly;
						particle.vel.y*=-1;
					}
				if(particle.pos.y<-half_Ly)
					{
						particle.pos.y=-half_Ly;
						particle.vel.y*=-1;
					}
				if(particle.pos.z>half_Lz)
					{
						particle.pos.z=half_Lz;
						particle.vel.z*=-1;
					}
				if(particle.pos.z<-half_Lz)
					{
						particle.pos.z=-half_Lz;
						particle.vel.z*=-1;
					}


			}
		
		float Shared_preasure(float densityA,float denistyB)
			{
				float pA=Density_to_preasure(densityA);
				float pB=Density_to_preasure(denistyB);
				return (pA+pB)*0.5f;
			}

		vector3D<double>  preasure_force(Particle * particles,vector3D<double> position)
			{
				vector3D<double> force(0,0,0);
				float dist=0;
				vector3D<double> slope(0,0,0);
				float density=0;
				vector3D<double> dir(0,0,0);

				for(int ii=0; ii<N_particles; ii++)
					{
						dist=(position-particles[ii].pos).norm();
						if(dist==0)
							{
								continue;
							}
						dir=unit(position-particles[ii].pos);
						slope=smoothing_derivative(position-particles[ii].pos);
						// density=Densities[ii];
						// float shared_preasure=Shared_preasure(Densities[ii]);
						force+= -Density_to_preasure(Densities[ii])*slope*particles[ii].mass/Densities[ii];
					}
				return force;
			}

		void time_step(Particle * particles, double dt)
				{	
					Calculate_properties(particles);
					for(int ii=0;ii<N_particles;ii++)
						{	
							bounding_collision(particles[ii]);
							particles[ii].vel+=particles[ii].acc*dt/2;	// kick
							particles[ii].pos+=particles[ii].vel*dt;	// drift

							particles[ii].acc=gravity_force(particles[ii].pos, ii);
							particles[ii].acc+=preasure_force(particles,particles[ii].pos)/Densities[ii];
							// particles[ii].vel.show();
							particles[ii].vel+=particles[ii].acc*dt/2;	// kick

						}
					
				}


				




	void show(sf::RenderWindow & window)
		{
			window.draw(bounding_box);
		}
};




int main()
	{	

		bool measure_time=false;
		int	screen_x=1600;
		int	screen_y=1200;
		
		if(measure_time)times.open("runtime.txt",std::ios::app);

		int t_current=0;

		int t_max=500;
		double frame_runtime[t_max];
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SFML Window");
		Crandom rand64(1);
		Particle particle[N_particles];

		for(int ii=0;ii<N_particles;ii++)
			{
				particle[ii].start(rand64,100,false);
			}
		

		Interact interaction(Lx,Ly,screen_x,screen_y,scale);
		while (window.isOpen())
		{	
			
			t_current++;
			
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();// "close requested" event: we close the window
				
			}


			if(t_current==t_max)
				{
					window.close();
				}

			// - Update here -
			auto start_time = std::chrono::high_resolution_clock::now();
			interaction.time_step(particle,dt);
			auto end_time = std::chrono::high_resolution_clock::now();
			// - end update -
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

			frame_runtime[t_current]=duration.count();
			window.clear(sf::Color::Black);// Clear the window with a black color

			// - Draw here -

			interaction.show(window); 			// Draw the bounding box
			for(int ii=0;ii<N_particles;ii++)	// Draw the particles
				{
					particle[ii].show_SFML(window,screen_x/2,screen_y/2,scale,10);
				}

			// - Stop drawing here -
			window.display();

			


		}

		double avg=0;
		if(measure_time){
			for(int ii=0;ii<t_max;ii++)
				{avg+=frame_runtime[ii];}
			}


		if(measure_time){times<<N_particles<<","<<avg/t_max<<"\n";}
		times.close();
	}