#include <cmath>
#include <iostream>
#include <fstream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>	// To threads in SFML
#include "Random64.h"
#include "Vector.h"
#include <chrono>

double h     = 1;  							// smoothing length
double h2	= h*h;							// h^2
double Cd   = 1.f/(4*M_PI*h2*h);			// Constant for the smoothing kernel

const int N_particles=140;

vector3D<double> down_vec(0,-1,0);

const double Lx = 500;
const double Ly = 500;
const double Lz = 500;
const float half_Lx = Lx/2;
const float half_Ly = Ly/2;
const float half_Lz = Lz/2;



float draw_radius=5;
float scale=1;

float grav_constant=10;

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
						vel.x = rand64.r()-0.5;
						vel.y = rand64.r()-0.5;
						vel.z = rand64.r()-0.5;
						vel.unit();
						vel=vel*max_vel;
					}

				mass=1;
			}

		void show_SFML(sf::RenderWindow & window,double x_screen, double y_screen,float scale, float radius)
			{
				shape.setFillColor(HSV_color(preasure*1000,0.8,0.8));
				// shape.setFillColor(sf::Color::White);
				shape.setPosition((pos.x*scale)+x_screen,-pos.y*scale+y_screen);
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
		
		void time_step(Particle * particles, double dt)
			{
				for(int ii; ii<N_particles; ii++)
					{	
						particles[ii].acc.load(0,0,0);
						particles[ii].acc=down_vec;
						particles[ii].vel+=particles[ii].acc*dt;
						particles[ii].pos+=particles[ii].vel*dt;
						
						bounding_collision(particles[ii]);
					}
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

		float Calculate_density(Particle * particles,vector3D<double> sample_point)
			{	
				float density=0;
				vector3D<double> distance(0,0,0);
				float influence=0;
				for(int ii=0;ii<N_particles;ii++)
					{
						vector3D<double> position=particles[ii].pos;
						distance=(sample_point-position);
						influence=smoothing_kernel(distance);
						density+=particles[ii].mass*influence;
					}

				return density;

			}
	void show(sf::RenderWindow & window)
		{
			window.draw(bounding_box);
		}
};


int main()
	{

		int	screen_x=1600;
		int	screen_y=1200;
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SPH_sim");
		float dt=0.01;
		int t_current=0;
		int t_max=100000;
		Crandom rand64(1);
		Interact interaction(Lx,Ly,screen_x,screen_y,scale);
		Particle particles[N_particles];

		for(int ii=0;ii<N_particles;ii++)
			{
				particles[ii].start(rand64,Lx,false,10);
			}

		while (window.isOpen())
			{	
				//---------------------Conditions and input---------------------
				t_current++;
				if(t_current==t_max)
					{window.close();}


				sf::Event event;
				while (window.pollEvent(event))
				{
					if (event.type == sf::Event::Closed)
						window.close();
					
				}

				// ---------------------Simulation---------------------
				interaction.time_step(particles,dt);



				// ---------------------Drawing---------------------
				window.clear(sf::Color::Black);
				interaction.show(window);
				for(int ii=0;ii<N_particles;ii++)
					{
						particles[ii].show_SFML(window,screen_x*0.5,screen_y*0.5,scale,draw_radius);
						if(ii==0)
							{
								particles[ii].show_pos();
							}
					}

				window.display();

				// SFML sleep
				sf::sleep(sf::milliseconds(1));

			}
			//---------------------End of loop---------------------





		return 0;

	}