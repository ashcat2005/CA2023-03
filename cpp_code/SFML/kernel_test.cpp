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
double h_1	= 1/h;							// 1/h
double h2	= h*h;							// h^2
double h2_1	= 1/h2;							// 1/h^2
double h3_1	= 1/(h*h*h);						// 1/h^3
double h6	= h2*h2*h2;						// h^6
double Cd   = 8/(M_PI*h2*h);				// Constant for the smoothing kernel
double spiky_constant=15/(M_PI*h6);			// Constant for the spiky kernel


const int N_particles=3000;


int mid_point=N_particles/2+std::sqrt(N_particles)/4;

int box_size=5;
const double Lx = box_size;
const double Ly = box_size;
const double Lz = box_size;
const float half_Lx = Lx/2;
const float half_Ly = Ly/2;
const float half_Lz = Lz/2;

bool print_cd=false;

float draw_radius=4;
float scale=100;

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
		vector2D<double> pos, vel, acc;
		double mass;
		double preasure;
		double density;
		double property=0;
		sf::CircleShape shape;

		void start(Crandom & rand64, float max_radius=100, bool random_vel=false, float max_vel=1, float max_acc=1)
			{
				pos.x = (rand64.r()-0.5)*max_radius; 	// max_radius is the spawing radius
				pos.y = (rand64.r()-0.5)*max_radius;
				// pos.y = 0;
				// pos.x = 0;

				shape.setRadius(draw_radius);
				shape.setOrigin(draw_radius,draw_radius);

				vel.load(0,0);
				acc.load(0,0);

				if(random_vel)
					{
						vel.x = rand64.r()-0.5;
						vel.y = rand64.r()-0.5;

						vel.unit();
						vel=vel*max_vel;
					}

				mass=1;
		}
		
		void show_SFML(sf::RenderWindow & window,double x_screen, double y_screen,float scale, float radius)
			{
				// shape.setFillColor(HSV_color(preasure*360/10,0.8,0.8));
				shape.setFillColor(sf::Color::White);
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
			}
		void time_step(Particle * particles, double dt)
			{	
				// Update_properties(particles);

				h_1	= 1/h;							// 1/h
				h2	= h*h;							// h^2
				h2_1	= 1/h2;							// 1/h^2
				h3_1	= 1/(h*h*h);						// 1/h^3
				h6	= h2*h2*h2;						// h^6
				Cd   = 8/(M_PI*h2*h);				// Constant for the smoothing kernel
				spiky_constant=15/(M_PI*h6);			// Constant for the spiky kernel


				// for(int ii; ii<N_particles; ii++)
				// 	{	
				// 		particles[ii].acc.load(0,0,0);
				// 		// particles[ii].acc=down_vec;

				// 		particles[ii].acc+= Calculate_preasure_force(particles,ii)/particles[ii].density;
				// 		if(ii==mid_point)
				// 			{	
				// 				print_cd=true;
				// 				std::cout<<"density: "<<Calculate_density(particles,particles[ii].pos)<<"\n";

				// 				print_cd=false;
				// 			}


				// 		particles[ii].vel=particles[ii].acc*dt;
				// 		bounding_collision(particles[ii]);
				// 	}


			}
		void show(sf::RenderWindow & window)
			{window.draw(bounding_box);}

		double smoothing_kernel(vector2D<double> vec)
			{
				double r=vec.norm();
				double value=0;
				double value2=0;
				double constant=10/(7*M_PI*h*h);
				double r_h=r/h;
				// std::cout<<"r/h: "<<r_h<<"\t";

				// if(0<r_h && r_h<0.5)
				// 	{

				// 		value=1-6*(r_h)*(r_h)+6*(r_h)*(r_h)*(r_h);
				// 		std::cout<<"se activo_1\n";
						
				// 	}
				// else if(r_h>=0.5 && r_h<1)
				// 	{
				// 		value=2*(1-r_h)*(1-r_h)*(1-r_h);
				// 	std::cout<<"se activo_2\n";}
				// else if(r_h>=1)
				// 	{value=0;}

				// else{
				// 		try
				// 			{
				// 				vec.show();
				// 				std::cout<<vec.norm()<<"\n";
				// 				throw std::runtime_error("Error in smothing_kernel: q out of range. r.norm()<0");
				// 			}
				// 		catch(const std::exception& e)
				// 			{std::cerr << e.what() << '\n';}
				// 	}
				/*
				double h2_r2=h*h-r*r;
				value2=h2_r2*h2_r2*h2_r2;
				return value2/(pow(h,8));
				! Esto funciona bien, por alguna razon
				double h2_r2=h*h-r*r;
				value2=h2_r2*h2_r2*h2_r2;
				if(r<h)
					{return value2/(pow(h,8));}
				else
					{return 0;}
				
				*/
				double h2_r2=h*h-r*r;
				value2=h2_r2*h2_r2*h2_r2;
				if(r<h)
					{return value2/(pow(h,8));}
				else
					{return 0;}
				// if(0<r && r<h)
				// 	{	
				// 		double q_2=(2-r/h);
				// 		double q_1=(1-r/h);
				// 		value=q_2*q_2*q_2 - 4*q_1*q_1*q_1;
				// 	}
				// if(h<r && r<2*h)
				// 	{
				// 		double q_2=(2-r/h);
				// 		value=q_2*q_2*q_2;
				// 	}
				// else
				// 	{value=0;}
				// double term_1=(1/(h*std::sqrt(M_PI)));
				// double term_2=std::exp(-r*r/(h*h));

				// return term_1*term_1* term_2;
				

				
				// return value*(10/(14*M_PI*h*h));

			}
		double Calculate_density(Particle * particles, vector2D<double> pos)
			{
				double density=0;
				vector2D<double> distance(0,0);
				float  influence=0;
				for(int ii=0; ii<N_particles; ii++)
					{	
						distance=(particles[ii].pos-pos);
						// std::cout<<"distance: "<<distance.norm()<<"\t";
						influence=smoothing_kernel(distance);
						density+=particles[ii].mass*influence;
						// std::cout<<"density:" <<density<<"\n";
					}
				return density;
			}
};
int main()
	{

		int	screen_x=1600;
		int	screen_y=1200;
		sf::CircleShape influence_circle;
		sf::Color color_influence(0,255,255,90);
		sf::Color color_influence_2(100,100,255,255);
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SPH_sim");
		float dt=0.01;
		int t_current=0;
		int t_max=100000;
		bool grid=true;
		float grid_closeness=0.1;
		Crandom rand64(1);
		Interact interaction(Lx,Ly,screen_x,screen_y,scale);
		Particle particles[N_particles];

		for(int ii=0;ii<N_particles;ii++)
			{
				particles[ii].start(rand64,Lx,false,10);
			}
		if(grid)
			{	int L =std::ceil(std::sqrt(N_particles));
				float x=0;
				float y=0;
				for(int ii=0;ii<N_particles;ii++)
					{	
						x=ii%L*grid_closeness-L*grid_closeness*0.5;
						y=std::floor((float)ii/L)*grid_closeness-L*grid_closeness*0.5;
						particles[ii].pos.load(x,y);
					}
				std::cout<<L<< "grid\n";
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
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A)) 
						
						{

							h+=0.1;
							std::cout<<"density: " <<interaction.Calculate_density(particles,particles[mid_point].pos)<<"\t h: "<<h<<"\n";
						}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Z)) 
						
						{
							h-=0.1;
							if(h<=0)
								{h=0.1;}

							std::cout<<"density: " <<interaction.Calculate_density(particles,particles[mid_point].pos)<<"\t h: "<<h<<"\n";

						}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) 
						
						{
							mid_point+=1;
							particles[mid_point].pos.show();
						}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) 
						
						{
							mid_point-=1;
							particles[mid_point].pos.show();
						}

				}

				// ---------------------Simulation---------------------
				interaction.time_step(particles,dt);

				// ---------------------Drawing---------------------
				window.clear(sf::Color::Black);
				interaction.show(window);
				for(int ii=0;ii<N_particles;ii++)
					{
						particles[ii].show_SFML(window,screen_x*0.5,screen_y*0.5,scale,draw_radius);

					}


				influence_circle.setRadius(h*scale);
				influence_circle.setOrigin(h*scale,h*scale);

				influence_circle.setFillColor(color_influence);
				influence_circle.setOutlineColor(color_influence_2);
				influence_circle.setOutlineThickness(5);

				influence_circle.setPosition((particles[mid_point].pos.x*scale)+screen_x*0.5,-particles[mid_point].pos.y*scale+screen_y*0.5);
				window.draw(influence_circle);

				window.display();

				// SFML sleep
				sf::sleep(sf::milliseconds(1));

			}
			//---------------------End of loop---------------------





		return 0;

	}