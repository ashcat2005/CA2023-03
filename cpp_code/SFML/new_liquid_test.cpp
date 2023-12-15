#include <cmath>
#include <iostream>
#include <fstream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>	// To threads in SFML
#include "Random64.h"
#include "Vector.h"
#include <chrono>

double h     = 0.6;  							// smoothing length
double h8	= pow(h,8);						// h^8
double h_1	= 1/h;							// 1/h
double h2	= h*h;							// h^2
double h2_1	= 1/h2;							// 1/h^2
double h3_1	= 1/(h*h*h);						// 1/h^3
double h6	= h2*h2*h2;						// h^6
double Cd   = 8/(M_PI*h2*h);				// Constant for the smoothing kernel
double spiky_constant=15/(M_PI*h6);			// Constant for the spiky kernel

vector2D<double> down(0,-1);
float gravity=10;
const int N_particles=500;


int mid_point=N_particles/2+std::sqrt(N_particles)/4;

int box_size=8;
const double Lx = box_size;
const double Ly = box_size;
const double Lz = box_size;
const float half_Lx = Lx/2;
const float half_Ly = Ly/2;
const float half_Lz = Lz/2;

bool print_cd=false;

float draw_radius=4;
float scale=100;

float target_density=1;
float preasure_multiplier=0.009;

const vector2D<double> right(1,0);

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


sf::ConvexShape arrow()
	{
		sf::ConvexShape arrow;
		arrow.setPointCount(7);
		arrow.setPoint(0, sf::Vector2f(10, 0));
		arrow.setPoint(1, sf::Vector2f(20, 10));
		arrow.setPoint(2, sf::Vector2f(15, 10));
		arrow.setPoint(3, sf::Vector2f(15, 20));
		arrow.setPoint(4, sf::Vector2f(5, 20));
		arrow.setPoint(5, sf::Vector2f(5, 10));
		arrow.setPoint(6, sf::Vector2f(0, 10));
		arrow.setPosition(50, 50);
		arrow.setFillColor(sf::Color::Red);
		arrow.setRotation(90);

		return arrow;
	}
float radians_to_degrees(float radians)
	{
		return radians*180/M_PI;
	}

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
				Update_properties(particles);

				h_1	= 1/h;							// 1/h
				h8	= pow(h,8);						// h^8
				h2	= h*h;							// h^2
				h2_1	= 1/h2;							// 1/h^2
				h3_1	= 1/(h*h*h);						// 1/h^3
				h6	= h2*h2*h2;						// h^6
				Cd   = 8/(M_PI*h2*h);				// Constant for the smoothing kernel
				spiky_constant=6/(M_PI*pow(h,4));			// Constant for the spiky kernel

				for(int ii=0;ii<N_particles;ii++)
					{	
						particles[ii].acc.load(0,0);
						particles[ii].acc+= - Calculate_preasure_gradient(particles,ii)/particles[ii].density;
						particles[ii].acc+= gravity*down;
						
						particles[ii].vel+=particles[ii].acc*dt;
						particles[ii].pos+=particles[ii].vel*dt;
						
						bounding_collision(particles[ii]);
					}	


			}
		void show(sf::RenderWindow & window)
			{window.draw(bounding_box);}
		double smoothing_kernel(vector2D<double> vec)
			{
				double r2=vec.norm2();
				double r=sqrt(r2);
				double value=0;
				double value2=0;
				double h2_r2=h2-r2;

				value=h2_r2*h2_r2*h2_r2;
				if(r<h)
					{return 4*value/(M_PI*h8);}
				else
					{return 0;}
			}
		vector2D<double> smoothing_kernel_gradient(vector2D<double> vec)
			{
				double r2=vec.norm2();
				double r=sqrt(r2);
				double value=0;
				double h2_r2=h2-r2;

				value=h2_r2*h2_r2*r;
				if(r<h)
					{return -24*value/(M_PI*h8)*unit(vec);}
				else
					{return 0*unit(vec);}
			}
		double spiky_kernel(vector2D<double> vec)
			{
				double value=0;
				double r=vec.norm();
				if(r<h)
					{value=(h-r)*(h-r);}
				else 
					{return 0;}
				

				return value/spiky_constant;
			}

		vector2D<double> spiky_gradient(vector2D<double> vec)
			{
				double value=0;
				double r=vec.norm();
		

				if(r<h)
					{value=(h-r);}
				else 
					{return 0*unit(vec);}

				return -2*value*spiky_constant*unit(vec);
			}

		double Calculate_density(Particle * particles, vector2D<double> pos)
			{
				double density=0;
				vector2D<double> distance(0,0);
				for(int ii=0; ii<N_particles; ii++)
					{density+=particles[ii].mass*spiky_kernel(particles[ii].pos-pos);}

				return density;
			}
		double Calculate_property(Particle * particles, vector2D<double> pos)
			{
				float property=0;

				for(int ii=0;ii<N_particles;ii++)
					{	
						property+=particles[ii].property*spiky_kernel(particles[ii].pos-pos)*
						particles[ii].mass/particles[ii].density;
					}
				return property;
			}
		vector2D<double> Calculate_preasure_gradient(Particle * particles, int particle_index)
			{
				vector2D<double> preasure_gradient(0,0);
				vector2D<double> pos=particles[particle_index].pos;
				for(int ii=0;ii<N_particles;ii++)
					{	
						if(ii==particle_index)
							{continue;}

						preasure_gradient+= - particles[ii].preasure*spiky_gradient(particles[ii].pos-pos)*
						particles[ii].mass/particles[ii].density;
					}

				return preasure_gradient;
			}
		double Convert_to_preasure(double density)//. place holder
			{
				return preasure_multiplier*(density-target_density);
			}
		void Update_properties(Particle * particles)
			{	
				bool xd=false;
				float density=0;
				for(int ii=0;ii<N_particles;ii++)
					{	
						
						Densities[ii]=Calculate_density(particles,particles[ii].pos);	
						particles[ii].density=Densities[ii];
						Preasures[ii]=Convert_to_preasure(Densities[ii]);
						particles[ii].preasure=Preasures[ii];
						
					}			
			}
};
int main()
	{

		int	screen_x=1600;
		int	screen_y=1200;
		sf::ConvexShape arrow1=arrow();
		sf::CircleShape influence_circle;
		sf::Color color_influence(0,255,255,90);
		sf::Color color_influence_2(100,100,255,255);
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SPH_sim");
		float dt=0.01;
		int t_current=0;
		int t_max=100000;
		bool grid=false;
		bool show_influence=false;
		float grid_closeness=0.1;
		Crandom rand64(1);
		Interact interaction(Lx,Ly,screen_x,screen_y,scale);
		Particle particles[N_particles];

		for(int ii=0;ii<N_particles;ii++)
			{
				particles[ii].start(rand64,Lx*0.5,false,10);
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
						{h+=0.1;
						std::cout<<"h: "<<h<<"\tdensity: "<<particles[mid_point].density<<"\n";}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Z)) 
						{h-=0.1;if(h<=0){h=0.1;}std::cout<<"h: "<<h<<"\tdensity: "<<particles[mid_point].density<<"\n";}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) 
						{mid_point+=1;}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) 
						{mid_point-=1;}
					if (event.type == sf::Event::KeyReleased) 
						{if(event.key.code==sf::Keyboard::Key::Q)
							{show_influence = !show_influence;
							std::cout<<"show_influence: "<<show_influence<<"\n";}
							
							
							std::cout<<"h: "<<h<<"\tdensity: "<<particles[mid_point].density<<"\n";
						}
					
					

					

				}

				// ---------------------Simulation---------------------
				interaction.time_step(particles,dt);

				// ---------------------Drawing---------------------
				window.clear(sf::Color::Black);
				interaction.show(window);
				for(int ii=0;ii<N_particles;ii++)
					{particles[ii].show_SFML(window,screen_x*0.5,screen_y*0.5,scale,draw_radius);}

				std::cout<<"mid_density: "<<particles[mid_point].density<<"\t"<<"mid_preasure: "<<particles[mid_point].preasure<<"\n";
				if(show_influence)
					{	



						float angle_a=0;
						try{
							if(particles[mid_point].acc.norm2()>0)
							{angle_a=angle(particles[mid_point].acc,right);}
						}
						catch(...){angle_a=0;}
						
						
						arrow1.setRotation(radians_to_degrees(angle_a));
						arrow1.setPosition((particles[mid_point].pos.x*scale)+screen_x*0.5,-particles[mid_point].pos.y*scale+screen_y*0.5);
						influence_circle.setRadius(h*scale);
						influence_circle.setOrigin(h*scale,h*scale);
						influence_circle.setFillColor(color_influence);
						influence_circle.setOutlineColor(color_influence_2);
						influence_circle.setOutlineThickness(5);
						influence_circle.setPosition((particles[mid_point].pos.x*scale)+screen_x*0.5,-particles[mid_point].pos.y*scale+screen_y*0.5);
						window.draw(influence_circle);
						window.draw(arrow1);		
					}
				

				window.display();

				// SFML sleep
				sf::sleep(sf::milliseconds(50));

			}
			//---------------------End of loop---------------------





		return 0;

	}
