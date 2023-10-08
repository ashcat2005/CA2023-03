#include <cmath>
#include <iostream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>	// To threads in SFML
#include "Random64.h"
#include "Vector.h"
using namespace std;

double G=50;


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

		void show_pos(void)
			{
				pos.show();
			}	
		void start(Crandom & rand64, float max_radius=100, bool random_vel=false,bool random_acc=false)
			{
				pos.x = (rand64.r()-1)*max_radius;
				pos.y = (rand64.r()-1)*max_radius;
				pos.z = (rand64.r()-1)*max_radius;
				pos.show();
				vel.load(0,0,0);
				acc.load(0,0,0);

				if(random_vel)
					{
						vel.x = rand64.r();
						vel.y = rand64.r();
						vel.z = rand64.r();
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

				// pos.show_2();
				// vel.show_2();
				// f.show();
				// std::cout<<"------\n";
				acc=f/mass;
				vel+=acc*dt;
				pos+=vel*dt;
			}
		void show_SFML(sf::RenderWindow & window,double x_screen, double y_screen,float scale)
			{
				sf::CircleShape shape(10.f*scale);
				// vel.show();
				shape.setFillColor(HSV_color(vel.norm()*10,0.8,0.8));
				shape.setPosition((pos.x*scale)+x_screen,pos.y*scale+y_screen);
				window.draw(shape);
			}



};



int main()
	{	
		int N_particles=1000;
		double t_max=1000;
		double dt=0.1;
		double N_steps=t_max/dt;
		int	screen_x=1600;
		int	screen_y=1200;
		int t_current=0;
		sf::RenderWindow window(sf::VideoMode(screen_x, screen_y), "SFML Window");

		Particle particle[N_particles];
		Crandom rand64(1);

		for(int ii=0; ii<N_particles; ii++)
			{
				particle[ii].start(rand64,500.0F,true);
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



			for(int ii=0; ii<N_particles; ii++)
				{	
					// particle.show_pos();
					particle[ii].move(dt);
					particle[ii].show_SFML(window,screen_x/2,screen_y/2,0.5);
				}

			// window.draw(...);


			window.display();
		}
	}