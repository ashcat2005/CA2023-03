import pygame
import sys
import math

# Constants
G = 6.674 * (10 ** -11)  # Gravitational constant
dt = 1000  # Time step in milliseconds

class Planet:
	def __init__(self, mass, position, velocity):
		self.mass = mass
		self.position = position
		self.velocity = velocity

	def update_position(self, time):
		# Update position based on current velocity
		self.position[0] += self.velocity[0] * time
		self.position[1] += self.velocity[1] * time

	def update_velocity(self, other, time):
		# Update velocity based on gravitational force
		distance = math.sqrt((other.position[0] - self.position[0]) ** 2 +
							(other.position[1] - self.position[1]) ** 2)
		force = G * (self.mass * other.mass) / (distance ** 2)

		# Calculate acceleration components
		acceleration_x = force / self.mass * (other.position[0] - self.position[0]) / distance
		acceleration_y = force / self.mass * (other.position[1] - self.position[1]) / distance

		# Update velocity components
		self.velocity[0] += acceleration_x * time
		self.velocity[1] += acceleration_y * time


# Initialize Pygame
pygame.init()

# Constants for the simulation
width, height = 800, 600
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Planet Orbit Simulation")

# Colors
black = (0, 0, 0)
white = (255, 255, 255)

# Create planets
sun = Planet(1.989 * (10 ** 30), [width // 2, height // 2], [0, 500])
earth = Planet(5.972 * (10 ** 24), [width // 2 + 300, height // 2], [0, 500])

# Main loop
clock = pygame.time.Clock()

while True:
	for event in pygame.event.get():
		if event.type == pygame.QUIT:
			pygame.quit()
			sys.exit()

	# Update planet positions and velocities
	earth.update_velocity(sun, dt / 1000)  # Convert time to seconds
	earth.update_position(dt / 1000)

	# Clear the screen
	screen.fill(black)

	# Draw the sun
	pygame.draw.circle(screen, white, (int(sun.position[0]), int(sun.position[1])), 30)

	# Draw the Earth
	print(int(sun.position[0]),int(sun.position[1]))
	print(int(earth.position[0]),int(earth.position[1]))
	pygame.draw.circle(screen, white, (int(earth.position[0]), int(earth.position[1])), 10)
	# pygame.draw.circle(screen,white,(100,100),10)

	# Update the display
	pygame.display.flip()

	# Cap the frame rate
	clock.tick(60)