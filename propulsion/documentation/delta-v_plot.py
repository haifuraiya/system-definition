import numpy as np
import matplotlib.pyplot as plt

# constants for calculation
mu = 3.986004418e14
geo_alt = 42371000
earth_radius = 6371000

# mission constants
mass = 10
# burn_period = 365*24*60*60
burn_period = 120*24*60*60

# pre-calculate the velocity in GEO
end_velocity = np.sqrt(mu / geo_alt)

# specify the parameter stepping for the calculation
number_steps = 40
inclination_range = np.linspace(0, 3*np.pi/5, number_steps)
radius_range = np.linspace(6771000, 10000000, number_steps)
delta_v = np.zeros([len(inclination_range), len(radius_range)])
required_thrust = np.zeros([len(inclination_range), len(radius_range)])

# loop through each radius
for r_x, r in enumerate(radius_range):

	# calculate the start velocity for the orbit
	start_velocity = np.sqrt(mu / r)

	# loop through each inclincation
	for i_x, i in enumerate(inclination_range):

		# calculate the total change in velocity required
		delta_v[r_x][i_x] = np.sqrt( start_velocity**2 + end_velocity**2 - 2*start_velocity*end_velocity*np.cos(i) )

		# calculate the minimum required thrust
		required_thrust[r_x][i_x] = delta_v[r_x][i_x] * mass / burn_period


# create a contour map of the delta V required
CS = plt.contour(180*inclination_range/np.pi, [(_ - earth_radius)/1000 for _ in radius_range], delta_v, cmap="viridis")
plt.clabel(CS, inline=1, fontsize=10, fmt='%1.0f')
plt.suptitle('Required $ \Delta v $ to raise to GEO as a function of starting orbit', fontsize=12, y=0.97)
plt.title("$ \Delta v $ in units of $m \cdot s^{-1}$, assumes circular start and finish orbits", fontsize=8)
plt.xlabel("Inclination ($^{\circ}$)")
plt.ylabel("Altitude (km)")
plt.savefig('deltav_contour.png')
plt.show()

# plot the delta v required at no inclincation
index = 0
plt.plot([(_ - earth_radius)/1000 for _ in radius_range], delta_v[:,index])
plt.suptitle('Required $ \Delta v $ to raise to GEO as a function of starting orbit inclincation', fontsize=12, y=0.97)
plt.title("assumes circular start and finish orbits with inclincation of $%1.1d^{\circ}$" % (180*inclination_range[index]/np.pi), fontsize=8)
plt.xlabel("Altitude (km)")
plt.ylabel("$ \Delta v $ ($m \cdot s^{-1}$)")
plt.grid()
plt.savefig('deltav_incl0.png')
plt.show()

# plot the delta v required at no inclincation
index = 19
plt.plot([(_ - earth_radius)/1000 for _ in radius_range], delta_v[:,index])
plt.suptitle('Required $ \Delta v $ to raise to GEO as a function of starting orbit inclincation', fontsize=12, y=0.97)
plt.title("assumes circular start and finish orbits with inclincation of $%1.1d^{\circ}$" % (180*inclination_range[index]/np.pi), fontsize=8)
plt.xlabel("Altitude (km)")
plt.ylabel("$ \Delta v $ ($m \cdot s^{-1}$)")
plt.grid()
plt.savefig('deltav_incl52.png')
plt.show()



# create a contour map of the delta V required
CS = plt.contour(180*inclination_range/np.pi, [(_ - earth_radius)/1000 for _ in radius_range], 1000*required_thrust, cmap="viridis")
plt.clabel(CS, inline=1, fontsize=10, fmt='%1.0f')
plt.suptitle('Required thrust to raise to GEO as a function of starting orbit', fontsize=12, y=0.97)
plt.title("Thrust in units of $mN$, assumes circular start and finish orbits, orbit elevation takes one year", fontsize=8)
plt.xlabel("Inclination ($^{\circ}$)")
plt.ylabel("Altitude (km)")
plt.savefig('thrust_contour.png')
plt.show()

# plot the delta v required at no inclincation
index = 0
plt.plot([(_ - earth_radius)/1000 for _ in radius_range], 1000*required_thrust[:,index])
plt.suptitle('Required $ \Delta v $ to raise to GEO as a function of starting orbit inclincation', fontsize=12, y=0.97)
plt.title("assumes circular start and finish orbits with inclincation of $%1.1d^{\circ}$" % (180*inclination_range[index]/np.pi), fontsize=8)
plt.xlabel("Altitude (km)")
plt.ylabel("Thrust ($mN$)")
plt.grid()
plt.savefig('thrust_incl0.png')
plt.show()

# plot the delta v required at no inclincation
index = 19
plt.plot([(_ - earth_radius)/1000 for _ in radius_range], 1000*required_thrust[:,index])
plt.suptitle('Required $ \Delta v $ to raise to GEO as a function of starting orbit inclincation', fontsize=12, y=0.97)
plt.title("assumes circular start and finish orbits with inclincation of $%1.1d^{\circ}$" % (180*inclination_range[index]/np.pi), fontsize=8)
plt.xlabel("Altitude (km)")
plt.ylabel("Thrust ($mN$)")
plt.grid()
plt.savefig('thrust_incl52.png')
plt.show()