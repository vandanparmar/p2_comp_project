import numpy as np
from scipy import integrate
from sim_plotting import *
colours = ["#FF0000","#cc0040","#990080","#6600BF","#3300FF","#FFFFFF","#000000"]


#creates elliptical orbits given a closest approach and a second mass (assuming mass1 = 1)
def elliptical_orbit(mass, q_tot):
	theta = np.pi/3
	#theta = 0
	mu = mass/(mass+1.0)
	q1 = mu*q_tot
	q2 = q1/mass
	r1 = q1/(1+np.cos(theta))
	r2 = -r1/mass
	v1 = -2*mass*np.sqrt(1/((1+mass)*(r1-r2)**2))
	v2 = -v1/mass
	print(r1-r2)
	print(v1,v2)
	mass1 = [r1*np.sin(theta/2),r1*np.cos(theta/2),v1,0.0,1.0]
	mass2 = [r2*np.sin(theta/2),r2*np.cos(theta/2),v2,0.0,mass]
	energy = v1**2/2.0+mass*v2**2/2.0 - 2.0*mass/(r1-r2)**2
	print(energy)
	return [mass1,mass2]

#creates parabolic orbits given a closest approach and a second mass (assuming mass1 = 1)
def parabolic_orbit(mass, q_tot):
	theta = np.pi/1.5
	#theta = 0
	mu = mass/(mass+1.0)
	q1 = 2.0*mu*q_tot
	r1 = -q1/(1+np.cos(theta))
	r2 = -r1/mass
	print(r1-r2)
	v1 = -mass*np.sqrt(2/((1+mass)*(r2-r1)))
	v2 = -v1/mass
	mass1 = [r1*np.cos(theta/2),-r1*np.sin(theta/2),0.0,v1,1.0]
	mass2 = [r2*np.cos(theta/2),-r2*np.sin(theta/2),0.0,v2,mass]
	print('initial energy')
	print(energy(np.append(mass1,mass2)))
	return [mass1,mass2]

#calculates the distance between a given pair of masses
def distance(masses):
	r1 = np.sqrt(masses[0]**2+masses[1]**2)
	r2 = np.sqrt(masses[5]**2+masses[6]**2)
	toReturn = r1+r2
	return toReturn

#calculates the energy of a given pair of masses
def energy(masses):
	v1 = np.sqrt(masses[2]**2+masses[3]**2)
	v2 = np.sqrt(masses[7]**2+masses[8]**2)
	#print(v1,v2)
	r1 = np.sqrt(masses[0]**2+masses[1]**2)
	r2 = np.sqrt(masses[5]**2+masses[6]**2)
	m1 = masses[4]
	m2 = masses[9]
	#print(m1,m2)
	toReturn = m1*v1**2/2.0+m2*v2**2/2.0 - m1*m2/(r1+r2)
	return toReturn

#creates individual ring of particles with specified relative density
def create_ring(radius, relative_density):
	no_of_particles = particle_density*relative_density
	velocity = np.sqrt(1/np.sqrt(radius**2+epsilon**2))
	angles = np.linspace(0,2*np.pi,no_of_particles)
	positions = radius*np.exp(1j*angles)
	velocities = velocity*1j*np.exp(1j*angles)
	coords = np.column_stack((np.real(positions),np.imag(positions),np.real(velocities),np.imag(velocities)))
	return [no_of_particles,coords]

#creates a set of rings given an array of radius,density pairs
def create_ring_set(radius_density_pairs,initial_conditions):
	rings = np.zeros((1,4))
	critical_values = []
	tot = 0
	#only looping over 6 values, so no need for numpy loop
	for [radius,density] in radius_density_pairs:
		ring = create_ring(radius,density)
		rings = np.append(rings,ring[1],axis=0)
		tot += ring[0]
		critical_values.append(tot)
	rings = np.delete(rings,0,0)
	rings[:,0] += initial_conditions[0]
	rings[:,1] += initial_conditions[1]
	rings[:,2] += initial_conditions[2]
	rings[:,3] += initial_conditions[3]
	return [critical_values,rings]

#calculates the differential step for the test masses
def g(masses, rings):
	x = rings[:,0]
	y = rings[:,1]
	toReturn = np.zeros_like(rings)
	for mass in masses:
		delta_x = x-mass[0]
		delta_y = y-mass[1]
		r = np.sqrt(np.power(delta_x,2)+np.power(delta_y,2)+np.power(epsilon,2)) #smoothing applied
		r3 = np.power(r,-3) 
		toReturn[:,2] -= np.multiply(r3,delta_x)*mass[4]
		toReturn[:,3] -= np.multiply(r3,delta_y)*mass[4]
	toReturn[:,0] = rings[:,2]
	toReturn[:,1] = rings[:,3]
	toReturn = np.nan_to_num(toReturn)
	return toReturn

#calculates the differential step for the galactic centers
def g_mass(masses):
	x = masses[:,0]
	y = masses[:,1]
	toReturn = np.zeros_like(masses)
	for mass in masses:
		delta_x = x-mass[0]
		delta_y = y-mass[1]
		r = np.sqrt(np.power(delta_x,2)+np.power(delta_y,2))
		r3 = np.power(r,-3)
		r3 = np.clip(r3,0.0,1.0e4) #to avoid infinities from self interactions
		toReturn[:,2] -= np.multiply(r3,delta_x)*mass[4]
		toReturn[:,3] -= np.multiply(r3,delta_y)*mass[4]
	toReturn[:,0] = masses[:,2]
	toReturn[:,1] = masses[:,3]
	toReturn = np.nan_to_num(toReturn)
	return toReturn

#the differential step for wrapper ODE solver
def ode_step(full_set,t,ring_no):
	rings = full_set[:ring_no]
	masses = full_set[ring_no:]
	rings = np.reshape(rings,(-1,4))
	masses = np.reshape(masses,(-1,5))
	dringsdt = g(masses,rings).flatten()
	dmassesdt = g_mass(masses).flatten()
	toReturn = np.append(dringsdt,dmassesdt)
	return toReturn

#differential step for the fixed ODE solver
def ode_step_rev(t,full_set,ring_no):
	rings = full_set[:ring_no]
	masses = full_set[ring_no:]
	rings = np.reshape(rings,(-1,4))
	masses = np.reshape(masses,(-1,5))
	dringsdt = g(masses,rings).flatten()
	dmassesdt = g_mass(masses).flatten()
	toReturn = np.append(dringsdt,dmassesdt)
	return toReturn

#performs a simulation with given initial conditions, can live plot or plot at a particular time, also save data to file.
def indiv_sim(masses,ring_set,totalTime,noOfSteps,timeToPlot, saveName):
	masses = np.reshape(masses,(-1,5))
	masses_f = masses.flatten()
	rings = ring_set[1]
	vals = ring_set[0]
	rings_f = rings.flatten()
	ring_no = len(rings_f)
	full_set_f = np.append(rings_f,masses_f) #creating joined masses and rings
	full_length = len(full_set_f)

	dt = totalTime/noOfSteps
	sol = [full_set_f.flatten()]


	integrator = integrate.ode(ode_step_rev).set_integrator("dopri5")
	integrator.set_initial_value(full_set_f,0.0)
	integrator.set_f_params((ring_no))
	while integrator.successful() and integrator.t<totalTime:
		sol = np.append(sol,[integrator.integrate(integrator.t+dt)],axis=1)
	sol = np.reshape(sol,(-1,full_length))
	ring_sol = sol[:,:ring_no]
	mass_sol = sol[:,ring_no:]
	# energies = []
	# distances = []
	# for i in range(0,noOfSteps):
	# 	energies.append(energy(mass_sol[i,:]))
	# 	distances.append(distance(mass_sol[i,:]))
	# plt.plot(energies)
	# plt.grid()
	# plt.show()
	t = np.linspace(0,totalTime,len(ring_sol))
	print(ring_sol[:,1])
	radii =[[]]
	for i in range(0,4):
		xs = ring_sol[:,0+4*i]
		ys = ring_sol[:,1+4*i]
		radii = np.sqrt(np.power(xs,2)+np.power(ys,2))-(2+i)
		plt.plot(t,radii,c=colours[i])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel('Time / Scaled Units')
	plt.ylabel('Error of Radius of Particle / Scaled Units')
	plt.title('The evolution of radii in circular orbits over a long period of time')
	plt.grid()
	plt.show()
	# # plt.plot(distances)
	# plt.show()
	if (timeToPlot ==0):
		plot_live_full(vals,ring_sol,mass_sol,dt)
	elif (timeToPlot<totalTime):
		index = round(timeToPlot/dt)
		plot_ring_set(vals,np.reshape(ring_sol[index],(-1,4)),mass_sol,timeToPlot,index)
	if (saveName!= ""):
		np.savetxt(saveName,sol,header=str(totalTime)+'\t' +str(noOfSteps)+'\t'+str(ring_no)+'\n'+str(vals))

epsilon = 0.001 #smoothing length
particle_density = 1
#masses = [[0.0,0.0,0.0,0.0,1.0],[0.0,20.0,0.31,0.0,2.0]]
#masses = [[0.0,0.0,0.0,0.0,1.0],[0.0,10.0,0.447,0.0,2.0]]
#masses = [[0.0,0.0,-0.15,0.0,1],[-30,-30,0.15,0.0,1]]
masses = [[0.0,0.0,0.0,0.0,1.0]]


totalTime = 10000
noOfSteps = 20000
timeToPlot = 100000
#fileName = str(masses)+'t='+str(totalTime)+'.txt'
#fileName = 'equal_mass_test.txt'
#masses = parabolic_orbit(1.0,9.0)
#ring_set = create_ring_set([[2,12],[3,18],[4,24],[5,30],[6,36]],masses[0][:4])
ring_set = create_ring_set([[2,1],[3,1],[4,1],[5,1],[6,1]],masses[0][:4])
indiv_sim(masses,ring_set,totalTime,noOfSteps,timeToPlot,'')


# qs = np.linspace(8.5,13.5,51)
# for q in qs:
# 	print(q)
# 	fileName = 'parabolic_equal_mass_q_'+str(int(q*10))+'.txt'
# 	masses = parabolic_orbit(1.0,q)
# 	ring_set = create_ring_set([[2,12],[3,18],[4,24],[5,30],[6,36]],masses[0][:4])
# 	indiv_sim(masses,ring_set,totalTime,noOfSteps,timeToPlot,fileName)



#print(sol[len(sol)-1])
#sol = integrate.odeint(ode_step,full_set_f,t,args=(ring_no,))