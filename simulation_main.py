import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


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

#plots a ring set with individual colours for each ring
def plot_ring_set(critical_values, rings, mass_sol,time,index):
	plot = []
	plt.gca().set_aspect('equal', adjustable='box')
	plt.title("Time = "+str(time))
	plt.grid(c="#cccccc")
	rings = np.split(rings,critical_values)
	for i in range(0,int(len(mass_sol[0])/5)):
		plt.plot(mass_sol[:,0+5*i],mass_sol[:,1+5*i],c="g")
		plt.plot(np.nan,np.nan)

	#only looping over 6 values so no need for numpy loop
	for i,ring in enumerate(rings):
		if(len(ring)):
			plot.append(plt.plot(ring[:,0],ring[:,1],".",c=colours[i])[0])
	masses = np.reshape(mass_sol[index],(-1,5))
	plot.append(plt.plot(masses[:,0],masses[:,1],"o",c="k")[0])
	plt.show()	
	#returns plot to allow live plotting
	return plot

#plots a ring set with individual colours for each ring
def plot_ring_live(critical_values,rings,masses,graph):
	rings = np.split(rings,critical_values)
	# #only looping over 6 values so no need for numpy loop
	for i,ring in enumerate(rings):
		if(len(ring)):
			graph[i].set_xdata(ring[:,0])
			graph[i].set_ydata(ring[:,1])
			graph[i].set_mfc(colours[i])
	graph[len(graph)-1].set_xdata(masses[:,0])
	graph[len(graph)-1].set_ydata(masses[:,1])
	plt.draw()
	plt.pause(0.0001)

#live plots the full simulation
def plot_live_full(vals,ring_sol,mass_sol, dt):
	this_rings = np.reshape(ring_sol[0],(-1,4))
	plt.gca().set_aspect('equal', adjustable='box')
	plt.ion()
	t = 0
	graph = plot_ring_set(vals,this_rings,mass_sol,t,0)
	plt.pause(0.3)
	for i in range(1,len(ring_sol)):
		plot_ring_live(vals,np.reshape(ring_sol[i],(-1,4)),np.reshape(mass_sol[i],(-1,5)),graph)
		t = t+dt
		plt.title("Time = "+str(t))

#calculates the differential step for the test masses
def g(masses, rings):
	x = rings[:,0]
	y = rings[:,1]
	toReturn = np.zeros_like(rings)
	for mass in masses:
		delta_x = x-mass[0]
		delta_y = y-mass[1]
		r = np.sqrt(np.power(delta_x,2)+np.power(delta_y,2)+epsilon**2) #smoothing applied
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
		r3 = np.clip(r3,0,1e4) #to avoid infinities from self interactions
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
	full_set_f = np.append(rings_f,masses_f)
	full_length = len(full_set_f)

	dt = totalTime/noOfSteps
	t = np.linspace(0,totalTime,noOfSteps)
	sol = [full_set_f.flatten()]


	integrator = integrate.ode(ode_step_rev).set_integrator("dopri5")
	integrator.set_initial_value(full_set_f,0.0)
	integrator.set_f_params((ring_no))
	while integrator.successful() and integrator.t<totalTime:
		sol = np.append(sol,[integrator.integrate(integrator.t+dt)],axis=1)
	sol = np.reshape(sol,(-1,full_length))
	ring_sol = sol[:,:ring_no]
	mass_sol = sol[:,ring_no:]
	if (timeToPlot ==0):
		plot_live_full(vals,ring_sol,mass_sol,dt)
	else:
		index = round(timeToPlot/dt)
		plot_ring_set(vals,np.reshape(ring_sol[index],(-1,4)),mass_sol,timeToPlot,index)

	if (saveName!= ""):
		np.savetxt(saveName,sol,header=str(totalTime)+'\t' +str(noOfSteps))


# epsilon = 0.1 #smoothing length
# colours = ["#FF0000","#cc0040","#990080","#6600BF","#3300FF","#FFFFFF","#000000"]
# particle_density = 5
# #masses = [[0.0,0.0,0.0,0.0,1.0],[0.0,20.0,0.31,0.0,2.0]]
# #masses = [[0.0,0.0,0.0,0.0,1.0],[0.0,10.0,0.447,0.0,2.0]]
# masses = [[0.0,0.0,-0.15,0.0,1],[-30,-30,0.15,0.0,3]]
# #masses = [[0.0,0.0,0.0,0.0,1.0]]
# ring_set = create_ring_set([[2,12],[3,18],[4,24],[5,30],[6,36]],masses[0][:4])

# totalTime = 300
# noOfSteps = 300
# timeToPlot = 1
# #fileName = str(masses)+'t='+str(totalTime)+'.txt'
# fileName = 'test.txt'


# indiv_sim(masses,ring_set,totalTime,noOfSteps,timeToPlot,fileName)

# #print(sol[len(sol)-1])
# #sol = integrate.odeint(ode_step,full_set_f,t,args=(ring_no,))