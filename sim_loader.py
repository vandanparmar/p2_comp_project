import numpy as np
from scipy import integrate
from sim_plotting import *

#calculates the distance between a given pair of masses
def index_min_distance(mass_sol):
	r1 = np.sqrt(np.power(mass_sol[:,0],2)+np.power(mass_sol[:,1],2))
	r2 = np.sqrt(np.power(mass_sol[:,5],2)+np.power(mass_sol[:,6],2))
	distances = r1+r2
	toReturn = np.argmin(distances)
	return toReturn

def closest_particle(mass_sol,ring_sol):
	x = mass_sol[5]
	y = mass_sol[6]
	ring_sol = np.reshape(ring_sol,(-1,4))
	delta_x = ring_sol[:,0]-x
	delta_y = ring_sol[:,1]-y
	r = np.sqrt(np.power(delta_x,2)+np.power(delta_y,2))
	toReturn = np.min(r)
	return toReturn


def file_loader(file_name):
	toReturn = {}
	with open(file_name,'r') as f:
		toReturn['totalTime'],toReturn['noOfSteps'],toReturn['ring_no'] = f.readline().split('\t')
		toReturn['totalTime'] = int(toReturn['totalTime'][2:])
		toReturn['noOfSteps'] = int(toReturn['noOfSteps'])
		toReturn['ring_no'] = int(toReturn['ring_no'][:-1])
		toReturn['vals'] = eval(f.readline()[2:-1])
		sol = np.loadtxt(f)
		toReturn['ring_sol'] = sol[:,:toReturn['ring_no']]
		toReturn['mass_sol'] = sol[:,toReturn['ring_no']:]
		toReturn['dt'] = toReturn['totalTime']/toReturn['noOfSteps']
	return toReturn

def plot_loaded_full(full_sol):
	plot_live_full(full_sol['vals'],full_sol['ring_sol'],full_sol['mass_sol'],full_sol['dt'])

def plot_loaded_time(full_sol,timeToPlot):
	index = int(round(timeToPlot/full_sol['dt']))
	plot_ring_set(full_sol['vals'],np.reshape(full_sol['ring_sol'][index],(-1,4)),full_sol['mass_sol'],timeToPlot,index)

def time_of_bridge(full_sol):
	t = np.linspace(0,full_sol['totalTime'],full_sol['noOfSteps'])
	t_closest = t[index_min_distance(full_sol['mass_sol'])]
	t_bridge = 0
	for i in range(0,full_sol['noOfSteps']):
		t_bridge = t[i]
		distance = closest_particle(full_sol['mass_sol'][i,:],full_sol['ring_sol'][i,:])
		if (distance<0.5):
			break
	toReturn = t_closest - t_bridge
	return [toReturn,t_closest,t_bridge]


number = 100
filename = 'parabolic_equal_mass/parabolic_equal_mass_q_'+str(int(number))+'.txt'
full_sol = file_loader(filename)
times = time_of_bridge(full_sol)

#plot_loaded_full(full_sol)
plot_loaded_time(full_sol,145)
#plot_loaded_time(full_sol,times[1]+2*np.abs(times[1]-times[2]))

# timeToPlot = 109

	
# qs = np.linspace(8.5,13.5,51)
# time_diff = []
# for q in qs:
# 	filename = 'parabolic_equal_mass/parabolic_equal_mass_q_'+str(int(q*10))+'.txt'
# 	full_sol = file_loader(filename)
# 	time_diff.append(time_of_bridge(full_sol))
# 	print(q)

# plt.plot(qs,time_diff)
# plt.show()


# f = open('parabolic_orbit_approach_bridge_time.txt','w')
# f.write('q \t time \n')
# for i in range(0,len(qs)):
# 	toWrite = str(qs[i]) +'\t' + str(time_diff[i])+'\n'
# 	f.write(toWrite)
# f.close






