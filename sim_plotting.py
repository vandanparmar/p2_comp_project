import numpy as np
import matplotlib.pyplot as plt

colours = ["#FF0000","#cc0040","#990080","#6600BF","#3300FF","#FFFFFF","#000000"]


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
