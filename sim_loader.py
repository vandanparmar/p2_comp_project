import numpy as np
from scipy import integrate
from sim_plotting import *


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


full_sol = file_loader('test.txt')

plot_live_full(full_sol['vals'],full_sol['ring_sol'],full_sol['mass_sol'],full_sol['dt'])