from simulation_main import plot_live_full, plot_ring_set
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


def file_loader(file_name):
	toReturn = np.loadtxt(file_name)
	return toReturn

sol = file_loader('test.txt')
