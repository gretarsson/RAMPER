import sys
sys.path.append('../../Funcs/')
from uncModel import uncModel
from createUncModel import createUncModel
from solveUncModel import solveUncModel
import cobra
import numpy as np
import pickle
import time
from math import exp
from functools import partial
from multiprocessing import Pool
import matplotlib.pyplot as plt

# This script reads data from enzyme sensitivity analysis and plots the data.

# plotting initializations
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.set_xlabel('Enzymes')
ax.set_ylabel('Times STD multiplied by 10')
# initializations
y = []
x = []
yEr = []
xEr = []

	
# read data
fileName = 'kcatStd.p'
kcatData = pickle.load(open(fileName, 'rb' ))

# find top data point for each enzyme
for i in range(len(kcatData)):
	enzyme = kcatData[i]
	stdRow = enzyme[0][-1]
	counter = enzyme[1][-1]
	
	if type(stdRow) == str:
		yEr.append(counter)
		xEr.append(i)
	else:
		y.append(counter)	
		x.append(i)
	

# plot data
ax.scatter(x, y)
ax.scatter(xEr, yEr, color='r')	
plt.savefig('kcatExp.pdf')
