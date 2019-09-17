
"""
d = MinCurvature(xyz2D, d2D, xyzLine)

Function to compute the minimum curvature interpolation from a grid of points
to a set of points



Created on Mon Feb 13 15:38:09 2017

@author: dominiquef
"""
import numpy as np
import scipy as sp
from scipy.sparse.linalg import bicgstab
import matplotlib.pyplot as plt

## USER INPUTS ##
# Subsample the grid for testing
n = 5

# Load some pre-computed data for test
grid2D = np.loadtxt('Tests\\Grid2D.dat')
line2D = np.loadtxt('Tests\\Line2D.dat')

# Reformat the data grid locations
indx = np.random.random_integers(0, int(grid2D.shape[0]-1), int(grid2D.shape[0]/3))

locX = grid2D[indx, 0]
locY = grid2D[indx, 1]
data = grid2D[indx, 3]

## MINIMUM CURVATURE ALGO STARTS HERE ####
ndat = locX.shape[0]
A = np.zeros((ndat,ndat))


for i in range(ndat):

    r = (locX[i] - locX)**2. + (locY[i] - locY)**2.
    A[i, :] = r.T * (np.log((r.T + 1e-8)**0.5) - 1.)

# Solve system for the weights
w = bicgstab(A, data, tol=1e-6)

# Compute new solution
# Reformat the line data locations but skip every n points for test
xx = line2D[::n, 0]
yy = line2D[::n, 1]
data_line = line2D[::n, 3]

d_i2d_MinCurv = np.zeros(len(xx))

for i in range(len(xx)):

    r = (xx[i] - locX)**2. + (yy[i] - locY)**2.
    d_i2d_MinCurv[i] = np.sum( w[0] * r.T * ( np.log( (r.T + 1e-8)**0.5 ) - 1. ))

#%% PRINT DATA AND SOLUTION #####
plt.figure(figsize=(12,10))
plt.subplot(2,1,1)
plt.scatter(grid2D[:,0], grid2D[:,1], 15, c=grid2D[:,3])
plt.scatter(locX, locY, 30, c='r', marker='o')
plt.plot(xx,yy,c='b', lw=3)

plt.subplot(2,1,2)
plt.plot(line2D[:,0],line2D[:,3])
plt.plot(xx,d_i2d_MinCurv, marker='o')
plt.legend(['True','Minimum Curvature: n='+str(n)])
plt.show()
