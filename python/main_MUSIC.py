# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:32:44 2019
Tester function to demonstrate capability of MUSIC algorithm


@author: wwhoward, wireless@vt
"""
import scipy.io
import DOA
import numpy as np
from numpy import random as rnd

#Load testing data as covariance matrix
X = scipy.io.loadmat('DATA/sig/manifold_matrix.mat')['manifold_matrix']
A = scipy.io.loadmat('DATA/sig/manifold_eig.mat')['manifold_eig']

K = 3; #number of incedent signals

freq = rnd.randint(0,1514); #This is specific to the format of the test data, pick a random frequency

u_true = rnd.permutation(36)[:K] #Pick random DoA's to test on


X_t = np.asmatrix(np.zeros([6,6]))

for angle in u_true:
    X_t = X_t + np.asmatrix(X[freq, angle, :,:])
j = DOA.est_1D(X_t, 'cov',  A[:,freq,:], K)


DOA.plot(j['B'],j['u_h'], u_true[:K])



#Now try loading testing data in the "rx" form: a 6xN vector
X = scipy.io.loadmat('DATA/sig/rx_data')['y']
A = np.squeeze(scipy.io.loadmat('DATA/sig/rx_data')['A']).transpose()

K = 1; #Set this to the actual number of singals present

u_true = [13];

X_rx = X[:,:,u_true]

j = DOA.est_MUSIC_1D(X_rx, 'rx', A, K)
 
DOA.plot(j['B'],j['u_h'], u_true) 

#Use the ML package to estimate DoA's. 
#NOTE: This requires a differently formatted dataset than the previous examples. 

DOA.est_ML_2D('DATA/ML/01')