# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:26:36 2019

est1D(X_t, input_type, A, K) Estimate azimuthal angle of arrival from recieved signal
INPUT:  X_t         6xN recieved signal data OR 6x6 covariance matrix
        input_type  string, = {'rx', 'cov'}
        A           RESx6 where RES is desired resolution
        K           Number of signals incedent on sensor

OUTPUT: ret         Dictionary of 'u_h' and 'B'
        u_h         Estimated DoA, as indicies of A (unsorted)
        B           Spatial Spectrum - output from MUSIC operation

Assumptions: 
        A RES of x corresponds to degree increments of 360/RES, i.e. RES = 36 => degree increments of 10
        Current support is for 1 degree of freedom, and the intent is for this to be azimuth
        

plot(B, u_h, u_true) Plot given data
INPUT:  B           MUSIC spatial spectrum
        u_h         Estimated DoA, as indicies of B
        u_true      True DoA, as indicies of B
        
OUTPUT: None

Assumptions: 
        None
        
@author: {wwhoward} @vt.edu
Wireless@VT
"""
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

import time
import os
import math
import hdf5storage as h5
os.environ["KMP_DUPLICATE_LIB_OK"] = "1"   # sometime make "1" for Mac 
import random
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import tensorflow as tf
import keras.backend.tensorflow_backend as KTF
from keras import regularizers
from keras.models import Model, Sequential
from keras.layers import Input, Dense, Lambda, Reshape, Conv1D, Conv2D,\
        AveragePooling2D,Flatten, Dropout, SimpleRNN, LSTM, concatenate, Layer
from keras import backend as K
from keras.layers.normalization import BatchNormalization
from keras.callbacks import EarlyStopping
import matplotlib.pyplot as plt
from numpy import ma
import scipy.io as sio
from matplotlib import cm as CM
import keras
from keras import callbacks 
from datetime import datetime
import scipy.io as sio
from keras.models import model_from_json

def est_MUSIC_1D(X_t, input_type, A, K):
    
    X_t = np.asmatrix(X_t)
    A   = np.asmatrix(A)
    
    #Decide what type of input data we're given, and output a covariance matrix
    if input_type == 'rx':
        Rxx = X_t * X_t.H 
    elif input_type == 'cov':
        Rxx = X_t
    else:
        print("Input type invalid. Current options are: 6xN recieved signal vector 'rx' OR 6x6 covariance matrix 'cov'. ")
        return
    

    w, v = la.eig(Rxx)              #Eigendecomposition of given covariance matrix
    
    w_max = np.argsort(w)[-K:]      #Find the indicies of the K largest eigenvalues in V
    
    Un = np.delete(v, w_max,1)      #Eliminate the K largest eigenvectors of V
        
    Unn = Un*Un.H                   #Construct noise subspace from (6-K) smallest eigenvectors of Rxx
    
    B = np.zeros(A.shape[0])        #Initialize a variable to hold MUSIC values. We'll later look for the peaks of this
    
    azi_n = A.shape[0]              #Amount of angles to search
    
    for azi in range(0, azi_n):     #Loop over possible DoA's, save results to B
        arr = np.asmatrix(A[azi]).T
        b = arr.H * Unn * arr
        B[azi] = 1/abs(b)
        
    u_h = np.argsort(B)[-K:]        #Sort B and save indicies of K largest values
        
    
    ret = dict()                    #Construct an object to hold out output values, currently outputting the estimated angles & "spatial spectrum"
    ret['u_h'] = u_h
    ret['B'] = B
    return ret

def est_MUSIC_2D(X_t, input_type, A, K):
    print("Two dimensional estimator not yet implemented. ")
    return
    
def plot(B, u_h, u_true):
    plt.figure()
    u_h.sort()                                          #Sort so we plot in the right order
    u_true.sort()
    
    for est, tru in zip(u_h, u_true):
        plt.plot(est, B[est], 'go', markersize = 12)
        plt.plot(tru, B[tru], 'yo', markersize = 6)
    plt.plot(B)
    plt.axis
    plt.xlabel("Azimuth in 10's of Degrees")
    plt.title('MUSIC Spatial Specrtrum')
    plt.legend(['Estimate','True'])
    plt.show()
    plt.close()
    return

def est_ML_1D(dataID="", plot_flag=1, train_flag=1, test_flag=1, eval_flag=1):
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    early_stop = EarlyStopping(monitor='val_loss', min_delta=0, patience=5, verbose=1, mode='auto')
    # ======================== global setting ======================== 
    #dataID   = '01'   
    methodID = '512X5/'
    Nsig     = 1          # Number of source signal
    INPUTLEN = 30
    print(Nsig)
    print(dataID)
    print(methodID)

    if not os.path.isdir( dataID):
        os.mkdir( dataID)
        if not os.path.isdir( dataID+methodID):
            os.mkdir( dataID+methodID)
    
    note = '1% batch size to all (default 1%/16. ~128)'   
    TEXT    = "low dataset" 
    epochs   = 400   

    # ------------- load dataset ------------------
    t = h5.loadmat(dataID+'Xtrain.mat')
    Xtrain = t['Xtrain'].T
    t = h5.loadmat(dataID+'Xtest.mat')
    Xtest = t['Xtest'].T
    t = h5.loadmat(dataID+'Ytrain.mat')
    Ytrain = t['Ytrain'].T
    t = h5.loadmat(dataID+'Ytest.mat')
    Ytest = t['Ytest'].T
    [sizeTrain,_] = np.shape(Xtrain)
    Nval = int(sizeTrain*0.30)
    Xval = Xtrain[0:Nval,:];
    Yval = Ytrain[0:Nval,:];
    Xtrain = Xtrain[Nval+1:,:];
    Ytrain = Ytrain[Nval+1:,:];
    Xtrain = Xtrain[:,:,np.newaxis]
    Xtest  = Xtest[:,:,np.newaxis]
    Xval   = Xval[:,:,np.newaxis]        
    Y_test  = Ytest
    Y_train  = Ytrain
    Y_val  = Yval
    errS = np.zeros((Nsig,2))
    rmseS = np.zeros((Nsig,2)) 
    
    # ---- config batch size regarding training size ----
    t = sizeTrain/100
    tn = math.log(t)/math.log(2)
    if tn > 6:
        batch_size = 2 ** math.floor(tn)
    else:
        batch_size = 2 ** 8
    batch_size = int(batch_size/16)  
    print(batch_size)
    
    
    tic = time.time()
    with open(dataID+methodID+"outPut_TF_DNN.txt", "a") as text_file:
        if train_flag:
            text_file.write( "=== " + note + " === \n" )
            text_file.write( "--- caseID %s  begin --- \n" %(dataID))
            text_file.write( "--- local time  " + dt_string + " --- \n" )    
            # ------------ config neural network ------------
            nn_input  = Input((INPUTLEN,1))                
            nn_output = Flatten()(nn_input)
            nn_output = Dense(512,activation='relu')(nn_output)
            nn_output = Dense(512,activation='relu')(nn_output)             
            nn_output = Dense(512,activation='relu')(nn_output) 
            nn_output = Dense(512,activation='relu')(nn_output)
            nn_output = Dense(512,activation='relu')(nn_output)
            nn_output = Dense(2*Nsig,activation='linear')(nn_output)  #----------- by default only estimate elevation and azimuth ---
            nn = Model(inputs=nn_input,outputs=nn_output)
            
            nn.compile(optimizer='adam', loss='mse',metrics=[dist])
            nn.summary()
        
            train_hist = nn.fit(x=Xtrain,y=Y_train,\
                                batch_size = batch_size ,epochs = epochs ,\
                                validation_data=(Xval, Y_val), \
                                shuffle=True,\
                                callbacks=[early_stop])
                
        
        # ----test -----------
        if test_flag:
            Ypred = nn.predict( Xtest)
            sio.savemat(dataID+methodID+'Ypred.mat', {'Ypred':Ypred})
        
        # ---- evaluation ----------
        if eval_flag: 
            errs = abs(Y_test-Ypred) * 180/math.pi
            err = np.mean(errs,axis=0)    
            errS = err
            rmse = np.sqrt(    np.mean( ( Y_test - Ypred) ** 2 , axis=0 )   )   * 180/math.pi
            rmseS = rmse    
            print( "--- MAE: --- %s" %(err))
            print( "--- RMSE: --- %s" %(rmse))
            text_file.write( " layer [ 512 *3 ] \n")     
            text_file.write( " batch size after D16 %d \n"%batch_size)  
            text_file.write( " test error %s (deg) \n" %(err))
            text_file.write( "rmse %s (deg) \n" %(rmse))
            toc =  time.time()
            timeCost = toc - tic
            print( "--- Total of %s seconds ---" %(timeCost))
            text_file.write( " timeCost %s \n" %(timeCost))
        
        #Plotters
        if plot_flag: 
            for i in range(0,2*Nsig):
                plt.figure(i)
                plt.scatter( Y_test[:,i], Ypred[:,i],s = 1, facecolors='none',edgecolors='b')
                plt.title(' est vs ground')
                plt.ylabel('est')
                plt.xlabel('ground')
                plt.grid(True)  
                plt.savefig(dataID+methodID + 'scatter-TF-%d.png'%(i) )
                plt.show()
                plt.clf()
                
            plt.figure(11)
            plt.plot(train_hist.history['dist'])
            plt.plot(train_hist.history['val_dist'])
            plt.title('distance')
            plt.ylabel('distance')
            plt.xlabel('epoch')
            plt.grid(True)  
            plt.legend(['train', 'validate'])
            plt.savefig(dataID+methodID+'hist_dist.png')
            plt.show()
            plt.clf()
        
            plt.figure(13)
            plt.plot(train_hist.history['loss'])
            plt.plot(train_hist.history['val_loss'])
            plt.title('model loss')
            plt.ylabel('loss')
            plt.xlabel('epoch')
            plt.grid(True)  
            plt.legend(['train', 'validate'])
            plt.savefig(dataID+methodID+'hist_loss.png')
            plt.show()
            plt.clf()
        
        # save model to file
        model_json = nn.to_json()
        with open(dataID+methodID+"model_TF.json", "w") as json_file:
            json_file.write(model_json)
            nn.save_weights(dataID+methodID+"model_TF.h5")
            print("Saved model to disk")
    return

def temporal_smooth(signal, blocks, snap, inter, delay, multi=[1,1], aoa=[np.pi/2, 3*np.pi/2]):
    # If you call this fn we assume that the prerequisites are met, i.e. signal is long enough. 
    # Assume the signal is already noisy
    rx = np.zeros(blocks, 6, snap)
    for b in range(blocks):
        for k in range(np.length(multi)):
            for p in range(multi(k)):
                if p != 1:
                    h = np.sqrt(0.5)*(np.random.randn + 1j*np.random.randn)
                else:
                    h = 1;
                [_1, _1, A] = manifold(aoa, [np.pi/2 * np.random.rand, 2*np.pi*np.random.rand - np.pi])
                sig = signal[k, 1+(b-1)*inter+(p-1)*delay : snap + (b-1)*inter+(p-1)*delay]
                rx[b,:,:] = h*A*sig + np.squeeze(rx[b,:,:])
    
    R = 1/snap * np.squeeze(rx[1,:,:]) * np.matrix.H(np.squeeze(rx[1,:,:])) #Normal covariance
    
    smooth = np.zeros((6,6))
    r = np.zeros((blocks, 6, 6))
    
    for b in range(blocks):
        r[b, :,:] = 1/snap * np.squeeze(rx[b,:,:]) * np.matrix.H(np.squeeze(rx[b,:,:]))
        smooth[:,:] = 1/blocks * np.squeeze(r[b,:,:]) + smooth[:,:] # Temporally smoothed covariance
    
    return smooth, R


#Local functions

def polar2cart(theta, phi):
    r1 = tf.math.sin(theta) * tf.math.cos(phi)
    r2 = tf.math.sin(theta) * tf.math.sin(phi)
    r3 = tf.math.cos(phi)
    return r1,r2,r3


# distance is defined as the error distance over a unit sphere
def dist(y_true, y_pred):
    Nsig = 1 ################# need to be reset  according  to number of source sginal ##########################
    S = 0
    for i in range(0,Nsig):
        R_true = polar2cart(y_true[:,2*i], y_true[:,2*i+1])
        R_pred = polar2cart(y_pred[:,2*i], y_true[:,2*i+1])
        s = tf.sqrt( tf.square (tf.abs(R_pred[0]-R_true[0])) + tf.square (tf.abs(R_pred[1]-R_true[1]))  + tf.square (tf.abs(R_pred[2]-R_true[2])))
        S = S + s       
    return   S

def manifold(aoa, pol):
    gtht = np.array([np.cos(aoa[1]) * np.cos(aoa[2]), -np.sin(aoa[2])], 
                    [np.cos(aoa[1]) * np.sin(aoa[2]),  np.cos(aoa[2])], 
                    [-np.sin(aoa[1])  ,  0], 
                    [-np.sin(aoa[2]), -np.cos(aoa[1]) * np.cos(aoa[2])], 
                    [np.cos(aoa[2]),  np.sin(aoa[1]) * np.sin(aoa[2])], 
                    [0,  np.sin(aoa[1])])
    gpol = np.array([np.sin(pol[1]) * np.exp(1j*pol[2])], [np.cos(pol[1])])
    
    g = gtht * gpol
    
    return gtht, gpol, g
    
    
    
    