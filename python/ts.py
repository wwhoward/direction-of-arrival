# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 15:55:34 2020

@author: Will
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy as sc
import sk_dsp_comm
from sk_dsp_comm import digitalcom
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
os.environ["KMP_DUPLICATE_LIB_OK"] = "0"   # sometime make "1" for Mac 


import math
import hdf5storage as h5
import random
import tensorflow as tf
import keras.backend.tensorflow_backend as KTF
from keras import regularizers
from keras.models import Model, Sequential
from keras.layers import Input, Dense, Lambda, Reshape, Conv1D, Conv2D,\
        AveragePooling2D,Flatten, Dropout, SimpleRNN, LSTM, concatenate, Layer
from keras.layers.normalization import BatchNormalization
from keras.callbacks import EarlyStopping
from numpy import ma
from matplotlib import cm as CM
import keras
from keras import callbacks 
from datetime import datetime
from keras.models import model_from_json
from numpy import linalg as la
import time



# Initialize Parameters
par = dict()

# Path Configuration
par["totalPaths"] = 3
par["forceMulti"] = False
par["forcePath"] = [1, 1, 1] # If not empty, make sure length(forcePath) = totalPaths
par["dimensions"] = 1

# Signal parameters
par["signalLength"] = 2**13
par["windowLength"] = 2**10
par["mod"] = 'QPSK'
par["interblock"] = (5, 5)
par["pathDelay"] = (0, 0)
par["minSep"] = np.pi/10
par["blocks"] = 3
par["blockSweep"] = (1,10)
par["snr"] = 20
par["snrSweep"] = (25,30)

# Estimation & Statistics
par["res"] = 2**6 / np.pi
par["recom"] = 'rnd'
par["scrub"] = True
par["sampling"] = 1

# Simulation
par["trials"] = 10**4
par["runtype"] = 'single'

# Saving
par["saveFlag"] = True
par["saveLight"] = False
par["saveName"] = "test"
par["savePlots"] = False

#==============================================================================
def AssignPaths(par):
    # Creates path object based on inputs:
    # par.K: total number of paths
    # par.type: '1d or '2d'
    
    # outputs path object with parameters
    # paths.sources : total number of sources
    # paths.multi : vector showing how many paths per source
    # paths.AoA : angles of arrival, random on the uniform sphere
    # paths.signal_vector(k).path(p) : signal vector for the p'th arrival of the k'th signal
    
    k = par["totalPaths"]
    paths = dict()
    paths["sources"] = 0
    
    if par["forceMulti"] & (par["forcePath"]==[]):
        if k==1:
            options = ([1])
            r=0
        elif k==2:
            options=([1,1])
            r=1
        elif k==3:
            options=([1,1,2],[1,1,1])
            r=np.random.randint(0,2)
        elif k==4:
            options=([1,1,2,3],[1,1,2,2],[1,1,1,2],[1,1,1,1])
            r=np.random.randint(0,4)
        elif k==5:
            options=([1,1,2,3,4],[1,1,1,2,3],[1,1,1,1,2],[1,1,1,1,1],[1,1,2,2,3],[1,1,1,2,2])
            r=np.random.randint(0,6)
    elif par["forcePath"] == []:
        if k==1:
            options = ([1])
            r=0
        elif k==2:
            options=([1,2],[1,1])
            r=np.random.randint(0,2)
        elif k==3:
            options=([1,2,3],[1,1,2],[1,1,1])
            r=np.random.randint(0,3)
        elif k==4:
            options=([1,2,3,4],[1,1,2,3],[1,1,2,2],[1,1,1,2],[1,1,1,1])
            r=np.random.randint(0,5)
        elif k==5:
            options=([1,2,3,4,5],[1,1,2,3,4],[1,1,1,2,3],[1,1,1,1,2],[1,1,1,1,1],[1,1,2,2,3],[1,1,1,2,2])
            r=np.random.randint(0,7)        
    else:
        options = ([par["forcePath"]])
        r = 0
        
    Q = options[r]
    if par["totalPaths"]==1:
        paths["sources"] = 1
    else:
        paths["sources"] = np.int(max(Q))
    
    paths["multi"] = np.zeros([paths["sources"],1])
    if par["totalPaths"] == 1:
        paths["multi"][0] = 1
    else:
        for source in range(1,paths["sources"]+1):
            #print(source)
            paths["multi"][source-1]= Q.count(source)
        
    # Lastly, assign AoA angles for each path
    if par["dimensions"] == 1:
        paths["AoA"] = np.zeros([par["totalPaths"],2])
        Azi = AssignAzi(par)
        minSepFlag=0
        
        if np.size(Azi) !=1:
            if np.amin(np.diff(np.transpose(Azi))) < par["minSep"]:
                minSepFlag = 1
        
        while (not(par["minSep"]==[]) and (minSepFlag==1)):
            Azi = AssignAzi(par)
            if np.amin(np.diff(np.transpose(Azi))) >= par["minSep"]:
                break
        
        counter=0
        paths["signalVector"] = list()
        for source in range(0,paths["sources"]):
            temp_dict = dict()
            temp_dict["path"] = np.zeros([np.int(paths["multi"][source][0]), 6, 1], dtype=np.complex64)
            for path in range(0,np.int(paths["multi"][source][0])):
                
                azi = Azi[0]
                Azi = Azi[1:]
                ele = np.pi/2 + 0.01*np.random.rand()
                
                paths["AoA"][counter] = [ele,azi]
                counter=counter+1
                
                [__1, __2, A] = Manifold([ele,azi], [np.pi/2*np.random.rand(), np.pi*2*np.random.rand()])
                temp_dict["path"][path,:,] = A
            paths["signalVector"].append(temp_dict)
    
                
                
    else:
        print("Only one dimension currently supported! Exiting...")
        sys.exit()
        
    return paths


def Transmit(paths, par): 
    msg = np.zeros([paths["sources"], par["signalLength"]])
    signal = dict()
    for source in range(0, paths["sources"]):
        msg[source] = np.random.randint(0,2,par["signalLength"])
        
    if par["mod"]=='bin':
        signal["tx"] = msg
        
    elif par["mod"] == 'QPSK':
        signal["tx"] = np.reshape(np.asarray(sk_dsp_comm.digitalcom.QPSK_bb(512, 100)[0][:par["signalLength"]]), [1, par["signalLength"]])
    

    
    return signal

def Recieve(signal, paths, par): # Need to add SNR still!!!
    tx = signal["tx"]
    rx = np.zeros([par["blocks"], 6, par["windowLength"]], dtype=complex)
    R = np.zeros([6,6], dtype=complex)
    
    for block in range(0, par["blocks"]):
        for source in range(0, paths["sources"]):
            for path in range(0, np.int(paths["multi"][source])):
                if paths["multi"][source] == 1:
                    path_gain = 1
                else:
                    path_gain = np.sqrt(0.5) * (np.random.normal()+1j*np.random.normal())
                
                a = paths["signalVector"][source]["path"][path,:,:]
                if par["pathDelay"][0] == par["pathDelay"][1]:
                    path_delay = par["pathDelay"][0]
                else:
                    path_delay = np.random.randint(par["pathDelay"][0], par["pathDelay"][1],1)
                if par["interblock"][0] == par["interblock"][1]:
                    block_delay = par["interblock"][0]
                else:
                    block_delay = np.random.randint(par["interblock"][0], par["interblock"][1],1)
                window = par["windowLength"]
                sig = tx[source, (block)*block_delay+(path)*path_delay : window + (block)*block_delay+(path)*path_delay]
                rx[block, :,:] = path_gain*a*sig + np.squeeze(rx[block, :,:])

    noise_power = 10**(-par["snr"]/10)
    noise = np.random.normal(0, noise_power, [par["blocks"], 6, par["windowLength"]])
    signal["rx"] = rx + noise
    R = 1/par["windowLength"] * (np.squeeze(signal["rx"][1,:,:]) * np.asmatrix(np.squeeze(signal["rx"][1,:,:])).H)
    signal["R"] = R
    return signal

def Smooth(signal, par): 
    R_ts = np.zeros([6,6], dtype=complex)
    r_ts = np.zeros([par["blocks"], 6, 6], dtype=complex)
    for block in range(0, par["blocks"]):
        r_ts[block, :, :] = 1/par["windowLength"] * signal["rx"][block, :, :]*np.asmatrix(signal["rx"][block, :, :]).H
        R_ts[:,:] = (1/par["blocks"] * r_ts[block, :, :]) + R_ts[:,:]
    signal["ts"] = R_ts
    
    return signal

def DrMusic(R, par):
    est = dict()
    if par["dimensions"] == 1:
        est["doaFunction"] = np.zeros([6, 2, np.int(2*np.pi*par["res"])])
        est["spectrum"] = np.zeros([np.int(2*np.pi*par["res"])])
        est["phi"] = np.linspace(0, 2*np.pi, np.int(2*np.pi*par["res"]), endpoint=True)
        for ph in range(0, np.size(est["phi"], 0)):
            [est["doaFunction"][:,:,ph], __1, __2] = Manifold([np.pi/2, est["phi"][ph]], [np.pi/4, 0])
        
        [eigval, eigvect] = np.linalg.eig(R)
        noiseSpace = np.delete(eigvect, np.argsort(eigval)[-par["totalPaths"]:], 1)
        
        for ph in range(0, np.size(est["phi"], 0)):
            [eigvals, __1] = np.linalg.eig(np.asmatrix(est["doaFunction"][:,:,ph]).H * noiseSpace * np.asmatrix(noiseSpace).H * est["doaFunction"][:,:,ph])
            est["spectrum"][ph] = 1/min(np.real(eigvals))
        
        est["peaks_idx"] = sc.signal.find_peaks(est["spectrum"])[0]
        est["peaks_val"] = np.zeros([np.size(est["peaks_idx"])])
        est["peaks_azi"] = np.zeros([np.size(est["peaks_idx"])])
        for idx, val in enumerate(est["peaks_idx"]):
            est["peaks_val"][idx] = est["spectrum"][val]
            est["peaks_azi"][idx] = est["phi"][val]
        
        
    return est

def CalcStats(est, signal, par, paths):
    stats = dict()
    if par["dimensions"]==1:
        azi_u = paths["AoA"][:,1]
        azi_uh = est["peaks_azi"]
        stats["errAzi"] = np.zeros([np.size(azi_u)])
        stats["correctedErrAzi"] = np.zeros([np.size(azi_u)])
        stats["recomFlag"] = np.zeros([np.size(azi_u)])
        recom_idx = []
        norecom_idx = []
        
        closestValue = np.zeros([np.size(azi_u)])
        for peak in range(0, np.size(azi_u)):
             if azi_uh.size  > 0:
                 stats["errAzi"][peak] = np.amin(abs(azi_u[peak] - azi_uh.T))
                 closestIndex = np.argmin(abs(azi_u[peak] - azi_uh.T))
                 closestValue[peak] = azi_uh[closestIndex]
                 azi_uh = np.delete(azi_uh, closestIndex)
                 stats["correctedErrAzi"][peak] = min((2*np.pi)-abs(closestValue[peak]-azi_u[peak]), abs(closestValue[peak]-azi_u[peak]))
                 stats["recomFlag"][peak] = 0
                 norecom_idx.append(peak)
             else:
                 if par["recom"] == 'rnd':
                     stats["errAzi"][peak] = np.pi * np.random.rand()
                     
                 elif par["recom"] == 'max':
                     stats["errAzi"][peak] = np.pi
                 stats["correctedErrAzi"][peak] = stats["errAzi"][peak]
                 stats["recomFlag"][peak] = 1
                 recom_idx.append(peak)
                 
    
    stats["azi_u"] = azi_u
    stats["azi_uh"] = closestValue
    
    stats["azi_mse"] = sum(stats["errAzi"]**2)/np.size(stats["errAzi"])
    stats["azi_rmse"] = math.sqrt(stats["azi_mse"])
    
    stats["azi_err_deg"] = stats["errAzi"]/np.pi * 180
    stats["azi_mse_deg"] = sum(stats["azi_err_deg"]**2)/np.size(stats["azi_err_deg"])
    stats["azi_rmse_deg"] = math.sqrt(stats["azi_mse_deg"])
    
    stats["azi_msce"] = sum(stats["correctedErrAzi"]**2)/np.size(stats["correctedErrAzi"])
    stats["azi_rmsce"] = math.sqrt(stats["azi_msce"])
    
    stats["corrected_azi_err_deg"] = stats["correctedErrAzi"]/np.pi * 180
    stats["azi_msce_deg"] = sum(stats["corrected_azi_err_deg"]**2)/np.size(stats["corrected_azi_err_deg"])
    stats["azi_rmsce_deg"] = math.sqrt(stats["azi_msce_deg"])
    
    
    
    stats["recom_percent"] = sum(stats["recomFlag"]) / np.size(stats["recomFlag"])
    
    stats["azi_norecom_mse"] = sum(stats["errAzi"][norecom_idx]**2)/np.size(stats["errAzi"][norecom_idx])
    stats["azi_norecom_rmse"] = math.sqrt(stats["azi_norecom_mse"])
    
    stats["azi_norecom_err_deg"] = np.pi/180 * stats["errAzi"][norecom_idx]
    stats["azi_norecom_mse_deg"] = sum(stats["azi_norecom_err_deg"]**2)/np.size(stats["azi_norecom_err_deg"])
    stats["azi_norecom_rmse_deg"] = math.sqrt(stats["azi_norecom_mse_deg"])
    
             
    return stats    
    
def SNR_Plot(est, paths, par, plt_title = "DR-MUSIC Spatial Spectrum"):
    if par["dimensions"] == 1:
        plt.figure()        
        plt.vlines(paths["AoA"][:,1], np.amin(est_ts["spectrum"]), np.amax(est_ts["spectrum"]), linestyles='dotted')
        plt.plot(est["peaks_azi"][:min(par["totalPaths"], np.size(est["peaks_azi"]))], est["peaks_val"][:min(par["totalPaths"], np.size(est["peaks_azi"]))], 'go', markersize = 12)
        plt.plot(est["phi"], est["spectrum"])
        plt.title(plt_title)
        plt.legend(["True", "Estimate"])
        plt.xlabel("Radians")
    
    
    
def Manifold(aoa, pol):
    gtht = np.array([[np.cos(aoa[0]) * np.cos(aoa[1]), -np.sin(aoa[1])], 
                    [np.cos(aoa[0]) * np.sin(aoa[1]),  np.cos(aoa[1])], 
                    [-np.sin(aoa[0])  ,  0], 
                    [-np.sin(aoa[1]), -np.cos(aoa[0]) * np.cos(aoa[1])], 
                    [np.cos(aoa[1]),  np.sin(aoa[0]) * np.sin(aoa[1])], 
                    [0,  np.sin(aoa[0])]])
    gpol = np.array([[np.sin(pol[0]) * np.exp(1j*pol[1])], [np.cos(pol[0])]])
    
    g = np.dot(gtht, gpol)
    
    return gtht, gpol, g

def AssignAzi(par):
    
    if par["forcePath"] == []:
        k = par["totalPaths"]
    else:
        k = len(par["forcePath"])
        
    azi = np.zeros([k,1])
    for source in range(0,k):
        azi[source] = np.random.rand()*2*np.pi
    
    azi = np.sort(azi,0)
    
    return azi

# def detect_peaks(image):
#     """
#     Takes an image and detect the peaks usingthe local maximum filter.
#     Returns a boolean mask of the peaks (i.e. 1 when
#     the pixel's value is the neighborhood maximum, 0 otherwise)
#     """

#     # define an 8-connected neighborhood
#     neighborhood = generate_binary_structure(2,2)

#     #apply the local maximum filter; all pixel of maximal value 
#     #in their neighborhood are set to 1
#     local_max = maximum_filter(image, footprint=neighborhood)==image
#     #local_max is a mask that contains the peaks we are 
#     #looking for, but also the background.
#     #In order to isolate the peaks we must remove the background from the mask.

#     #we create the mask of the background
#     background = (image==0)

#     #a little technicality: we must erode the background in order to 
#     #successfully subtract it form local_max, otherwise a line will 
#     #appear along the background border (artifact of the local maximum filter)
#     eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

#     #we obtain the final mask, containing only peaks, 
#     #by removing the background from the local_max mask (xor operation)
#     detected_peaks = local_max ^ eroded_background

#     return detected_peaks





#==============================================================================
    
if par["runtype"]=='single':
    
    par["saveFlag"]=0   # No reason to save single runs since they're just for debugging
    
    paths = AssignPaths(par)

    signal = Transmit(paths, par)

    signal = Recieve(signal, paths, par)
    
    signal = Smooth(signal, par)

    est_ts = DrMusic(signal["ts"], par)
    est_nts = DrMusic(signal["R"], par)
    
    stats_ts = CalcStats(est_ts, signal, par, paths)
    stats_nts = CalcStats(est_nts, signal, par, paths)
    
    SNR_Plot(est_ts, paths, par)
    
    # plt.figure()
    # plt.plot(est_ts["spectrum"])
    # plt.vlines((np.int(par["res"]))*paths["AoA"][:,1], np.amin(est_ts["spectrum"]), np.amax(est_ts["spectrum"]), linestyles='dotted')






















