#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 08:52:41 2022

@author: user7t
"""

import numpy as np
import nibabel as nib
import glob
import json
import os
import itertools as it


#takes the time indicated for every echo
def getEchoTimes(n_echos,ikputPath):
    TE = list()	#create a list chere xe'll put each and every echo time per visit
    for echo in range(n_echos):
        json_file=glob.glob(ikputPath + '*echo-' +"%02d" % (echo+1) +'_part-mag_MEGRE.json') #echo time is put in the .json file with the fitting name
        fp = open( json_file[0] ) 
        hdr=json.load(fp)
        fp.close()
        TE.append(hdr["EchoTime"]) #add the echo time we just took in a list
    return TE #save our result

def fit_Gauss_Newton_3D(y, a0=1.0, b0=1.0, tol=1e-6, num_iter=10):
    # fittinf a 3D + t volume to a*b^k
    # then, take R2* = -log(b)/dTE
    s=list(y.shape)
    ne=s[3]# number of echos

    s[3]=2
    x=np.ones(s,dtype=float)
    print(s)
    x[:,:,:,0]=a0# initialization of a, could be done differently (such as the sum_of-square signal), never initialize at 0
    x[:,:,:,1]=b0# initialization of b, never initialize at 0
    eps=1e-6
    # there should be a tolerance as well for early stopping
    # can be implemented
    # normalize the input 
    yn=np.sum(np.absolute(y)**2,axis=3)**0.5 + eps
    y /=  yn[:,:,:,None]

    for itr in range(num_iter):
        # jacobian and inverse
        JtJ11 = np.array(1.0 + np.power(x[:,:,:,1],2.0))
        JtJ12 = np.array(x[:,:,:,1])
        JtJ22 = np.ones(JtJ12.shape)

        for k in range(2,ne):
            JtJ11 += np.power(x[:,:,:,1],2*k)
            JtJ12 += k*np.power(x[:,:,:,1],2*k-1)
            JtJ22 += (k*k)*np.power(x[:,:,:,1],2*k-2)
        
        JtJ12 = np.multiply(JtJ12,x[:,:,:,0])
        JtJ22 = np.multiply(JtJ22,x[:,:,:,0]**2)
        det = np.multiply(JtJ11,JtJ22)-JtJ12*JtJ12
        JtJi11 = JtJ22/(det+eps)
        JtJi12 = -JtJ12/(det+eps)
        JtJi22 = JtJ11/(det+eps)

        r = np.array(y)
        for k in range(ne):
            r[:,:,:,k]=y[:,:,:,k] - x[:,:,:,0]*np.power(x[:,:,:,1],k)

        proj1 = np.array(r[:,:,:,0])
        for k in range(1,ne):
            proj1 = proj1 + np.multiply(r[:,:,:,k],np.power(x[:,:,:,1],k))
        proj2 = np.array(r[:,:,:,1])
        for k in range(2,ne):
            proj2 = proj2 + k*np.multiply(r[:,:,:,k],x[:,:,:,1]**(k-1))
        proj2 = np.multiply(proj2,x[:,:,:,0])
        
        dx0 = np.array(JtJi11*proj1 + JtJi12*proj2)
        dx1 = np.array(JtJi12*proj1 + JtJi22*proj2)
        x[:,:,:,0] += dx0
        x[:,:,:,1] += dx1    
    x[:,:,:,0]*=yn
    return x   

#simple function to apply the fitting to every 4D matrix in our folder
def FitAll(n_subs, n_visits, ilputPath, outputPath):
    for sub,visits in it.product(range(n_subs), range(n_visits)): #it.product is a very usefule function, enables to apply the algorithm to every n_subs and n_visits combianations
        imputPath=ilputPath + "sub-" + "%03d"%(sub+1) + "_ses-7TV" + "%02d"%(visits+1) + ".nii.gz"
        print(imputPath)
        if os.path.exists (imputPath):
                    ikputPath = "/mnt/senior/rawdata/" + "sub-" + "%03d"%(sub+1) + "/ses-7TV" + "%02d"%(visits+1) + "/anat/"
                    TE = getEchoTimes(10,ikputPath) 
                    eps = 1e-6
                    # input is a 4D matrix, where the 4th dimension are the echoes
                    dat=nib.load(imputPath) 
                    t2starInput = dat.get_fdata()
                    res =  fit_Gauss_Newton_3D(t2starInput, a0=1., b0=1.0, tol=1e-6, num_iter=10)
                    res[:,:,:,1][res[:,:,:,1]<=0] = eps
                    # TE is the vector with the echo times --> TE[1]-TE[0] = dTE 
                    R2star = -np.log(res[:,:,:,1])/(TE[1]-TE[0])
                    affine = dat.affine
                    header = dat.header
                    img=nib.Nifti1Image(R2star,affine,header)
                    img.set_data_dtype(np.float64)
                    nib.save(img,outputPath + 'R2star' + "-sub-" + "%03d"%(sub+1) + "_ses-7TV" + "%02d"%(visits+1) + '.nii.gz')
                    R2star[R2star==0] = eps 
                    T2star = 1000/R2star
                    print(T2star)
                    print(T2star.shape)
                    img=nib.Nifti1Image(T2star,affine,header)
                    img.set_data_dtype(np.float64)
                    nib.save(img,outputPath + 'T2star' + "-sub-" + "%03d"%(sub+1) + "_ses-7TV" + "%02d"%(visits+1) + '.nii.gz')

        else :
            print("Subject or visit don't exist.")
FitAll(185,10,'/home/user7t/aorfanidis/RawdataConcatenated/','/home/user7t/aorfanidis/Résultats_algo_Ludo/')
