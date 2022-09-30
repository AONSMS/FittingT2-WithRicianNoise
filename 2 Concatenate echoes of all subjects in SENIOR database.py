#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 10:10:02 2022
Take all the echos from the rawdata of every subject and compile them in one file per subject per visit.
@author: Orfanidis Alexandre
"""
import nibabel as ni
import glob #glob allow to ignore some parts of a file name
import os
import numpy as np
from nilearn import plotting
from nilearn import image
import matplotlib.pyplot as plt

def concatenateSub(n_sub,n_visits,n_echos,ipputPath):
    for sub in range (n_sub):
        imputPath=ipputPath + 'sub-' + '%03d' % (sub+1) + '/'
        print(imputPath) #follow the algorithm easily
        if os.path.exists(imputPath): #find if a file/folder exist
            
             for visits in range(n_visits):
                 inputPath = imputPath + 'ses-7TV' + "%02d" % (visits+1) + '/anat/'
                 print(inputPath)
                 if os.path.exists(inputPath):
                     data=False #data is false in order to be reshaped
                     #Get the data
                     
                     for echo in range(n_echos):
                         nii = glob.glob(inputPath + '*echo-' +"%02d" % (echo+1) +'_part-mag_MEGRE.nii.gz') #create the path where we'll fetch the file
                         dat=ni.load(nii[0]) #load the path we just input dat=img
                         datOut=dat.get_fdata()  #datOut=Matrix
                         #Reshape the data in 4D (length, width, depth, n_echos)
                         if data is False:
                             sh=datOut.shape+(1,) #add one dimension to datOut
                             data=np.reshape(datOut,sh)
                         else:
                            data=np.concatenate([data, np.reshape(datOut,sh)], axis=3)
                            affine = dat.affine #both affine and header are transformations applied on DatOut.
                            header = dat.header
                        #save the newly acquired data - change the path
                            img = ni.Nifti1Image(data, affine, header) #define image
                            ni.save(img, '/home/user7t/aorfanidis/RawdataConcatenated/'+'sub-' + '%03d' % (sub+185) + '_ses-7TV'+ "%02d" % (visits+1)+".nii.gz") #save the img        
                 else:
                    print('visit does not exist')
            
        else:
            print('Subject does not exist')
        
    
concatenateSub(185,10,10,"/mnt/senior/rawdata/")

