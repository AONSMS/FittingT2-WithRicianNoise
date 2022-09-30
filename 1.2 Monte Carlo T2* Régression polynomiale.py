#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:11:27 2022

@author: user7t
"""

import numpy as np
import glob
import json
import matplotlib.pyplot as plt 

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
MT2starmeanMean=[]
MT2starMeanStDMean=[]
TT2_star_list=[]
COUNT=0
while COUNT<500:
    FinalMT2=[]
    TrueMeasuredSignal=[]
    v = 4
    s = 9
    N = 1# the standard deviation of the normal distribution
    S0=500      #Base signal
    T2_star_list=list(range(1,11))
    TT2_star_list.append(T2_star_list)
    TE = getEchoTimes(10,"/mnt/senior/rawdata/sub-007/ses-7TV08/anat/") #We get a list with all the TE we need
    listofTE=TE*10
    

    for a in range(len(T2_star_list)): #Apply all the algorithm for different T2 value
            
            S0_list=[]
            theoS=[]#theoS will be the list where we'll put all the datas we need for our figure
            for b in range(len(TE)):       
                S=S0*np.exp(-TE[b]*1000/T2_star_list[ a ]) #ms #equation giving the signal value
                S0_list.append(S0)#create a list with a n_S0 equel to n_TE
                theoS.append(S) # add the result to theoS  
                # plt.figure(a) #create a figure for each T2_star
            count=0
                
            List_theoSandNoise=[]
            theoSandNoise=[]
            n_iter=10
            while count<n_iter: #number of iteration we want (objectif=100000)
                for c in range(len(theoS)):    
                    theoSandNoise.extend(theoS[c] + np.linalg.norm(np.random.normal(scale=s, size=(N, 2)) + [[v,0]], axis=1))    #add noise to our theorical values
                    #print(theoSandNoise)
                        

                w=np.asarray(theoSandNoise)     
                List_theoSandNoise.append(theoSandNoise)
                theoSandNoise=[]
                count+=1
                
            MeasuredSignal=np.polyfit(TE,List_theoSandNoise,deg=1)
            MeasuredSignal=MeasuredSignal.tolist()
                
            TrueMeasuredSignal.append(MeasuredSignal[1])
            FinalMT2.append(List_theoSandNoise)
                
                #Create the abscisse list of same size than MT2star to create the scatter plot, contain T2* Theorical values :
    abscisse=[]
    for a in range(len(T2_star_list)):
        Préabsci = [a] * n_iter
        #print(Préabsci)
        abscisse.append(Préabsci)
            
    TrueMeasuredSignal=np.asarray(TrueMeasuredSignal)
    TE = [item *(-1000) for item in TE]
    theoS=np.asarray(theoS)
    MT2star=TE/np.log(TrueMeasuredSignal/500)
    #print(MT2star)
    
    t2starInput = np.asarray(FinalMT2)
    MT2star2=[]
    MT2starMean=[]
    for i in range(len(MT2star)):
                    x=3*i+20
                    MT2star2.append(MT2star[i].clip(0,x))
                    #MT2star2.append(MT2star[i][MT2star[i]>=0])
                    MT2starMean.append(np.average(MT2star2[i]))
                    MT2starMeanStD=np.std(MT2star2,1)#Calculate standard deviation of MT2star (orange on the plot)
    MT2starmeanMean.append(MT2starMean)     
    MT2starMeanStDMean.append(MT2starMeanStD)#doing a mean of that is problably non sense          
         
            #create the figure:
    #plt.figure()
    #fig, ax = plt.subplots()
   
    plt.xlabel('T2* Théorique')
    plt.xlim([0, 10])
    plt.ylabel('T2*mesuré')
    plt.ylim([0, 10])
    plt.title('T2* mesuré et théorique'+' : '+'100000'+' itérations')
    plt.scatter(abscisse,MT2star,c='blue')
    #plt.errorbar(T2_star_list,MT2starMean,MT2starMeanStD,c='orange')
    #plt.plot(T2_star_list,MT2starMean,c='orange')
    plt.plot(T2_star_list,T2_star_list,'-', c='red')
    #plt.show
    
    
    COUNT+=1

MT2starmeanMean=np.mean(MT2starmeanMean, axis=0)
MT2starmeanMean=np.asarray(MT2starmeanMean)
MT2starMeanStDMean=np.mean(MT2starMeanStDMean, axis=0)
MT2starMeanStDMean=np.asarray(MT2starMeanStDMean)
#MT2starmeanMean=np.average(MT2starmeanMean)

plt.errorbar(T2_star_list,MT2starmeanMean,MT2starMeanStDMean,c='orange')
plt.plot(T2_star_list,MT2starmeanMean,c='red')
plt.show
#plt.savefig(str(n_iter*COUNT) + 'T2linearfittingGaussianNoiseTriche'+ '[0,10]')
"""

plt.figure()
fig, ax = plt.subplots()
plt.xlabel('T2* théorique (ms)')
plt.xlim([0, 200])
plt.ylabel('DeltaT2* mesuré / T2* mesuré (%)')
plt.ylim([0, 100])
plt.title("Importance de l'incertitude en fonction du T2* théorique"+' : '+str((n_iter*COUNT))+' itérations')
plt.plot(T2_star_list,100*MT2starMeanStDMean/MT2starmeanMean,'-', c='red')
plt.show
plt.savefig(str((n_iter*COUNT))+'ImportanceIncertitudeT2ricienfitLinéaireTriche' + '2')
"""