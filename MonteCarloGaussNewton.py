#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 14:48:21 2022
MonteCarlo simulation for the GaussNewton algorithm 
@author: Orfanidis Alexandre, Ludovic de Rochefort made the Gauss Newton part
"""
import numpy as np
import glob
import json
import matplotlib.pyplot as plt 

def fit_Gauss_Newton_3D(y, a0=1.0, b0=1.0, tol=1e-6, num_iter=10):
    # fittinf a 3D + t volume to a*b^k
    # then, take R2* = -log(b)/dTE
    s=list(y.shape) #dim of the matrix y (in our case : (200,100,10))


    ne=s[2]# number of echos

    s[2]=2
    x=np.ones(s,dtype=float)

    x[:,:,0]=a0# initialization of a, could be done differently (such as the sum_of-square signal), never initialize at 0
    x[:,:,1]=b0# initialization of b, never initialize at 0
    eps=1e-6
    # there should be a tolerance as well for early stopping
    # can be implemented
    # normalize the input 
    yn=np.sum(np.absolute(y)**2,axis=2)**0.5 + eps
    

    y /=  yn[:,:,None]

    for itr in range(num_iter):
        # jacobian and inverse
        JtJ11 = np.array(1.0 + np.power(x[:,:,1],2.0))
        JtJ12 = np.array(x[:,:,1])
        JtJ22 = np.ones(JtJ12.shape)

        for k in range(2,ne):
            JtJ11 += np.power(x[:,:,1],2*k)
            JtJ12 += k*np.power(x[:,:,1],2*k-1)
            JtJ22 += (k*k)*np.power(x[:,:,1],2*k-2)
        
        JtJ12 = np.multiply(JtJ12,x[:,:,0])
        JtJ22 = np.multiply(JtJ22,x[:,:,0]**2)
        det = np.multiply(JtJ11,JtJ22)-JtJ12*JtJ12
        JtJi11 = JtJ22/(det+eps)
        JtJi12 = -JtJ12/(det+eps)
        JtJi22 = JtJ11/(det+eps)

        r = np.array(y)
        for k in range(ne):
            r[:,:,k]=y[:,:,k] - x[:,:,0]*np.power(x[:,:,1],k)

        proj1 = np.array(r[:,:,0])
        for k in range(1,ne):
            proj1 = proj1 + np.multiply(r[:,:,k],np.power(x[:,:,1],k))
        proj2 = np.array(r[:,:,1])
        for k in range(2,ne):
            proj2 = proj2 + k*np.multiply(r[:,:,k],x[:,:,1]**(k-1))
        proj2 = np.multiply(proj2,x[:,:,0])
        
        dx0 = np.array(JtJi11*proj1 + JtJi12*proj2)
        dx1 = np.array(JtJi12*proj1 + JtJi22*proj2)
        x[:,:,0] += dx0
        x[:,:,1] += dx1    
    x[:,:,0]*=yn
    
    return x   

FinalMT2=[] #(n_sub=200,n_iter,n_echos)
sigma = 1
mu=0
S0=30      #Base signal
#S0=500      #Base signal
T2_star_list=list(range(1,201))
TE = (0.005,0.01,0.015,0.020,0.025,0.030)

for a in range(len(T2_star_list)): #Apply all the algorithm for different T2 value

    theoS=[]#theoS will be the list where we'll put all the datas we need for our figure
    for b in range(len(TE)): #loop to obtain the associated signal depending of the T2 in this iteration      
        S=S0*np.exp(-TE[b]*1000/T2_star_list[a]) #ms #equation giving the signal value
        theoS.append(S) # add the result to theoS  
       
    n_iter=100000 #number of iteration we want (objectif=100000)
    count=0
    List_Meas_T2star=[] # (n_iter,n_echos)
    while count<n_iter: #loop to create our of 100000 noised datas
      theoSandNoise=[] #(10=n_echos)
      for c in range(len(theoS)):  
          noise = np.random.normal(mu, sigma, 1) + 1j*np.random.normal(mu, sigma, 1) #our way of generating rician noise
          theoSandNoise.extend(np.abs(theoS[c] + noise))    #add noise to our theorical values
          #theoSandNoise.extend(theoS[c] + np.linalg.norm(np.random.normal(scale=s, size=(N, 2)) + [[v,0]], axis=1))    #add noise to our theorical values
      List_Meas_T2star.append(theoSandNoise)
      theoSandNoise=[]#reset our list for the next iteration
      count+=1
        
    FinalMT2.append(List_Meas_T2star)

t2starInput = np.asarray(FinalMT2)
eps = 1e-6
res =  fit_Gauss_Newton_3D(t2starInput, a0=1., b0=1.0, tol=1e-6, num_iter=10)
res[:,:,1][res[:,:,1]<=0] = eps
# TE is the vector with the echo times --> TE[1]-TE[0] = dTE 
MR2star = -np.log(res[:,:,1])/(TE[1]-TE[0])
MT2star = 1000/MR2star

MT2starMean=[]
for i in range(len(MT2star)):
   MT2starMean.append(np.average(MT2star[i]))#to obtain the mean of our measured T2*, thanks to which we will see the result of our algorithm compared to the theorical values
MT2starMeanStD=np.std(MT2star,1)#Calculate standard deviation of MT2star (orange on the plot)
MT2starMean=np.asarray(MT2starMean)
#Create the abscisse list of same size than MT2star to create the scatter plot, contain T2* Theorical values :
abscisse=[]
for a in range(len(T2_star_list)):
    Préabsci = [a] * n_iter
    #print(Préabsci)
    abscisse.append(Préabsci)
  
#show MonteCarlo results for our GaussNewton algorithm:
plt.figure()
fig, ax = plt.subplots()
plt.xlabel('T2* Théorique')
plt.xlim([0, 200]) #set the limit in the physiological values
plt.ylabel('T2*mesuré')
plt.ylim([0, 200]) #set the limit in the physiological values
plt.title('T2* mesuré et théorique'+' : '+str(n_iter)+' itérations')
plt.scatter(abscisse,MT2star) #points obtained during Monte Carlo
plt.errorbar(T2_star_list,MT2starMean,MT2starMeanStD,c='orange') #measured T2* curve 
plt.plot(T2_star_list,T2_star_list,'-', c='red') #theorical curve
plt.show
plt.savefig(str(n_iter)+'MonteCarloResultsGaussNewton')

#Give the figure for relative uncertainty
plt.figure()
fig, ax = plt.subplots()
plt.xlabel('T2* théorique (ms)')
plt.xlim([0, 200]) #set the limit in the physiological values
plt.ylabel('DeltaT2* mesuré / T2* mesuré (%)')
plt.ylim([0, 40]) #limit arbitrarly set to give a readable image
plt.title("Importance de l'incertitude en fonction du T2* théorique"+' : '+str(n_iter)+' itérations')
plt.plot(T2_star_list,100*MT2starMeanStD/MT2starMean,'-', c='red') 
plt.show
plt.savefig(str(n_iter)+'IncertitudeRelativeGaussNewton')

bias=np.polyfit(T2_star_list,MT2starMean,1)
print ('BIAS : '+ str(bias))
