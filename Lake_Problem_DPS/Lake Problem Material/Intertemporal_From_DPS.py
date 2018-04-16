# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:31:31 2018

@author: Rohini

#DPS Formulation 

#Objectives:
#1) Maximize expected economic benefit
#2) Minimize worst case average P concentration 
#3) Maximize average inertia of P control policy
#4) Maximize average reliability 

#Constraints: 
#Reliability has to be <85%

#Decision Variables 
#vars: vector of size 3n 
#n is the number of radial basis functions needed to define the policy
#Each has a weight, center, and radius (w1, c1, r1..wm,cn,rn)

#Time Horizon for Planning, T: 100 Years
#Simulations, N: 100 

#Known Lake Parameters 

#b: 0.42
#q: 2.0

#Initial Lake State: 0

#Economic Benefit Parameters 

#alpha: 0.40
#delta: 0.98

"""

import numpy as np
import os
from sys import *
from math import *
from scipy.optimize import root
import scipy.stats as ss
from borg import *


#Lake Parameters
b = 0.42
q = 2.0

#Natural Inflow Parameters 
mu=0.03
sigma=np.sqrt(10**-5)

#Economic Benefit Parameters

alpha=0.4
delta=0.98

# Set the number of decision variables, objectives and constraints 

nvars = 100
nobjs = 4
nYears = 100
nSamples = 100
nSeeds = 10
nconstrs = 1

#Set Thresholds

reliability_threshold = 0.85
inertia_threshold = -0.02

######################################## Main Lake Problem Model ####################################################

def LakeProblemDPS(seed, *vars): 

    #Solve for the critical phosphorus level
    def pCrit(x):
        return [(x[0]**q)/(1+x[0]**q) - b*x[0]]
      
    sol = root(pCrit, 0.5)
    critical_threshold = sol.x
    
    #Initialize arrays
    average_annual_P = np.zeros([nYears])
    discounted_benefit = np.zeros([nSamples])
    yrs_inertia_met = np.zeros([nSamples])
    yrs_Pcrit_met = np.zeros([nSamples])
    lake_state = np.zeros([nYears+1])
    objs = [0.0]*nobjs
    constrs = [0.0]*nconstrs
    

    #Generate nSamples of nYears of natural phosphorus inflows 
    natFlow = np.zeros([nSamples,nYears])
    for i in range(nSamples):
        np.random.seed(seed+i) 
        natFlow[i,:] = np.random.lognormal(mean=log(mu**2/np.sqrt(mu**2+sigma**2)), sigma=np.sqrt(log((sigma**2+mu**2)/mu**2)), size=nYears)
        
     
    #Run model simulation 
    for s in range(nSamples):
        lake_state[0] = 0
        
        
        for i in range(nYears):
            lake_state[i+1] = lake_state[i]*(1-b) + (lake_state[i]**q)/(1+(lake_state[i]**q)) + vars[i] + natFlow[s,i]
            average_annual_P[i] = average_annual_P[i] + lake_state[i+1]/nSamples
            discounted_benefit[s] = discounted_benefit[s] + alpha*vars[i]*delta**i
            
            if i>=1 and ((vars[i] - vars[i-1]) > inertia_threshold):
                yrs_inertia_met[s] = yrs_inertia_met[s] + 1
                
            if lake_state[i+1] < critical_threshold:
                yrs_Pcrit_met[s] = yrs_Pcrit_met[s] + 1

        
                
 # Calculate minimization objectives (defined in comments at beginning of file)
    objs[0] = -1*np.mean(discounted_benefit) #average economic benefit
    objs[1] = np.max(average_annual_P) #minimize the max average annual P concentration
    objs[2] = -1*np.sum(yrs_inertia_met)/((nYears-1)*nSamples) #average percent of transitions meeting inertia thershold
    objs[3] = -1*np.sum(yrs_Pcrit_met)/(nYears*nSamples) #average reliability
    
    constrs[0] = max(0.0, reliability_threshold - (-1*objs[3]))

    return (objs,constrs)               
     
#########################################################################################################################
 
#Set up Borg Interface 

for j in range(nSeeds):
    borg = Borg(nvars, nobjs, nconstrs, LakeProblemDPS, j)
    borg.setBounds(*[[0.01, 0.1]]*(nvars))
    borg.setEpsilons(0.01, 0.01, 0.0001, 0.0001)

    result = borg.solve({"maxEvaluations":50000})
    
    f = open(os.getcwd() + '/setsI/' + \
        str(j+1) + '.set','w')

    f.write('#Borg Optimization Results\n')
    f.write('#First ' + str(nvars) + ' are the decision variables, ' \
        'last ' + str(nobjs) + ' are the objective values\n')

    for solution in result:
        line = ''
        for i in range(len(solution.getVariables())):
            line = line + (str(solution.getVariables()[i])) + ' '

        for i in range(len(solution.getObjectives())):
            line = line + (str(solution.getObjectives()[i])) + ' '

        f.write(line[0:-1]+'\n')

    f.write("#")
        
    f.close()
      
 #####################################################################################################################       
        
        