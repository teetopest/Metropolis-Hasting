#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 15:00:00 2023

@author: toheebibrahim
"""

import numpy as np
import random
from matplotlib import pyplot as plt 




# length of row of lattice also equal to the length of column 
L = 10 # int(input('enter the number of row or column of lattice')) 

p = 3 #int(input('enter the initial configuration 1 for +ve, -1 for -ve, other number for random '))

# Ferromagnetic constant chosen to be equal to the Boltzmann constant
J = kb = 1 

# Temperature
T = 0.4


beta_l = [0.4, 0.5]

N = L**2

N_sw = [N, 2*N, 3*N, 4*N, 5*N, 6*N, 7*N, 8*N, 9*N] #[2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 180000 ]



max_sw = max(N_sw)

Time = np.linspace(0,N,N+1)

n = 10000


def gen_ones(L):
    return np.ones((L, L))

def gen_minus_ones(L):
    return np.full((L, L), -1)

def gen_rand(L):
    return np.random.choice([-1, 1], size=(L, L), p=[0.5, 0.5])


# Choice of initial configuration      
def init_config(config):
    if config == 1:
        return gen_ones(L)
    elif config == -1:
        return gen_minus_ones(L)
    else:
        return gen_rand(L)

# Neigbours
def nhb(M,i,j):
    res = 0
    if i+1 >= L and j+1 >= L:
        res += M[i+1-L][j]+ M[i-1][j] +  M[i][j+1-L] + M[i][j-1]
    elif i+1 >= L and j < L:
        res += M[i + 1 - L][j] + M[i - 1][j] +  M[i][j+1] + M[i][j-1]
    elif i < L and j + 1 >= L:
        res += M[i+1][j]+ M[i-1][j] +  M[i][j+1-L] + M[i][j-1]
    else:
        res += M[i+1][j]+ M[i-1][j] +  M[i][j+1] + M[i][j-1]
    return res
  



# Hamiltonian
def hamil(M,l):
    sum_up = 0
    for i in range(l):
        for j in range(l):
            sum_up += M[i][j]*nhb(M,i,j)
    sum_up = sum_up * -(J/2) 
    #print(sum)
    return sum_up  


'''
# Spin configuration: multiple times spin
def spin_n(config,l):
    count = 1
    config_flip = config
    mag_l = [(np.sum(config_flip) * 1/N)] 
    energy_l = [hamil(config_flip, l) * 1/N]
    for beta in beta_l:
        #plt.subplots(3,3)
        for n in range(0, max_sw+1, 1):
            rand_row = np.random.randint(0,l)
            rand_col = np.random.randint(0,l)
            config_flip[rand_row, rand_col] = -config_flip[rand_row][rand_col]
            ham_config_flip = hamil(config_flip,l)
            ham_config = hamil(config, l)
            res = ham_config_flip - ham_config
            r = min(1, np.exp(-beta * res ))    #Computing Prob. of acc: Metropolis rule
            if random.random() <= r:            
                config_flip = config_flip
            else:
                config_flip = config
            
        
            if n <= N-1:
                energy_1 = hamil(config_flip,l) * 1/N # energy
                energy_l.append(energy_1)
                mag = np.sum(config_flip) * 1/N # magnetisation
                mag_l.append(mag)
            if n== N-1:
                plt.plot(Time, mag_l)
                plt.ylabel('Magnetization')
                plt.xlabel('Time')
                plt.title("Magnetization vs Time for beta = {}".format(beta))
                plt.show()
                plt.plot(Time, energy_l)
                plt.title("Energy vs Time for beta = {}".format(beta))
                plt.ylabel('Energy')
                plt.xlabel('Time')
                plt.show()
            
            if (n % N == 0) and (n!=0):
                plt.subplot(3, 3, count)
                plt.imshow(config_flip, cmap= 'binary')
                count = count + 1    
                #plt.show()
        plt.show()        
        count = 1
        mag_l = [np.sum(config)]
        energy_l = [hamil(config, l) * 1/N]
        
    return
           
'''

def mean(config,l):
    config_flip = config
    
    mag_l = [(np.sum(config_flip) * 1/N) ] 
    
    energy_l = [hamil(config_flip, l) * 1/N]
    
    mag_sum = ((np.sum(config_flip) * 1/N))/n
    
    energy_sum = (hamil(config_flip, l) * 1/N)/n
    
    mag_sum_list = []
    
    energy_sum_list = []
    
    mag2_sum = (((np.sum(config_flip) * 1/N))**2)/n
    
    energy2_sum = (((np.sum(config_flip) * 1/N))**2)/n
    
    Suscp_list = []
    
    Spec_heat_list = []
    
    for beta in beta_l:
        for i in range(0, n, 1): 
            rand_row = np.random.randint(0,l)
            rand_col = np.random.randint(0,l)
            config_flip[rand_row, rand_col] = -config_flip[rand_row][rand_col]
            ham_config_flip = hamil(config_flip,l)
            ham_config = hamil(config, l)
            res = ham_config_flip - ham_config
            r = min(1, np.exp(-beta * res ))    #Computing Prob. of acc: Metropolis rule
            if random.random() <= r:            
                config_flip = config_flip
            else:
                config_flip = config
            
            # energy
            energy_1 = hamil(config_flip,l) * 1/N 
            energy2 = energy_1 ** 2 
            energy_l.append(energy_1)
            
            # magnetisation
            mag = np.sum(config_flip) * 1/N 
            mag2 = mag**2
            mag_l.append(mag)
            
            # mean magnetisation 
            mag_sum = mag_sum + mag/n       # mean magnetisation
            mag2_sum = mag2_sum + (mag2 / n) # mean of square magnetisation
            
            # mean energy
            energy_sum = energy_sum + energy_1/n # mean energy
            energy2_sum = energy2_sum + energy2/n # mean of square energy
            
            
        # list of mean magnetisation for each beta   
        mag_sum_list.append(mag_sum)
        
        # list of mean energy for each beta
        energy_sum_list.append(energy_sum)
        
        # Susceptibility
        susc1 = beta  * N * ( mag2_sum - (mag_sum**2) )
        Suscp_list.append(susc1)
        
        # Specific heat
        spec_heat1 = beta ** 2 * N * (energy2_sum - (energy_sum ** 2))
        Spec_heat_list.append(spec_heat1)
        
    plt.plot(beta_l,energy_sum_list)
    plt.title("Mean Energy vs Temperature for varying")
    plt.ylabel('Mean Energy')
    plt.xlabel('Temp')
    plt.show()
    
    plt.plot(beta_l, mag_sum_list)
    plt.title("Mean magnetisation vs Time for varying beta")
    plt.ylabel('Mean magnetization')
    plt.xlabel('Temp')
    plt.show()
    
    plt.plot(beta_l, Suscp_list)
    plt.title("Susceptibility vs temperature")
    plt.ylabel('Susceptibility')
    plt.xlabel('Temp')
    plt.show()
    
    plt.plot(beta_l, Spec_heat_list)
    plt.title("Specific heat vs temperature")
    plt.ylabel('Specific heat')
    plt.xlabel('Temp')
    plt.show()
    
        
        
            











def metropolis_hasting(p):
    init = init_config(p)
    mean(init,L)
    #spin_n(init,L)
     #print(init)
     #spin2 = spin(init)
     #print(init[0][1])
     #print(L)
     #d = delta(init, L)
     #Pacc(init, L)
     #mag_den(init, L)
     #plt.plot(Time,den(init, L).mag_l) 
     #plt.show()
     #den(init, L)
     
     
     #hamil(init, L)

    
metropolis_hasting(p) 
                 
        
    
