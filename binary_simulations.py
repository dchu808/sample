# !/usr/bin/env python

import numpy as np
import pandas as pd
import pylab as plt
import matplotlib
from astropy import units as u
from astropy.stats import LombScargle
from astropy.time import Time
from astropy.table import Table
import os,scipy, subprocess
from scipy import stats
from scipy.stats import f
import shutil, logging, datetime
import MySQLdb as mdb
from tqdm import tqdm
from math import pi, cos, sin, fabs,sqrt, floor

font = {        'size'   : 20}
matplotlib.rc('font', **font)

'''
This code is a collection of simulations for binary star properties.

'''

'''
Below is the equation for the amplitude, defined as K = 2 pi a sin i / (P(1-e**2)^(1/2))

'''
def k_equation(a, i, per, e):
    ##a is in AU, period is in days, i should be radians
    k = 2* np.pi* a * np.sin(i) / (per * np.sqrt(1 - e**2))
    return k

    '''
    Below is the equation for obtaining a semi-major axis (a) from a period.
    Using Kepler's Laws, use relation of log a \propto 3/2 log p

    '''

def loga_from_logperiod(log_period):
    log_a = (3./2)*log_period
    return log_a

    '''
    Below is the equation minimum separation
    a is semi-major axis, e is eccentricity

    '''
    
def rmin(a,e):
    r = a*(1-e)
    return r
    
    '''
    The following is code for simulating values of K.
    K = 2 pi a sin i / (P(1-e**2)^(1/2))
    Draw from distribution of values, then see what is distribution of K values
    power values came from Sana+ 2012, see supplimental section B.3
    log period pi = -0.55
    eccentricity eta = -0.45
    kappa = -0.55
    '''
    
def kvalue_simulation(ntrials=100000,popsize = 100000,pi_exponent = -0.55, eta_exponent = -0.45, min_separation = 0.004, make_plot=True,load_sim = False):

    ##create distribution of parameters
    ##inclination distribution, assuming it's flat
    i_dist = np.random.uniform(low=0.,high=np.pi,size=popsize)
    ##check i_distribution
    # plt.figure()
    # plt.hist(i_dist)
    # plt.xlabel('i (radians)')
    # plt.show()
    ##period distribution from power law
    ##maximum power value is 1000 days
    logper_dist = np.log10(365) * np.random.power(pi_exponent + 1,popsize)
    ##check peirod_distribution
    # plt.figure()
    # plt.hist(logper_dist)
    # plt.xlabel('log P')
    # plt.ylabel('Number')
    # plt.show()
    ##eccentricity distribution
    ecc_dist = np.random.power(eta_exponent + 1,popsize)
    ##check e_distribution
    # plt.figure()
    # plt.hist(ecc_dist)
    # plt.xlabel('e')
    # plt.ylabel('Number')
    # plt.show()
    ##test
    i_random = np.random.choice(i_dist,replace=True)
    logper_random = np.random.choice(logper_dist,replace=True)
    e_random = np.random.choice(ecc_dist,replace=True)
    # loga_random = np.random.choice(logper_dist,replace=True)
    loga_random = (3/2)*logper_random
    per_random = 10**logper_random
    a_random = 10**loga_random
    print(i_random)
    print(per_random)
    print(a_random)
    print(e_random)
    print(k_equation(a_random, i_random, per_random, e_random))
    ##draw from the distributions to get K value
    k_values = np.empty(ntrials)
    i_values = np.empty(ntrials)
    e_values = np.empty(ntrials)
    a_values = np.empty(ntrials)
    logp_values = np.empty(ntrials)
    if load_sim == True:
        kval_sort = np.load('simulations/kvalue_simulation.npy')
        y_array = np.arange(kval_sort.size)
        s = float(kval_sort.size)
        #this way the y-axis goes from 0 - 1.
        y_array_norm = y_array/s
        ##looking for 95% confidenc level
        print('95percent CL value')
        idx95cl = np.where(y_array_norm > 0.95)[0]
        print(kval_sort[idx95cl[0]])
    else:
        i = 0
        while i < ntrials:
            i_random = np.random.choice(i_dist,replace=True)
            logper_random = np.random.choice(logper_dist,replace=True)
            # loga_random = np.random.choice(logper_dist,replace=True)
            loga_random = loga_from_logperiod(logper_random)
            per_random = 10**logper_random
            a_random = 10**loga_random
            e_random = np.random.choice(ecc_dist,replace=True)
            ##check to see if the binaries are separated, otherwisebinaries will collide
            r_min = rmin(a_random,e_random)
            if r_min > min_separation:
                k_value_it = k_equation(a_random, i_random, per_random,e_random)
                k_values[i] = k_value_it
                i_values[i] = i_random
                e_values[i] = e_random
                logp_values[i] = logper_random
                a_values[i] = a_random
                i+=1
            else:
                continue
        ##find the 95% K value
        print(len(k_values))
        kval_sort = np.sort(k_values)
        y_array = np.arange(kval_sort.size)
        s = float(kval_sort.size)
        #this way the y-axis goes from 0 - 1.
        y_array_norm = y_array/s
        ##looking for 95% confidenc level
        print('95percent CL value')
        # idx95cl = np.where(y_array_norm > 0.95)[0]
        # print(kval_sort[idx95cl[0]])
        print(np.percentile(k_values,95))
        np.save('simulations/kvalue_simulation',k_values)
    if make_plot == True:
        plt.figure()
        plt.hist(k_values)
        plt.xlabel('K (km/s)')
        plt.savefig('plots/kvalue_simulation_hist.pdf')
        plt.close()
        ##make cdf plot
        plt.figure()
        plt.plot(kval_sort,y_array_norm)
        plt.xlabel('K (km/s)')
        plt.ylabel('CDF')
        plt.savefig('plots/kvalue_simulation_cdf.pdf')
        plt.close()
        plt.figure()
        plt.hist(e_values)
        e_range = np.linspace(0,1,25)
        plt.xlabel('Eccentricity')
        plt.show()
        plt.figure()
        plt.hist(logp_values)
        plt.xlabel('Log P')
        plt.show()
        plt.figure()
        plt.hist(a_values)
        plt.xlabel('s')
        plt.show()

        '''
        The following is code for simulating values of K for given periods.
        K = 2 pi a sin i / (P(1-e**2)^(1/2))
        Draw from distribution of values, then see what is distribution of K values
        power values came from Sana+ 2012, see supplimental section B.3
        eccentricity eta = -0.45
        '''

def k_simulation_per_period(ntrials=1000,popsize = 10000,eta_exponent = -0.45, min_separation = 0.004,number_of_kbins = 50, per_start=2.,per_end=200.,step=0.05):
    ##determining the periods to run for
    n_loop=int(floor((per_end-per_start)/step) + 1)
    per_loop = np.linspace(per_start,per_end,n_loop)
    bin_edges = np.linspace(0,100,num=number_of_kbins,endpoint=True)
    ##population distributions
    i_dist = np.random.uniform(low=0.,high=np.pi,size=popsize)
    ecc_dist = np.random.power(eta_exponent + 1,popsize)
    ##full histogram to run 2d hist on
    full_hist = np.empty((n_loop,ntrials))
    full_period_array = np.empty((n_loop,ntrials))
    ##loop through the periods
    for i in tqdm(range(n_loop)):
        period_for_this_loop = per_loop[i]
        per_exp_for_this_loop = np.log10(period_for_this_loop)
        loga = loga_from_logperiod(per_exp_for_this_loop)
        a = 10**loga
        k_array_temp = np.empty(ntrials)
        period_array_temp = np.empty(ntrials)
        for j in range(ntrials):
            ##distributions to draw from
            i_random = np.random.choice(i_dist,replace=True)
            e_random = np.random.choice(ecc_dist,replace=True)
            k_value_it = k_equation(a, i_random, period_for_this_loop,e_random)
            k_array_temp[j] = k_value_it
            period_array_temp[j] = period_for_this_loop
        full_hist[i] = k_array_temp
        full_period_array[i] = period_array_temp
    ##make a 2d histogram
    ##flatten the histogram and period array
    hist_flat = full_hist.flatten()
    per_flat = full_period_array.flatten()
    print(len(hist_flat))
    print(len(per_flat))
    plt.figure(figsize=(12,12))
    plt.hist2d(per_flat,hist_flat,range=[[2,200],[0,150]],normed=True,bins=[200,number_of_kbins],cmap='viridis')
    # plt.imshow(full_hist,cmap='viridis',aspect='auto')
    plt.xlabel('Period (Days)')
    plt.ylabel('K Amplitude (km/s)')
    # per_ticks = np.linspace(0,full_hist.shape[1],10)
    # per_ticks_labels = np.linspace(per_start,per_end,10)
    # plt.xticks(per_ticks,(per_ticks_labels))
    # k_ticks = np.linspace(0,number_of_kbins,5)
    # k_ticks_labels = np.linspace(100,0,5)
    # plt.yticks(k_ticks,(k_ticks_labels))
    plt.colorbar()
    plt.savefig('/u/devinchu/Research/s-star_rv_periodicity/plots/simtest1.pdf')
    plt.close()

    '''
    The following is code for simulating values of K for given periods from chains.
    K = 2 pi a sin i / (P(1-e**2)^(1/2))
    Draw from distribution of values, then see what is distribution of K values
    power values came from Sana+ 2012, see supplimental section B.3
    eccentricity eta = -0.45
    '''
            
def plot_period_k_2d(star_name = 'S0-1_only_orbit_mod_keck',number_of_kbins = 50,overplot_95=False):
    chain_dir = '/localcompute_devin/s-star_periodicity_chains/'+star_name+'/'
    # chain_dir = '/Volumes/Files/s-star_periodicity_chains/'+star_name+'/'
    ##get the number of chain directories
    # number_of_chain_dirs = call(['ls',chain_dir,'| wc -l'])
    # number_of_chain_dirs = call('ls '+chain_dir+' | wc -l',shell=True)
    number_of_chain_dirs = check_output('ls '+chain_dir+' | wc -l',shell=True)
    # number_of_chain_dirs = os.system('ls '+chain_dir+' | wc -l')
    print(number_of_chain_dirs)
    cmd = Popen('ls '+chain_dir+' | wc -l',shell=True,stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output,err = cmd.communicate()
    print(output)
    ##temporary fix, hopefully can generalize this
    per_start=2.
    per_end=200.
    step=0.05
    n_loop=int(floor((per_end-per_start)/step) + 1)
    per_loop = np.linspace(per_start,per_end,n_loop)
    full_hist = np.empty(0)
    full_weights = np.empty(0)
    full_per_array = np.empty(0)
    # bin_edges = np.linspace(0,100,num=number_of_kbins,endpoint=True)
    for i in range(n_loop):
        ch = np.loadtxt(chain_dir+'chains_'+str(i)+'/.txt')
        full_hist = np.append(full_hist,ch[:,2])
        full_weights = np.append(full_weights,ch[:,0])
        period_of_this_loop = per_loop[i]
        per_temp = np.full_like(ch[:,2],period_of_this_loop)
        full_per_array = np.append(full_per_array,per_temp)
    plt.figure(figsize=(12,12))
    plt.hist2d(full_per_array,full_hist,weights=full_weights,bins=[200,number_of_kbins],cmap='viridis')
    plt.xlabel('Period (Days)')
    plt.ylabel('K Amplitude (km/s)')
    plt.colorbar()
    if overplot_95 == True:
        ##load 95 limit
        result95 = np.loadtxt('/u/devinchu/Research/s-star_rv_periodicity/multinest_results/'+star_name+'_maxPer')
        plt.plot(result95[:,0],result95[:,1],color='C1')
    plt.savefig('/u/devinchu/Research/s-star_rv_periodicity/plots/'+star_name+'_khist_2d.pdf')
    plt.close()

    '''
    For some stars, the mass needs to come from an isochrone generated from evolutionary model
    Read the isochrone file, get a mass
    '''
def get_mass_from_isochrone(iso_file):
    iso = Table.read(iso_file)
    iso.info
    
        
    