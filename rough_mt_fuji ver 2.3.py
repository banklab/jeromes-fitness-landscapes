#%%
"""
Rough Mt. Fuji-Type fitness landscapes as described in Aita T, Hidefumi U, 
inaoka T, Nakajima M, Kokubo T and Husimi Y (2000) Analysis of a Local Fitness Landscape with a Model of
the Rough Mt. Fuji-Type Landscape: Application to Prolyl Endopeptidase and Thermolysin.  
Biopolymers 54: 64–79

Jérôme Stäheli 17-111-220

version 2.3
30.11.2021
"""

#%% Packages
"""
Load Packages (matplotlib.pyplot for plots, Numpy for normal distribution, random for random numbers 
and intertools for all possible gene combinations)
"""
import random as rd
from itertools import product
import numpy as np
import matplotlib.pyplot as plt


#%% Fitness influence of locus

# Creates fitness influence value, these are taken from a normal distribution
def fit_infl(site_num, stand_dev_fit, dist_mean_fit):
    
    return np.random.normal(dist_mean_fit, stand_dev_fit, site_num)


#%% Optimal genome
"""
generates random numbers for every site to generate the optimal genome 
with the help of a defined number of sites and a number of possible variations
"""
def opt_gen(site_num, max_res_var):
    
    # Empty list to save every position of the optimal gene
    opt_gen = []
    
    # generates the optimal genotype randomly with the help of the site number and the number of possible 
    # replacements
    for num in range(0, site_num):
        opt_gen.append(rd.randint(1, max_res_var))
    
    print(opt_gen, "optimal gene")
    return opt_gen


#%% Fitness calculation
"""
Calculates the fitness of a given genotype and generates the fitness with the help of
a fitness influence list, the optimal genotype and its fitness
"""
def add_fit_calc(fit_infl, opt_fit, opt_gen, cur_gen):
    
    # Fitness of the gene is calculated from the fitness of the optimal genotype
    cur_fit = opt_fit
    
    # Checks every site if it is the same as the optimal genotype and if it is not the fitness influence is 
    # substracted
    for num in range(0, site_num):
        
        if cur_gen[num] != opt_gen[num]:
            cur_fit -= float(fit_infl[num])
        
    print(cur_fit)
    return cur_fit


#%% Rough Mt. Fuji-Type fitness landscape
"""
Generates the Rough Mt Fuji-Type fintnes landscape with the fitness influence of every site, the number
of sites, how many diffrent things can be expressed on one site, the mean and the standart deviation 
of the normal distribution for non additive effects, the optimal genotype and its fittnes
"""
def rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen):
    
    # generates all possible combination with the help of site number and the number of possibilities of 
    # each site
    pos_comb = list(product({1, max_res_var},repeat = site_num))
    
    # generates the fitness of all possible genotypes and saves them in a list 
    #(with the help of the fitness calculation)
    fit_list = []
    
    for num2 in range(0, len(pos_comb)):
        
        fit_list.append(float(np.random.normal(dist_mean, stand_dev, 1) + 
                              add_fit_calc(fit_infl, opt_fit, opt_gen, pos_comb[num2])))

    
    # generates a list with the number of diffrent sites compared to a all non mutated genotype
    # this is just used for the creation of a plot
    mut_list = []

    for num3 in range(0, len(pos_comb)):
        
        num_mut  = 0
        
        for num4 in range(0, len(pos_comb[num3])):
            
            if pos_comb[num3][num4] != 1:
                num_mut += 1
    
        mut_list.append(num_mut)
        
    print(pos_comb)

    
    # adds lines in the plot that connect all the dots of each genotype
    for num4 in range(0, len(pos_comb)):
        
        for num5 in range(0, len(pos_comb)):
            
            sim = site_num
            
            for num6 in range(0, site_num):
            
                if pos_comb[num4][num6] == pos_comb[num5][num6]:
                    sim -= 1
            
            if sim == 1:
                plt.plot([mut_list[num4], mut_list[num5]],
                         [fit_list[num4], fit_list[num5]],
                         'k-')

    
    # plots the genotypes with the number of mutations and the lines that indicates to which genotype they
    # can further mutate
    plt.plot(mut_list, fit_list,  'ro')
    print(fit_list, "list with final fitnesses")
    print(mut_list, "list with numbers of mutation")
    
    # retun of the list of possible combinations and the corresponding fitnesses
    return fit_list
    return pos_comb
    

#%% Parameters

# Number of sites that have influence on fitness and can mutate
site_num        = 5
# Number of possible aminoacides on each site
max_res_var     = 2
# Fitness of the optimal Genotype
opt_fit         = 2
# Standart deviation of non additive effects
stand_dev       = 0.1
# Distribution mean of non additive effects
dist_mean       = 0
# Standart distribution for the fitness influence of the site
stand_dev_fit   = 0.1
# Distribution mean of the fitness influence of the site
dist_mean_fit   = 0


#%% Use funtions

# fitness influence (fit_infl) and optimal genotype (opt_gen) are needed for the Rough Mt Fuji function
# might later on include it in the funtion
fit_infl = fit_infl(site_num, stand_dev_fit, dist_mean_fit)
opt_gen =  opt_gen(site_num, max_res_var)

# using funtion
rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen)


