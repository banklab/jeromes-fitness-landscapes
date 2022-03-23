#%%
"""
Rough Mt. Fuji-Type fitness landscapes as described in Aita T, Hidefumi U, 
inaoka T, Nakajima M, Kokubo T and Husimi Y (2000) Analysis of a Local Fitness Landscape with a Model of
the Rough Mt. Fuji-Type Landscape: Application to Prolyl Endopeptidase and Thermolysin.  
Biopolymers 54: 64–79

Jérôme Stäheli 17-111-220

version 2.4
30.11.2021
"""

#%% Packages
"""
Load Packages (matplotlib.pyplot for plots, Numpy for normal distribution, random for random numbers 
and intertools for all possible gene combinations)
"""

from itertools import product
import numpy as np
from numpy import zeros
from numpy import asarray
from numpy import full
from numpy import array
import matplotlib.pyplot as plt
rng = np.random.default_rng()


#%% Fitness influence of locus

# Creates fitness influence value, these are taken from a normal distribution
def fit_infl(site_num, stand_dev_fit, dist_mean_fit):
        
    return rng.normal(dist_mean_fit, stand_dev_fit, site_num) 


#%% Optimal genome
"""
generates random numbers for every site to generate the optimal genome 
with the help of a defined number of sites and a number of possible variations
"""
def opt_gen(site_num, max_res_var):
    
    # Empty list to save every position of the optimal gene
    return rng.integers(1, max_res_var, size=site_num)


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
    for site in range(0, site_num):
        
        if cur_gen[site] != opt_gen[site]:
            cur_fit -= fit_infl[site]
    
    """printing diffrent things to see what they are
    print(cur_fit)
    """
    return cur_fit


#%%

def pos_comb_list(max_res_var, site_num):
    
    # generates all possible combination with the help of site number and the number of possibilities of 
    # each site
    
    return array(list(product({1, max_res_var}, repeat = site_num)))


#%% Rough Mt. Fuji-Type fitness landscape
"""
Generates the Rough Mt Fuji-Type fintnes landscape with the fitness influence of every site, the number
of sites, how many diffrent things can be expressed on one site, the mean and the standart deviation 
of the normal distribution for non additive effects, the optimal genotype and its fittnes
"""
def rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen, print_plot = "No"):
    
    # generates all possible combination with the help of site number and the number of possibilities of 
    # each site
    pos_comb = pos_comb_list(max_res_var, site_num)
        
    # generates the fitness of all possible genotypes and saves them in a list 
    #(with the help of the fitness calculation)
    #fit_list = zeros(len(pos_comb), dtype = float)
    
    fit_list = rng.normal(dist_mean, stand_dev, len(pos_comb))
        
    for comb in range(0, len(pos_comb)):
        
        fit_list[comb] += add_fit_calc(fit_infl, opt_fit, opt_gen, pos_comb[comb])
    
    # generates a list with the number of diffrent sites compared to a all non mutated genotype
    # this is just used for the creation of a plot
    mut_list = zeros(len(pos_comb), dtype = int)

    for comb in range(0, len(pos_comb)):
        
        num_mut  = 0
        
        for site in range(0, len(pos_comb[comb])):
            
            if pos_comb[comb][site] != 1:
                num_mut += 1
    
        mut_list[comb] = num_mut
    
    """
    printing diffrent things to see what they are
    print(pos_comb, "list of all possible combinations")
    """
    
    if print_plot == "yes":
    
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
        plt.show()
    
    """
    printing diffrent things to see what they are
    print(fit_list, "list with final fitnesses")
    print(mut_list, "list with numbers of mutation")
    print(len(mut_list), "list with possible mutants")
    """
    
    # retun of the list of possible combinations and the corresponding fitnesses
    return fit_list, pos_comb, mut_list
    

#%% Parameters for fit_infl, opt_gen and rought_mount_fuji

# Number of sites that have influence on fitness and can mutate
site_num        = 3
# Number of possible aminoacides on each site
max_res_var     = 2
# Fitness of the optimal Genotype
opt_fit         = 1
# Standart deviation of non additive effects
stand_dev       = 0.1
# Distribution mean of non additive effects
dist_mean       = 0
# Standart distribution for the fitness influence of the site
stand_dev_fit   = 0.1
# Distribution mean of the fitness influence of the site
dist_mean_fit   = 0


#%% Use funtions
"""
# fitness influence (fit_infl) and optimal genotype (opt_gen) are needed for the Rough Mt Fuji function
# might later on include it in the funtion
fit_infl_rmf = fit_infl(site_num, stand_dev_fit, dist_mean_fit)
opt_gen_rmf =  opt_gen(site_num, max_res_var)

# using funtion
rough_mt_fuji(fit_infl_rmf, 
              opt_fit, 
              site_num, 
              max_res_var, 
              dist_mean, 
              stand_dev, 
              opt_gen_rmf, 
              "yes")

"""