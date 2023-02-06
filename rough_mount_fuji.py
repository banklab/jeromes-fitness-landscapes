# -*- coding: utf-8 -*-
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
from numpy import array
import matplotlib.pyplot as plt
rng = np.random.default_rng()


#%% Fitness influence of locus

# Creates fitness influence value, these are taken from a normal distribution
def fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit):
    return rng.normal(dist_mean_fit, stand_dev_fit, number_of_loci)


#%% Optimal genome
"""
generates random numbers for every locus to generate the reference genome
with the help of a defined number of biallelic loci
"""
def reference_genotype(number_of_loci):
    # Empty list to save every position of the optimal gene
    return rng.integers(2, size = number_of_loci)


def genotypes_list(number_of_loci):
    # generates all possible combinations with the help of number of loci
    return array(list(product({0, 1}, repeat = number_of_loci)))


#%% Rough Mt. Fuji-Type fitness landscape
"""
Generates the Rough Mt Fuji-Type fitness landscape with the fitness influence of every locus, the number
of loci, the mean and the standart deviation
of the normal distribution for non additive effects, the optimal genotype and its fitness
"""
def rough_mt_fuji(fitness_influence, reference_fitness, number_of_loci, dist_mean, stand_dev, reference_genotype, name_str = "", print_plot=False, save_direct = ".\\"):
    # generates all possible genotypes
    pos_comb = genotypes_list(number_of_loci)

    recoded_pos_comb = np.abs(reference_genotype - pos_comb)

    # generates the fitness of all possible genotypes and saves them in a list
    #(with the help of the fitness calculation)
    fitness_list = rng.normal(dist_mean, stand_dev, len(pos_comb)) - np.dot(recoded_pos_comb, fitness_influence) + reference_fitness

    """
    for comb in range(0, len(pos_comb)):
        fitness_list[comb] += add_fit_calc(fitness_influence, reference_fitness, reference_genotype, pos_comb[comb])
    """

    if print_plot:
        plot_fitness_landscape(fitness_list, number_of_loci, pos_comb, save_direct + "RMF_" +  name_str + ".pdf")

    fitness_list = np.exp(fitness_list)

    # return the list of possible combinations and the corresponding fitnesses
    return fitness_list


def plot_fitness_landscape(fitness_landscape, number_of_loci, pos_comb, filename):
    # generates a list with the number of diffrent sites compared to a all non mutated genotype
    # this is just used for the creation of a plot
    mut_list = zeros(len(pos_comb), dtype = int)

    for comb in range(0, len(pos_comb)):
        num_mut  = 0
        for site in range(0, len(pos_comb[comb])):
            if pos_comb[comb][site] != 1:
                num_mut += 1

        mut_list[comb] = num_mut

    plt.figure(figsize=(8, 4))

    # adds lines in the plot that connect all the dots of each genotype
    for genotype1 in range(0, len(pos_comb) - 1):
        for genotype2 in range(genotype1 + 1, len(pos_comb)):
            if np.abs(pos_comb[genotype1] - pos_comb[genotype2]).sum() == 1:
                plt.plot([mut_list[genotype1] - 0.2, mut_list[genotype2] + 0.2],
                         [fitness_landscape[genotype1], fitness_landscape[genotype2]],
                         'k-')

    # plots the genotypes with the number of mutations and the lines that indicates to which genotype they
    # can further mutate
    plt.plot(mut_list, fitness_landscape, 'ro', color = "white")

    plt.xlim(-0.3, number_of_loci + 0.3)

    for points in range(len(pos_comb)):
        plt.annotate(pos_comb[points],
                     (mut_list[points] - 0.14,
                      fitness_landscape[points] - 0.005))

    plt.rcParams.update({"axes.titlesize" : 20})
    plt.rcParams.update({"axes.labelsize" : 17})
    plt.rcParams.update({"lines.linewidth" : 3})
    plt.rcParams.update({"xtick.labelsize" : 14})
    plt.rcParams.update({"ytick.labelsize" : 14})

    plt.xlabel("Genotypes")
    plt.ylabel("Fitness")
    plt.title("Rough Mt. Fuji")

    plt.savefig(filename)

    plt.show()
    plt.close("all")



#%% Parameters for fitness_influence, reference_genotype and rought_mount_fuji
'''
# Number of loci that have influence on fitness and can mutate
number_of_loci  = 3
# Fitness of the reference genotype
reference_fitness = 1.
# Standart deviation of non additive effects
stand_dev       = 0.1
# Distribution mean of non additive effects
dist_mean       = 0.
# Standart distribution for the fitness influence of the site
stand_dev_fit   = 0.1
# Distribution mean of the fitness influence of the site
dist_mean_fit   = 0.
'''

#%% Use funtions

'''
# fitness influence (fitness_influence) and optimal genotype (reference_genotype) are needed for the Rough Mt Fuji function
# might later on include it in the funtion
fitness_influence_rmf = fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit)
reference_genotype_rmf =  reference_genotype(number_of_loci)

# using funtion
rough_mt_fuji(fitness_influence_rmf,
              reference_fitness,
              number_of_loci,
              dist_mean,
              stand_dev,
              reference_genotype_rmf,
              print_plot=True,
              save_direct='')
'''
