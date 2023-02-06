# -*- coding: utf-8 -*-
"""
Master project, September 2021

Jérôme Stäheli 17-111-220

version 4.2
23.01.2023
"""
# import third party modules
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import math as math
import os
import random as rd
import platform as pf

# import custom modules
import rough_mount_fuji as rmf


#%% Generating a population

# Generate a random initial population of a fixed size, distributed uniformly
# across the genotypes
def initial_population_uniform(population_size, number_of_loci):
    number_of_genotypes = 2**number_of_loci
    population = zeros(number_of_genotypes)
    for _ in range(population_size):
        population[np.random.randint(number_of_genotypes)] += 1
    return population


#%% Generating amount of resources added per generation from a poisson distribution
def res_list_fun(res_dist_mean, res_num):
    return np.random.poisson(res_dist_mean, res_num)

#%% Fitness effect of phenotype on resource uptake
def get_phenotype(number_of_resources, fitness_influence, reference_fitness, number_of_loci, dist_mean, stand_dev, reference_genotype, name_str="", print_plot=False, save_directory=".\\"):
    phenotype_list = zeros((number_of_resources, len(rmf.genotypes_list(number_of_loci))), dtype = float)

    for resource in range(number_of_resources):
        # generate list of fitness influenced by different parameters you can choose, rmf.rough_mt_fuji described in module rough_mt_fuji
        phenotype_list[resource] = rmf.rough_mt_fuji(fitness_influence[resource], reference_fitness, number_of_loci, dist_mean, stand_dev, reference_genotype[resource], name_str, print_plot, save_directory)

    return phenotype_list


#%% Fitness calculation depending on the resource abundance, resource use of the phenotype and the number of individuals
def geno_fit_calc(resource_list, cur_pop, pop_list, phenotype):
    geno_fit = 0

    # absolute fitness of every population
    for resource in range(len(resource_list)):
        whole_pop = 0

        # calculation of the mean resource uptake of a given resource
        for genotype in range(len(pop_list)):
            # resource uptake of the whole population
            whole_pop += pop_list[genotype] * phenotype[resource][genotype]

        # Fitness calculation based on the uptake ability of a genotype compared to the whole population
        geno_fit += (phenotype[resource][cur_pop] / whole_pop) * resource_list[resource]

    return geno_fit

#%% Calculate relative fitness for every population

def cur_fit_calc(resource_list, pop_list, phenotype):
    # repeats the fitness calculations for every genotype
    cur_fit_list = zeros(len(pop_list), dtype = float)

    for genotype_pop in range(len(pop_list)):
        cur_fit_list[genotype_pop] = geno_fit_calc(resource_list, genotype_pop, pop_list, phenotype)

    # print(cur_fit_list, "fitness of the current generation")
    return cur_fit_list


#%%

# this function converts the sequence of a genotype to its index in pos_comb
def sequence_to_index(sequence, powers):
    return np.dot(sequence, powers)

# this implementation assumes that an individual gets at most one mutation per
# generation and that the mutation rate per genome per generation
def population_after_mutation(current_population, pos_comb, mutation_rate, number_of_loci):
    next_gen_pop = current_population.copy()
    population_size = np.sum(current_population)

    # auxiliary quantity to calculate indices
    powers = 2**np.flip(np.arange(number_of_loci))

    # calculate the total number of mutations that occur in the current generation
    number_of_mutations = np.random.binomial(population_size, mutation_rate)

    # choose which individuals will mutate
    individuals = zeros(population_size, dtype=int)
    idx = 0
    for genotype, n in enumerate(current_population):
        individuals[idx:idx+n] = genotype
        idx += n

    # choose which individuals will mutate
    individuals_to_mutate = np.random.choice(individuals, number_of_mutations, replace=False)

    # and mutate each of them
    for individual in individuals_to_mutate:
        # remove the individual that mutated
        next_gen_pop[individual] -= 1

        # calculate its new sequence and index
        sequence = pos_comb[individual].copy()
        locus_to_mutate = np.random.randint(number_of_loci)
        sequence[locus_to_mutate] = 1 - sequence[locus_to_mutate]

        index = sequence_to_index(sequence, powers)

        # add the individual back in the population
        next_gen_pop[individual] += 1

    return next_gen_pop


#%% Generation of the number of individuals reproducing for the next time step, with the list of resource production, list of genotype individual numbers and theire ability to take up the resources

def next_gen_calc(res_list, population, phenotype, pos_comb, mutation_rate, number_of_loci):
    # calculate the number of individuals entering the reproductive state
    fitness = cur_fit_calc(res_list, population, phenotype)

    current_fitness = fitness * population

    # make the sum of the current fitness = 1
    current_fitness = current_fitness / np.sum(current_fitness)
    next_generation = np.random.multinomial(np.sum(population), current_fitness)

    # add mutation
    next_generation_after_mut = population_after_mutation(next_generation, pos_comb, mutation_rate, number_of_loci)

    return next_generation_after_mut



#%% Repeat the calculation of the numbers of the individuals of each genotype

def consumer_resource_model(number_of_generations, res_list, pop_list, phenotype, number_of_loci, pos_comb, mut_rate):
    all_pop = zeros([number_of_generations, len(pop_list)], dtype = int)
    cur_pop_temp = pop_list[:]
    all_pop[0] = pop_list

    for cur_generation in range(number_of_generations):
        cur_pop_temp = next_gen_calc(res_list, cur_pop_temp, phenotype, pos_comb, mut_rate, number_of_loci)

        all_pop[cur_generation] = cur_pop_temp

    return all_pop

#%% Single run function

# check what system the code is running on
def get_directory():
    save_directory = "model_output"
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    if pf.system() == "Windows":
        save_directory = save_directory + "\\"
    else:
        save_directory = save_directory + "/"

    return save_directory


def get_reference_genomes(correlation_type, number_of_loci):
    # generate refrence genome based on the correlation type
    reference_genotype_1, reference_genotype_2 = None, None
    if "neg_cor" in correlation_type:
        # opposite gene for negative correlation
        reference_genotype_1 = rmf.reference_genotype(number_of_loci)
        reference_genotype_2 = np.abs(1 - reference_genotype_1)

    elif "pos_cor" in correlation_type:
        # Same optimum gene for positive correlation
        reference_genotype_1 = rmf.reference_genotype(number_of_loci)
        reference_genotype_2 = reference_genotype_1

    elif "no_cor" in correlation_type:
        # both regrences generated sepretly
        reference_genotype_1 = rmf.reference_genotype(number_of_loci)
        reference_genotype_2 = rmf.reference_genotype(number_of_loci)

    return reference_genotype_1, reference_genotype_2

#%% Save the population numbers externaly

def save_array(every_pop, name_str):
    file_name = name_str + "_every_pop"
    np.save(file_name, every_pop)



#%% Shannon index calculation

def shannon_diversity(every_population):
    # Variable to save the index in
    shannon_index = zeros(len(every_population), dtype=float)

    # Since the sum is needed repeat for every population
    for i, population in enumerate(every_population):
        total_population = np.sum(population)

        for n_genotype in population:
            # Existing populations do not partake in the calculation for the shannon index
            if n_genotype != 0:
                # Calculation of a part the sum in the shannon index consists of
                f = n_genotype / total_population
                shannon_index[i] += f * math.log(f)

    # Since the whole term for calculating is negative multiplicate with -1
    shannon_index *= -1

    # return the calculated shannon index
    return shannon_index

#%% Calculates the haplotype diversity and returns it
def haplotype_diversity(every_pop):
    hap_div = zeros(len(every_pop), dtype = float)

    for i, generation in enumerate(every_pop):
        hap_freq = 0
        total_population = np.sum(generation)

        for gen_type in generation:
            hap_freq += (gen_type / total_population) ** 2

        hap_div[i] = (total_population / (total_population - 1)) * (1 - hap_freq)

    return hap_div

#%% Nucleotide diversity calculation and returns it
def nucleotide_diversity(every_pop, pos_comb):
    nuc_div = zeros(len(every_pop), dtype = float)

    for generation in range(len(every_pop)):
        geno_site_sum = 0

        total_population = np.sum(every_pop[generation])

        for genotype1 in range(len(every_pop[0]) - 1):
            for genotype2 in range(genotype1 + 1, len(every_pop[0])):
                freq_genotype_1 = every_pop[generation][genotype1] / total_population
                freq_genotype_2 = every_pop[generation][genotype2] / total_population
                nuc_dif = np.sum(abs(pos_comb[genotype1] - pos_comb[genotype2]))

                geno_site_sum += 2 * freq_genotype_1 * freq_genotype_2 * nuc_dif

        nuc_div[generation] = (np.sum(every_pop) / (np.sum(every_pop) - 1)) * geno_site_sum

    return nuc_div


#%% Plot single runs

def plot_single_run(every_pop, save_as_string, pos_comb, save_directory):
    # Make a plot that shows the development of each genotype and also the sum of individuals
    #plt.plot(np.sum(every_pop,axis=1),".")
    plt.plot(every_pop)

    # Lables of population development plot
    plt.title("Development of Genotypes")
    plt.xlabel("Generations")
    plt.ylabel("Number of individuals")

    # Save and show fig
    plt.savefig(save_directory + save_as_string + "_plot.pdf")
    plt.show()
    plt.close("all")

    # Calculate the shannon index and add it to new plot
    shannon_ind = shannon_diversity(every_pop)
    plt.plot(shannon_ind, label = "Shannon-Index")

    # Calculate the haplotype diversity and add it to plot
    hap_div = haplotype_diversity(every_pop)
    plt.plot(hap_div, label = "haplotype diversity")

    # Calculate the nucleotype diversity and add it to plot
    nuc_dif = nucleotide_diversity(every_pop, pos_comb)
    plt.plot(nuc_dif, label = "nucleotide diversity")

    # Lables of diversity plot
    plt.title("Meassures of diversity in whole Population")
    plt.xlabel("Generations")
    plt.ylabel("Degree of diversity")

    # add Legends
    plt.legend(loc= "upper right")

    # save plot to pdf and show it
    plt.savefig(save_directory + save_as_string + "_stats.pdf")
    plt.show()
    plt.close("all")

#%% Plot multi runs with the shannon index, haplotype diversity and the nucleotide diversity and the standard diviation of each

def plot_multi_run(num_generation, shannon_index_list, hap_div_list, nuc_div_list, save_as_string):
    
    plt.fill_between(range(num_generation),
                     np.mean(shannon_index_list, axis = 0) + np.std(shannon_index_list, axis = 0),
                     np.mean(shannon_index_list, axis = 0) - np.std(shannon_index_list, axis = 0),
                     alpha = 0.3)
    plt.fill_between(range(num_generation),
                     np.mean(hap_div_list, axis = 0) + np.std(hap_div_list, axis = 0),
                     np.mean(hap_div_list, axis = 0) - np.std(hap_div_list, axis = 0),
                     alpha = 0.3)
    plt.fill_between(range(num_generation),
                     np.mean(nuc_div_list, axis = 0) + np.std(nuc_div_list, axis = 0),
                     np.mean(nuc_div_list, axis = 0) - np.std(nuc_div_list, axis = 0),
                     alpha = 0.3)

    plt.xlim(0, num_generation - 1)
    plt.ylim(0, np.max(shannon_index_list) + 0.5)

    plt.plot(np.mean(shannon_index_list, axis = 0))
    plt.plot(np.mean(hap_div_list, axis = 0))
    plt.plot(np.mean(nuc_div_list, axis = 0))

    # Lables of population development plot
    plt.title("Development of Genotypes")
    plt.xlabel("Generations")
    plt.ylabel("Number of individuals")

    plt.legend(["Shannon-Index", "Haplotype diversity", "Nucleotide diversity"])
    plt.savefig(save_as_string + "_multirun.pdf")
    plt.show()
    plt.close("all")

#%% Function that plots and runs a single run with the chosen parameter and saves it in the correlation_type string


def single_run(correlation_type, pop_list, res_list, num_generation, mut_rate, reference_fitness, number_of_loci, dist_mean, stand_dev, stand_dev_fit, dist_mean_fit):
    print_plot = True

    save_directory = get_directory()
    name_str = "CRM_" + "single_run_" + correlation_type

    reference_genotype_1, reference_genotype_2 = get_reference_genomes(correlation_type, number_of_loci)

    # Fitness influence of a site
    fitness_influence_1 = rmf.fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit)
    fitness_influence_2 = rmf.fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit)

    # using funtion to get survivability of the different
    genotypes_list = rmf.genotypes_list(number_of_loci)

    # Generate the ability to take up a resource by the different genotypes the first index since it creates a array of array. Generates array of arrays if the other variables are all the same
    phenotype = get_phenotype(2, [fitness_influence_1, fitness_influence_2], reference_fitness, number_of_loci, dist_mean, stand_dev, [reference_genotype_1, reference_genotype_2], name_str + "_1", print_plot, save_directory)

    # print phenotype correlation
    print('phenotype correlation: ', np.corrcoef(phenotype[0], phenotype[1])[0,1])

    # Run the model for num_generation generation
    every_pop = consumer_resource_model(num_generation, res_list, pop_list, phenotype, number_of_loci, genotypes_list, mut_rate)

    save_array(every_pop, save_directory + name_str)

    # Plot single run with negative correlated resource uptake
    plot_single_run(every_pop, name_str, genotypes_list, save_directory)


#%% Multi run function

def multi_run(correlation_type, pop_list, res_list, number_of_loci, num_generation, mut_rate, stand_dev_fit, dist_mean_fit, reference_fitness, dist_mean, stand_dev, rep_sim, print_plot = False):
    reference_genotype_1 = [0] * number_of_loci

    # empty array for saving
    shannon_index_list = zeros([rep_sim, num_generation])
    hap_div_list       = zeros([rep_sim, num_generation])
    nuc_div_list       = zeros([rep_sim, num_generation])
    some_pop_list      = zeros([rep_sim, num_generation, 2 ** number_of_loci])
    phenotype_list       = zeros([rep_sim, 2, 2 ** number_of_loci])

    genotypes_list = rmf.genotypes_list(number_of_loci)
    
    # run the model
    for repeat in range(rep_sim):
        save_directory = get_directory()
        name_str = "CRM_" + "multi_run_" + correlation_type

        # Fitness influence of a site
        fitness_influence_1 = rmf.fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit)
        fitness_influence_2 = rmf.fitness_influence(number_of_loci, stand_dev_fit, dist_mean_fit)
        
        # refrence genome for fitness calculation
        reference_genotype_1, reference_genotype_2 = get_reference_genomes(correlation_type, number_of_loci)
        
        # rough Mt. fuji generated phenotype
        phenotype = get_phenotype(2, [fitness_influence_1, fitness_influence_2], reference_fitness, number_of_loci, dist_mean, stand_dev, [reference_genotype_1, reference_genotype_2], name_str + "", print_plot, save_directory)

        ### print phenotype correlation
        print('phenotype correlation: ', np.corrcoef(phenotype[0], phenotype[1])[0,1])

        some_pop = consumer_resource_model(num_generation, res_list, pop_list, phenotype, number_of_loci, genotypes_list, mut_rate)
        
        # save the results inside the array
        some_pop_list[repeat]      = some_pop
        phenotype_list[repeat]     = phenotype
        nuc_div_list[repeat]       = nucleotide_diversity(some_pop, genotypes_list)
        hap_div_list[repeat]       = haplotype_diversity(some_pop)
        shannon_index_list[repeat] = shannon_diversity(some_pop)
        
    # saves the arrays in the directory with a name that is determined by the correlation type and the parameters chosen
    save_array(phenotype_list, save_directory + name_str + "_phenotype_list")
    save_array(some_pop_list, save_directory + name_str + "_pop_list")
    save_array(nuc_div_list, save_directory + name_str + "_nuc_div_list")

    plot_multi_run(num_generation, shannon_index_list, hap_div_list, nuc_div_list, save_directory + name_str)

#%% Parameters for for model

"""
# Number of sites that can differ
number_of_loci = 5
# Population distribution mean
pop_dist_mean  = 20
# population size
population_size = 200
# List of population size
#pop_list       = pop_list_fun(pop_dist_mean, number_of_loci)
pop_list = initial_population_uniform(population_size, number_of_loci)
# Number of resources
res_num        = 2
# Ressours distribution mean
res_dist_mean  = 50
# ressurce use list
res_list       = res_list_fun(res_dist_mean, res_num)
# Fitness of the optimal Genotype
reference_fitness        = 3
# Standart deviation of non additive effects
stand_dev      = 0.2
# Distribution mean of non additive effects
dist_mean      = 0
# Standart distribution for the fitness influence of the site
stand_dev_fit  = 0.2
# Distribution mean of the fitness influence of the site
dist_mean_fit  = -0.2
# Mutation rate, chance for a individual to mutate to a other genotype
mutation_rate       = 4 * 10 ** (- 3)
# threshhold of individual with which a population counts as maintained
num_generation = 500
# Number of repeats
rep_sim        = 20

correlation_type = 'neg_cor'
"""

"""
single_run(correlation_type, pop_list, res_list, num_generation, mutation_rate, reference_fitness, number_of_loci, dist_mean, stand_dev, stand_dev_fit, dist_mean_fit)
"""
#%%
