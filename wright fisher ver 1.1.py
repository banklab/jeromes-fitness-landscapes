#%%
"""
Wright Fisher model based on the describtion of Tran TD, Hofrichter J and Jost J (2013) 
An introduction to the mathematical structure of the Wright–Fisher model of population genetics. 
Theory Bioscience 132: 73–82

Jérôme Stäheli 17-111-220

version 1.1
30.11.2021
"""

#%% Packages

# Load Packages
import numpy as np
import matplotlib.pyplot as plt


#%% Parameters

# Allele Frequencies
allele_freq_1 = 0.5
# Population size
pop_start = 500
# Number of generations
gen_time = 5000

#%% Faculty

# Calculate the faculty
def fac(fac_num):
    fac = 1
    for num_fac in range(1, fac_num + 1):
        fac = fac * num_fac
    return fac


#%% Probability of inheritance

# Calculate probability with the binominal coefficient
def prob_pop(pop_size, all_freq_1):
    
    # Create empty list for saving every chance to pass on certain amount of copies of one allel
    prob_list = []
    
    # Fisher calculation for the chance to pass on copies of a allele
    for num_prob in range(0, (2 * pop_size) + 1):
        prob = (((fac(2 * pop_size))
                / (fac(num_prob) * fac(2 * pop_size-num_prob)))
                * (all_freq_1 ** num_prob)
                * ((1 - all_freq_1) ** (2 * pop_size-num_prob)))
        prob_list.append(prob)

    return prob_list


#%% Allele frequency of the next generation

# Calculate allele frequencies of next generation
def next_gen(pop_size, all_freq_1):

    # Choose a amount of alleles from one type of alleles with the calculated probability of prob_pob function
    next_gen = np.random.choice(list(range(0, (2 * pop_start) + 1)), 
                                1, 
                                p = prob_pop(pop_start, all_freq_1))
    all_freq_new = next_gen / (pop_size * 2)
    
    return all_freq_new[0]


#%% Whright Fisher model

# Repeat the calculation of the new allele frequencies a chosen number of times
def wright_fisher_model(gen_num, pop_size, all_freq_1):
    
    # Repeat the next_gen function to generate new allele frequencies from the previous allele frequencies
    all_freq = [all_freq_1]
    
    for rep_calc in range(1, gen_num):
        
        # Add new allele frequency to frequency list
        all_freq.append(next_gen(pop_size, all_freq[rep_calc-1]))
        
        # Terminate function when reaching 1 or 0 since it does not change anymore (to save time)
        if all_freq[rep_calc-1] == 0 or all_freq[rep_calc-1] == 1:
            break
        
    print(all_freq)
    return all_freq


#%% Checking if the funktion is working

allele_frequencies = wright_fisher_model(gen_time, pop_start, allele_freq_1)
plt.plot(range(0, len(allele_frequencies)), allele_frequencies)

