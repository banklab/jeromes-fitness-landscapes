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
from numpy import zeros
from numpy import array
import matplotlib.pyplot as plt
from math import factorial


#%% Parameters

# Allele Frequencies
allele_freq_1 = 0.5
# Population size
pop_start = 100
# Number of generations
gen_time = 5000


#%% Probability of inheritance

# Calculate probability with the binominal coefficient
def prob_pop(pop_size, all_freq_1):
    
    # Create empty list for saving every chance to pass on certain amount of copies of one allel
    prob_list = zeros((2 * pop_size) + 1, dtype = float)
    
    
    # Fisher calculation for the chance to pass on copies of a allele
    for num_prob in range(0, (2 * pop_size) + 1):
        prob_list[num_prob] = (((factorial(2 * pop_size))
                               / (factorial(num_prob) * factorial(2 * pop_size-num_prob)))
                               * (all_freq_1 ** num_prob)
                               * ((1 - all_freq_1) ** (2 * pop_size-num_prob)))

    return prob_list


#%% Allele frequency of the next generation

# Calculate allele frequencies of next generation
def next_gen(pop_size, all_freq_1):

    # Choose a amount of alleles from one type of alleles with the calculated probability of prob_pob function
    next_gen = np.random.choice(range(0, (2 * pop_size) + 1), 
                                1, 
                                p = prob_pop(pop_size, all_freq_1))
    all_freq_new = next_gen / (pop_size * 2)
    
    # print(all_freq_new)
    
    return all_freq_new[0]


#%% Whright Fisher model

# Repeat the calculation of the new allele frequencies a chosen number of times
def wright_fisher_model(gen_num, pop_size, all_freq_1):
    
    # Repeat the next_gen function to generate new allele frequencies from the previous allele frequencies
    all_freq = zeros(gen_num, dtype = float)
    all_freq[0] = all_freq_1
    
    
    for rep_calc in range(1, gen_num):
        
        # Add new allele frequency to frequency list
        all_freq[rep_calc] = next_gen(pop_size, all_freq[rep_calc-1])
        
        # Stop looping when population reaches 0 or 1
        if all_freq[rep_calc] == 0 or all_freq[rep_calc] == 1:
            break
        
    all_freq = all_freq[all_freq != 0]
        
      
    # print(all_freq)
    return all_freq


#%% Checking if the funktion is working
"""

allele_frequencies = wright_fisher_model(gen_time, pop_start, allele_freq_1)

plot_range = array(range(0, len(allele_frequencies)), dtype = int)
plt.plot(plot_range, allele_frequencies)

"""