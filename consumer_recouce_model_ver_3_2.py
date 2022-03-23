
"""
The REAL seascape project, September 2021

Jérôme Stäheli 17-111-220

version 0.8
27.01.2021
"""
# import third party modules
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import math as math

# import custom modules
import rough_mt_fuji_ver_3_1 as rmf
import wright_fisher_ver_3_1 as wfm



#%% Generating an population from a poisson distribution

def pop_list_fun(pop_dist_mean, site_num):
    
    return np.random.poisson(pop_dist_mean, 2 ** site_num)

#%% Generating amount of recources added per generation from a poisson distribution

def res_list_fun(res_dist_mean, res_num):
    
    return np.random.poisson(res_dist_mean, res_num)

#%% Fitness effect of phenotype on recource uptake

def res_use(res_num, fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen):
    
    res_use_list = zeros((res_num, len(rmf.pos_comb_list(max_res_var, site_num))), dtype = float)
    
    for resource in range(0, res_num):
        
        # generate list of fitness influenced by diffrent parameters you can chose, rmf.rough_mt_fuji described in module rough_mt_fuji
        res_use_list[resource] = rmf.rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen)[0]
        
    # print(res_use_list)
    
    return res_use_list


#%% mean fitness calculation

def mean_fit_calc(res_list, pop_list):
    
    mean_fit = np.sum(res_list) / np.sum(pop_list)
    
    # print(mean_fit, "mean fitness of all genotypes")
    
    return mean_fit



#%% Fitnes calculation depending on the resource abundance, resurce use of the phenotype and the number of indivituals

def geno_fit_calc(res_list, cur_pop, pop_list, res_use_list):
    
    geno_fit = 0
    
    # absolute fitness of every population
    for resource in range(len(res_list)):
        
        whole_pop = 0
        
        # calculation of the mean recource uptake of a given resource
        for genotype in range(len(pop_list)):
            
            # resource uptake of the whole population
            whole_pop += pop_list[genotype] * res_use_list[resource][genotype]
        
        # Fittnes calculation based on the uptake ability of a genotype compared to the whole population
        geno_fit += (res_use_list[resource][cur_pop] / whole_pop) * res_list[resource]

    
    # print(geno_fit, "fitness of the current genotype")
    
    return geno_fit

#%% Relative fitness calculation

def rel_fit_calc(res_list, cur_pop, pop_list, res_use_list):
    
    # absoluta fitness calculation
    geno_fit = geno_fit_calc(res_list, cur_pop, pop_list, res_use_list)
    
    # mean fitness of the population
    mean_fit = mean_fit_calc(res_list, pop_list)
    
    # mean fitness calculation by dividing the absolute fitness by the mean fitness
    rel_fit = geno_fit / mean_fit

    # print(rel_fit, "relative fitness")

    return rel_fit

#%% Callculate relative fitness for every population

def cur_fit_calc(res_list, pop_list, res_use_list):
    
    
    # repeats the fitness calculations for every genotype
    cur_fit_list = zeros(len(pop_list), dtype = float)
    
    for genotype_pop in range(len(pop_list)):
        
        cur_fit_list[genotype_pop] = rel_fit_calc(res_list, genotype_pop, pop_list, res_use_list)
    
    # print(cur_fit_list, "fitness of the current generation")
    
    return cur_fit_list


#%% Content of the population after mutation events


def mut_direction_calc(pos_comb, mut_rate, site_num):
    
    
    mut_direction = zeros((len(pos_comb), site_num), dtype = int)
    
    
    for comb1 in range(0, len(pos_comb)):
        
        
        # serging for one mutation step between population
        for comb2 in range(0, len(pos_comb)):
            
            count3 = 0
            
            for site in range(0, len(pos_comb[comb2])):
                
                if pos_comb[comb1][site] == pos_comb[comb2][site]:
                    
                    count3 += 1
            
            if count3 == 2:
                
                mut_direction[comb1][site] = comb2
                
        return mut_direction

#%%
        
def pop_after_mut(cur_pop_list, pos_comb, mut_rate, site_num, mut_direction):

    next_gen_pop = cur_pop_list[:].copy()
        
    for comb1 in range(0, len(pos_comb)):
        
        
        # generates the new population after mutation events
        for mut_pop in range(len(mut_direction[comb1])):
            
            count2 = 0
            
            for genotype_pop in range(0, cur_pop_list[comb1]):
                
                count2 += int(np.random.choice([0, 1], p = [1 - mut_rate, mut_rate]))
            
            next_gen_pop[comb1] -= count2
            next_gen_pop[mut_direction[comb1][mut_pop]] += count2
            
    #print(cur_pop_list, pos_comb_2, mut_rate, next_gen_pop)
    
    return next_gen_pop
    
#pop_after_mut([200, 200, 200, 200, 200, 200, 200, 200], [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)], 0.01)


#%% Generation of the number of individuals reproducing for the next time step, with the list of resource production, list of genotype individual numbers and theire ability to take up the resources

def next_gen_calc(res_list, pop_list, res_use_list, pos_comb, mut_rate, site_num, mut_direction):
    
    # calculate the number of individuals entering the reproductive state
    cur_fit_list = cur_fit_calc(res_list, pop_list, res_use_list)
    
    # print(cur_fit_list, "current fitness list")
    # print(pop_list, "current population size list")
    
    # new fitness list
    cur_fit = zeros(len(cur_fit_list), dtype = float)
    
    # make the sum of the current fitness = 1
    for gen_fit in range(len(cur_fit_list)):
        
        cur_fit[gen_fit] = cur_fit_list[gen_fit] / sum(cur_fit_list)
    
    #print(cur_fit, sum(cur_fit), "fitnesses of the current generation, should be fit = 1")
    
    
    
    
    cont_next_gen = zeros(len(pop_list), dtype = int)
    
    # Add the number of individuals to a list whith all genotype individual numbers of the current time step
    for pop in range(len(pop_list)):
        
        new_disc_pop = 0
        
        
        # each individual has the calculated chance to survive after resource consumption
        for individual in range(pop_list[pop]):
        
            new_disc_pop += int(np.random.choice([0, 1], p = [1 - cur_fit[pop], cur_fit[pop]]))
    
        cont_next_gen[pop] = new_disc_pop                    
        
    # print(cont_next_gen, "gives out the survivors of the resouces competition")
    
    
    
    
    
    allel_freq_next_gen = []
    
    # Put the new number of Individuals that survived of each genotype through the wright fisher model
    for num12 in range(len(cont_next_gen)):
                
        allel_freq_next_gen.append(wfm.next_gen(sum(cont_next_gen), cont_next_gen[num12]/sum(cont_next_gen)))
    
    # make that frequencies add up to one
    for num16 in range(len(allel_freq_next_gen)):
        
        allel_freq_next_gen[num16] = allel_freq_next_gen[num16] / sum(allel_freq_next_gen)
    
    # print(allel_freq_next_gen, "allel freq, share of every genotyp in the next gen")
    
    
      
    
    # print(sum(res_list), sum(pop_list), "sum of all populations and resources")
    
    
    cont_next_gen_after_fisher = []
    
    # now determine next generation with the help of the allele frequency and the number of individuals of the previous generation
    for num14 in range(len(pop_list)):
        
        next_gen_pop_count = 0
        
        for num16 in range(sum(res_list)):
        
            next_gen_pop_count += int(np.random.choice([0, 1], p = [1 - allel_freq_next_gen[num14], allel_freq_next_gen[num14]])) 
            
        
        cont_next_gen_after_fisher.append(next_gen_pop_count)
    
    # print(cont_next_gen_after_fisher, "after fisher")
    
    
    
    # add mutation in
    
    cont_next_gen_after_fisher = pop_after_mut(cont_next_gen_after_fisher, pos_comb, mut_rate, site_num, mut_direction)
    
    
    
    return cont_next_gen_after_fisher
    



    
    



#%% Repeat the calculation of the numbers of the individuals of each genotype

def consumer_recource_model(gen_time, res_list, pop_list, res_use_list, site_num, pos_comb_3, mut_rate, mut_direction):
    
    all_pop = []
    cur_pop_temp = pop_list[:]
    all_pop.append(pop_list)
    
    
    for num6 in range(1, gen_time):
        
        cur_pop_temp = next_gen_calc(res_list, cur_pop_temp, res_use_list, pos_comb_3, mut_rate, site_num, mut_direction)
        
        all_pop.append(cur_pop_temp)
        
        
        ext_genotype = 0
        
        for num17 in range(len(cur_pop_temp)):
            
            if cur_pop_temp[num17] == 0:
                
                ext_genotype +=1
        
        
        if ext_genotype == 2 ** site_num - 1:
            
            return all_pop
            
    
                
    
    
    
    # print(all_pop)
    
    return all_pop
    

#%% Plot the genotype individuals numbers

def plot_mft_crm(every_pop):

    sort_pop = []
    
    for num7 in range(len(every_pop[0])):
        
        sort_cur_pop = []
        
        for num8 in range(len(every_pop)):
            
            sort_cur_pop.append(every_pop[num8][num7])
            
        sort_pop.append(sort_cur_pop)
        
    # print(len(sort_pop), len(sort_pop[1]))
    
    for num9 in range(len(sort_pop)):
        plt.plot(sort_pop[num9])

#%% plot the whole population to show that it does not change

def mean_mft_crm(every_pop):
    
    pop_mean = []
    
    for num10 in range(len(every_pop)):
        
        cur_whole_pop = 0
        
        for num11 in range(len(every_pop[num10])):
            
            cur_whole_pop += every_pop[num10][num11]
        
        pop_mean.append(cur_whole_pop)
            
    
    # print(pop_mean)
    plt.plot(pop_mean, ".")
    return pop_mean

#%% Shannon index calculation

def shannon_calc(pop_sizes):
    
    # Vaiable to save the index in
    shannon_index = 0
    
    # Since the sum is needed repeat for every population
    for num25 in pop_sizes:
        
        # Extingt populations do not partake in the calculation for the shannon index
        if num25 != 0:
            
            # Callculation of a part the sum in the shannon index consists off
            shannon_index += num25 / sum(pop_sizes) * math.log(num25 / sum(pop_sizes))
       
    # Since the whole term for calculating is negative multiplicate with -1
    shannon_index *= -1
    
    # return the calculated shannon index
    return shannon_index
    
    
    

#%% 

#Parameters for fit_infl, opt_gen and rought_mount_fuji


site_num = 3
# Population distribution mean
pop_dist_mean = 200
# Number of resources
res_num = 2
# Ressours distribution mean
res_dist_mean = 700
# Number of possible aminoacides on each site
max_res_var     = 2
# Fitness of the optimal Genotype
opt_fit         = 1
# Standart deviation of non additive effects
stand_dev       = 0.2
# Distribution mean of non additive effects
dist_mean       = 0
# Standart distribution for the fitness influence of the site
stand_dev_fit   = 0.2
# Distribution mean of the fitness influence of the site
dist_mean_fit   = -0.5
# Mutation rate, chance for a individual to mutate to a other genotype
mut_rate = 4.5 * 10 ** (- 9)
# threshhold of individual with which a population counts as maintained
ind_thresh = 20


# fitness influence (fit_infl) and optimal genotype (opt_gen) are needed for the Rough Mt Fuji function
# might later on include it in the funtion
fit_infl = rmf.fit_infl(site_num, stand_dev_fit, dist_mean_fit)
opt_gen =  rmf.opt_gen(site_num, max_res_var)


# empty list to fill with everything
fit_pos_comb_mut_list = []

# using funtion to get survivability of the diffrent
pos_comb = rmf.rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen)[1]
mut_list = rmf.rough_mt_fuji(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen)[2]

print(mut_list, "list with numbers of mutation")
print(pos_comb, "list of all possible combinations")




res_use_list = [[], []]
res_use_list[0] = (res_use(1, fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, [1, 1, 1])[0])
res_use_list[1] = (res_use(1, fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, [2, 2, 2])[0])
print(res_use_list, type(res_use_list))



res_list     = res_list_fun(res_dist_mean, res_num)
pop_list     = pop_list_fun(pop_dist_mean, site_num)

mut_direction = mut_direction_calc(pos_comb, mut_rate, site_num)



every_pop = consumer_recource_model(10, res_list, pop_list, res_use_list, site_num, pos_comb, mut_rate, mut_direction)

print(every_pop[-1])
print(mut_list)
print(shannon_calc(every_pop[-1]))

mean_mft_crm(every_pop)
plot_mft_crm(every_pop)


#%% run the it multiple times and look at how many times it survived



"""
count4 = 0

shannon_index_list = []

for num23 in range(0, 200):
    
    count5 = 0
    
    res_use_list = [[], []]
    res_use_list[0] = (res_use(1, fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, [1, 1, 1])[0])
    res_use_list[1] = (res_use(1, fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, [2, 2, 2])[0])
    
    res_list     = res_list_fun(res_dist_mean, res_num)
    pop_list     = pop_list_fun(pop_dist_mean, site_num)
    
    some_pop = consumer_recource_model(500, res_list, pop_list, res_use_list, site_num, pos_comb, mut_rate)[-1]
    
    for num24 in some_pop:
        
        if num24 > ind_thresh:
            
            count5 += 1
            
    if count5 > 1:
        
        count4 += 1
    
        
    shannon_index_list.append(shannon_calc(some_pop))
    
    
        
print(count4)
"""          



