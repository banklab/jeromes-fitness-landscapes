
import random as rd
from itertools import product
import numpy as np
import matplotlib.pyplot as plt

site_num    = 4
max_res_var = 2
min_fit     = max_res_var / 2
opt_fit     = 2
fit_infl    = min_fit / site_num
stand_dev   = 0.1
dist_mean   = 0


def opt_gen(site_num, max_res_var):
    
    opt_gen = []
    for num in range(0, site_num):
        opt_gen.append(rd.randint(1, max_res_var))
    
    print(opt_gen, "optimal gene")
    return opt_gen


def add_fit_calc(fit_infl, opt_fit, opt_gen, cur_gen):
    
    cur_fit = opt_fit
    
    for num in range(0, site_num):
        
        if cur_gen[num] != opt_gen[num]:
            cur_fit -= fit_infl
        
    print(cur_fit)
    return cur_fit


def fit_calc(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen):
    
    pos_comb = list(product({1, max_res_var},repeat = site_num))
    
    fit_list = []
    
    for num2 in range(0, len(pos_comb)):
        
        fit_list.append(float(np.random.normal(dist_mean, stand_dev, 1) + 
                              add_fit_calc(fit_infl, opt_fit, opt_gen, pos_comb[num2])))


    mut_list = []

    for num3 in range(0, len(pos_comb)):
        
        num_mut  = 0
        
        for num4 in range(0, len(pos_comb[num3])):
            
            if pos_comb[num3][num4] != 1:
                num_mut += 1
    
        mut_list.append(num_mut)
        
    print(pos_comb)



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

    
        
    plt.plot(mut_list, fit_list,  'ro')
    print(fit_list, "list with final fitnesses")
    print(mut_list, "list with numbers of mutation")
    return fit_list
    return pos_comb
    

    

fit_calc(fit_infl, opt_fit, site_num, max_res_var, dist_mean, stand_dev, opt_gen(site_num, max_res_var))


