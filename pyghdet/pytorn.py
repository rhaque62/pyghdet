import numpy as np
import math
from scipy.stats import cauchy
import phyde as hd
import pandas as pd
from itertools import combinations
from typing import NamedTuple


## A class to keep the detailed results
class result_det(NamedTuple):
    """
    A class to hold the result of the global hybrid detection test.
    
    Example:
    .. code:: py
      result_det(global_pv, sig_res)  
      
    """
    
    p_value : float
    detailed : None
    __slot__ = ()
    def __repr__(self):
        return f"\np_value: {self.p_value}\n\ndetailed:\n{self.detailed}"


## A class to keep the p-value of the global test
class result_pv(NamedTuple):
    """
    A class to hold the p_value of the global hybrid detection test.
    
    Example:
    .. code:: py
      result_pv(global_pv)  
      
    """
    p_value : float
    __slot__ = ()
    def __repr__(self):
        return f"\np_value:{self.p_value}"


## A function to get the combinition of species
def spcomb(species_list, sus_species = None):
    """
    A function that takes the list of species (without the outgroup) and the 
    list of suspected hybrid species and return all possible individual test
    setup with two parents and a hybrid species.
    
    Example:
    .. code:: py
      spcomb(['sp1', 'sp2', 'sp3', 'sp4'], ['sp2']) 
      
    """
    if set(sus_species).issubset(species_list):
        unq_species_new = list(set(species_list).difference(sus_species))
    else:
       return f"Error:The provided suspected hybrid/s {sus_species} is/are not in the list of species in the data!"
    comb = []
   
    for sh in range(len(sus_species)):
        hyb = sus_species[sh]
        nhy = list(set(sus_species).difference([hyb]))
        new_unq_species = unq_species_new + nhy
        com = combinations(new_unq_species,2)
        for item in com:
            c_item = list(item)
            c_item.insert(1,hyb)
            comb.append(c_item)
    return(comb)



## cauchy combination test codes
def cct(pvals, weights = None):
    """
    A function to perform the Cauchy combination test. It takes a list of
    p-values and a list of weights and return the global p-value.
    
    Example:
    .. code:: py
      import pyghdet as ghd
      ghd.cct([0.01,0.05,0.55, 0.99, 0.02])
      
    """
    
    pv_arr = np.array(pvals)
    
    ## check if there is any non-numeric values in the p-vals
    if not all(type(e) in (int, float) for e in pvals):
        return "Warning: The individual tests produced p-values containing non-numeric character! Failed to test the global null hypothesis"
    
    ## check if all the p-vals are between 0 and 1
    
    if sum(pv_arr<0) > 0 or sum(pv_arr>1) > 0:
        return "Warning: All the individual p-values must be between 0 and 1! Failed to test the global null hypothesis"
    
    ## check if there are p-values that are exactly 0 or 1
    
    is_zero = sum(pv_arr==0)
    is_one = sum(pv_arr==1)
    
    if(is_zero>0 and is_one>0):
        return "Error: cannot have both 0 and 1 p-values!"
    elif(is_zero>0):
        print("Warning: there are p-values that are exactly zero")
        return 0
    elif(is_one>0):
        print("Warning: there are p-values that are exactly one")
        return 1 
    

    ## check the weights
    if weights==None:
        weights = [1/len(pvals)]*len(pvals)
    else:
        if len(pvals) != len(weights):
            return "Error: weights and pvlaues should be same length!"
        else:
            w_ar = np.array(weights)
            if sum(w_ar<0)>0:
                return "Error: All the weights must be positive!"
            else:
                weights = [ w/sum(weights) for w in weights]
                
    
    ## create two groups for p-values: very small and large
        
    small_p = []
    small_w = []
    large_p = []
    large_w = []
    
    
    i = 0
    for item in pvals:
        if item < 1e-16:
            small_p.append(item)
            small_w.append(weights[i])
            i = i+1
        else:
            large_p.append(item)
            large_w.append(weights[i])
    
    if(len(small_p)==0):
        cct_stat = []
        j = 0
        for item in large_p:
            cct_stat.append(large_w[j]*math.tan((0.5 - item)*math.pi))
            j = j+1
        cct_stat = sum(cct_stat)
    else:
        cct_small = []
        cct_large = []
        k = 0
        for item in small_p:
            cct_small.append((weights[k]/item)/math.pi)
            k = k+1
        cct_small = sum(cct_small)
        
        l = 0
        for item in large_p:
            cct_large.append(large_w[l]*math.tan((0.5 - item)*math.pi))
            l = l+1
        cct_large = sum(cct_large)
        
        cct_stat = cct_small + cct_large
    
    ## calculate the p-value for the global test
        
    if(cct_stat>1e+15):
        pval = (1/cct_stat)/math.pi
    else:
        pval = 1 - cauchy.cdf(cct_stat)
    
    return(float(pval))



## code for MCM test
def mcm(pval):
    """
    A function to perform the MCM test. It takes a list of 
    p-values and return the global p-value.
    
    Example:
    .. code:: py
      import pyghdet as ghd
      ghd.mcm([0.01,0.05,0.55, 0.99, 0.02])
      
    """
    p_min = min(1,len(pval)*min(pval))
    p_mcm = min(1, 2*min(cct(pval),p_min))
    return p_mcm


## code for CMC test
def cmc(pval):
    """
    A function to perform the CMC test. It takes a list of 
    p-values and return the global p-value.
    
    Example:
    .. code:: py
      import pyghdet as ghd
      ghd.cmc([0.01,0.05,0.55, 0.99, 0.02])
      
    """
    p_min = min(1,len(pval)*min(pval))
    p_cmc = cct([cct(pval), p_min])
    return p_cmc



## combination test for individuals
def comb_indiv(infile, mapfile, outgroup, nindiv, ntaxa, nsite, sus_hyb = None, alpha = 0.05, remove_amb_site = False):
    
    """
    Main method for testing the global null hypothesis: there is no hybrid 
    individual in the data. It is also possible to provide a set of suspected 
    hybrid species
   
    
    Arguments
    ---------
 
        - infile         <string> : name of the DNA sequence data file.
        - mapfile        <string> : name of the taxon map file.
        - outgroup       <string> : name of the outgroup.
        - nindiv            <int> : number of sampled individuals.
        - ntaxa             <int> : number of sampled taxa/populations.
        - nsites            <int> : number of sampled sites.
        - sus_hyb         <string>: list of suspected hybrid species.
        - alpha            <float>: intended level of significance.
        - ignore_amb_sites <flag> : ignore missing/ambiguous sites.
        
        
    Example(No suspected hybrid):
    .. code:: py
      import pyghdet as ghd
      res = ghd.comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000)
      
    
    Example(with suspected hybrid):
    .. code:: py
      import pyghdet as ghd
      res = ghd.comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000, ['sp1'])
    """
    
    ## creating hyde data file and reading the map file
    if remove_amb_site:
        dat = hd.HydeData(infile, mapfile, outgroup, nindiv, ntaxa, nsite, True)
    else:
        dat = hd.HydeData(infile, mapfile, outgroup, nindiv, ntaxa, nsite)
    
    mapf = pd.read_csv(mapfile, delimiter='\t', header=None)
    
    
    ## get the name of all the species as list
    species_all = list(mapf.iloc[:, 1])

    ## remove the outgroup from the list of species

    species_all = [ sp for sp in species_all if sp!= outgroup ]

    ## select all the unique species

    unq_species =[]

    for sp in species_all:
        if sp not in unq_species:
            unq_species.append(sp)
    
    if sus_hyb == None:
        comb = spcomb(unq_species, unq_species)
    else:
        if set(sus_hyb).issubset(unq_species):
            comb = spcomb(unq_species, sus_hyb)
        else:
           return f"Error:The provided suspected hybrid/s {sus_hyb} is/are not in the list of species in the data!"
    
    p_val=[]
    Z_score=[]
    gamma = []
    parent1=[]
    hybrid=[]
    parent2=[]
    
    for item in comb:
        p1 = item[0]
        h = item[1]
        p2 = item[2]
        
        res1 = dat.test_individuals(p1, h, p2)
        
        for item in res1:
            res_each = res1[item]
            p_val.append(res_each["Pvalue"])
            Z_score.append(res_each["Zscore"])
            gamma.append(res_each["Gamma"])
            parent1.append(p1)
            parent2.append(p2)
            hybrid.append(item) 
            

    
    result = pd.DataFrame(list(zip(parent1, hybrid, parent2, gamma, Z_score, p_val)),
                          columns=["Parent1", "Hybrid", "Parent2","Gamma","Z_score", 
                                   "P_value"])

    
    sig_res = result[result["P_value"]<alpha]
    sig_res = sig_res[sig_res["Gamma"] <= 1]
    sig_res = sig_res[sig_res["Gamma"]>= 0]
    
    ## making the p-value ready for cct
    pvs = []
    
    for p in p_val:
        if p==0:
            pvs.append(1e-17)
        else:
            if p==1:
                pvs.append(0.99)
            else:
                pvs.append(p)
    
    
    ## running the mcm test
    global_pv = mcm(pvs)
    
    
    ## returning the significant results if global null is rejected
    if global_pv <= alpha:
        return result_det(global_pv, sig_res) 
    else:
        return result_pv(global_pv)





## combination test for species
def comb_species(infile, mapfile, outgroup, nindiv, ntaxa, nsite, sus_hyb = None, alpha = 0.05, remove_amb_site = False):
    
    """
    Main method for testing the global null hypothesis: there is no hybrid 
    species in the data. It is also possible to provide a set of suspected 
    hybrid species
   
    
    Arguments
    ---------
 
        - infile         <string> : name of the DNA sequence data file.
        - mapfile        <string> : name of the taxon map file.
        - outgroup       <string> : name of the outgroup.
        - nindiv            <int> : number of sampled individuals.
        - ntaxa             <int> : number of sampled taxa/populations.
        - nsites            <int> : number of sampled sites.
        - sus_hyb         <string>: list of suspected hybrid species.
        - alpha            <float>: intended level of significance.
        - ignore_amb_sites <flag> : ignore missing/ambiguous sites.
        
        
    Example(No suspected hybrid):
    .. code:: py
      import pyghdet as ghd
      res = ghd.comb_species("data.txt", "map.txt", "out", 16, 4, 50000)
      
    
    Example(with suspected hybrid):
    .. code:: py
      import pyghdet as ghd
      res = ghd.comb_species("data.txt", "map.txt", "out", 16, 4, 50000, ['sp1'])
    """
    
    ## creating hyde data file and reading the map file
    if remove_amb_site:
        dat = hd.HydeData(infile, mapfile, outgroup, nindiv, ntaxa, nsite, True)
    else:
        dat = hd.HydeData(infile, mapfile, outgroup, nindiv, ntaxa, nsite)
    
    
    mapf = pd.read_csv(mapfile, delimiter='\t', header=None)
    
    
    ## get the name of all the species as list
    species_all = list(mapf.iloc[:, 1])
    

    ## remove the outgroup from the list of species
    species_all = [ sp for sp in species_all if sp!= outgroup ]


    ## select all the unique species
    unq_species =[]

    for sp in species_all:
        if sp not in unq_species:
            unq_species.append(sp)
    
    
    if sus_hyb == None:
        comb = spcomb(unq_species, unq_species)
    else:
        if set(sus_hyb).issubset(unq_species):
            comb = spcomb(unq_species, sus_hyb)
        else:
           return f"Error:The provided suspected hybrid/s {sus_hyb} is/are not in the list of species in the data!"
        
    
    p_val=[]
    Z_score=[]
    gamma = []
    parent1=[]
    hybrid=[]
    parent2=[]
    
    
    
    for item in comb:
        p1 = item[0]
        h = item[1]
        p2 = item[2]
        
        res1 = dat.test_triple(p1, h, p2)
        parent1.append(p1)
        parent2.append(p2)
        hybrid.append(h)
        p_val.append(res1["Pvalue"])
        Z_score.append(res1["Zscore"])
        gamma.append(res1["Gamma"])
    

    result = pd.DataFrame(list(zip(parent1, hybrid, parent2, gamma, Z_score, p_val)),
                          columns=["Parent1", "Hybrid", "Parent2","Gamma", "Z_score", 
                                   "P_value"])
    
    sig_res = result[result["P_value"]< alpha]
    
    ## making the p-value ready for cct
    pvs = []
    
    for p in p_val:
        if p==0:
            pvs.append(1e-17)
        else:
            if p==1:
                pvs.append(0.99)
            else:
                pvs.append(p)
    
    
    ## running the mcm test
    global_pv = mcm(pvs)
    
    
    ## returning the significant results if global null is rejected
    
    if global_pv <= alpha:
        return result_det(global_pv, sig_res)
    else:
        return result_pv(global_pv)
