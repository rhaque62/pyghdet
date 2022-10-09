#!/usr/bin/env python
# -*- coding: utf-8 -*-

# indiv_level.py
# Rejuan
"""
<<hdet_indiv.py>>

Run a global hybridization detection analysis or test for any hybrid individuals.

Arguments
---------
For more details on script arguments, type: run_hyde.py -h
 
    - infile         <string> : name of the DNA sequence data file.
    - mapfile        <string> : name of the taxon map file.
    - outgroup       <string> : name of the outgroup.
    - nindiv            <int> : number of sampled individuals.
    - ntaxa             <int> : number of sampled taxa/populations.
    - nsites            <int> : number of sampled sites.
    - sus_hyb         <string>: list of suspected hybrid species.
    - alpha            <float>: intended level of significance.
    - ignore_amb_sites <flag> : ignore missing/ambiguous sites.
        
        
Output
------
Return the p-value that test the null hyphothesis that there is no hybrid individual in the data.
"""
import pyghdet
import argparse

    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Options for individual_hyde.py",
                                     add_help=True)

    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--infile', action="store", type=str, required=True,
                          metavar='\b', help="name of the data input file")
    required.add_argument('-m', '--map', action="store", type=str, required=True,
                          metavar='\b', help="map of individuals to taxa")
    required.add_argument('-o', '--outgroup', action="store", type=str, required=True,
                          metavar='\b', help="name of the outgroup (only one accepted)")
    required.add_argument('-n', '--num_ind', action="store", type=int, required=True,
                          metavar='\b', help="number of individuals in data matrix")
    required.add_argument('-t', '--num_taxa', action="store", type=int, required=True,
                          metavar='\b', help="number of taxa (species, OTUs)")
    required.add_argument('-s', '--num_sites', action="store", type=int, required=True,
                          metavar='\b', help="number of sites in the data matrix")

    additional = parser.add_argument_group("additional arguments")
    additional.add_argument('-q', '--quiet', action="store_true",
                            help="supress printing to stdout")
    additional.add_argument('-sus_hyb', '--sus_hyb', action="store", type=str,
                            help="comma seperated list of suspected hybrids")
    additional.add_argument('-a', '--alpha', action="store", type=float,
                            help="Chosen level of significance")
    additional.add_argument('--ignore_amb_sites', action="store_true",
                            help="ignore missing/ambiguous sites")

    args             = parser.parse_args()
    infile           = args.infile
    mapfile          = args.map
    outgroup         = args.outgroup
    nind             = args.num_ind
    ntaxa            = args.num_taxa
    nsites           = args.num_sites
    quiet            = args.quiet
    sus_hyb          = args.sus_hyb
    alpha            = args.alpha
    ignore_amb_sites = args.ignore_amb_sites

    
    if not quiet: print("\nRunning hdet_indiv.py")

    if sus_hyb != None:
        sus_hyb = list(sus_hyb.split(","))
    if alpha == None:
        alpha = 0.05
    
    pyghdet.comb_indiv(infile, mapfile, outgroup, nind, ntaxa, nsites, sus_hyb, alpha, ignore_amb_sites)
    