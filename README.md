# GHDet: Global Hybridization Detection using Phylogenetic Invariants
---------------------------------------------------------------------------------
<meta name="google-site-verification" content="reUTkch_kKqgtm2_N-6pEpu5RpqzLSdojBIX2klaie8" />

GHDet is a python package that can detect if there is any hybridization event given a set of taxa. 
The main module of this package is called ``pyghdet`` (**Py**\ thonic **G**\ lobal **H**\ ybridization **Det**\ ection).
Using this package, we can test whether there has been any hybridization event among a set of taxa, helping us decide whether or not to perform tree or network analysis.

# Installation
------------

Requirements:

-  Python 3.6 or above
-  Python Modules:

   -  Cython
   -  Numpy
   -  Matplotlib
   -  Seaborn
   -  Multiprocess
   -  math
   -  scipy
   -  phyde
   -  pandas
   -  itertools
   -  typing

-  C++ compiler

```bash
    # To install dependencies
    pip install cython numpy matplotlib seaborn multiprocess math scipy phyde pandas itertools typing

    # To install GHDet
    pip install pyghdet
```

# Usages
-------------

The package has two main modules: ``comb_species`` and ``comb_indiv``. The ``comb_species`` module performs global hypothesis test to detect if any of the species is a hybrid or not. Similarly, ``comb_indiv`` module performs global hypothesis test to detect if there is any hybrid individual in the data. The ``comb_species`` and ``comb_indiv`` modules has the following arguments:

```bash
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
```

# Examples
-----------------
## Detect if any species is hybrid:

```python
    ## No suspected hybrid:
    
      import pyghdet as ghd
      res = ghd.comb_speices("data.txt", "map.txt", "out", 16, 4, 50000)
      res.p_value
      res.detailed
      
    
    ## with suspected hybrid:
      import pyghdet as ghd
      res = ghd.comb_species("data.txt", "map.txt", "out", 16, 4, 50000, ['sp1', 'sp2'])
      res.p_value
      res.detailed
```

## Detect if any individual is hybrid:

```python
    ## No suspected hybrid:
    
      import pyghdet as ghd
      res = ghd.comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000)
      res.p_value
      res.detailed
      
    
    ## with suspected hybrid:
      import pyghdet as ghd
      res = ghd.comb_indiv("data.txt", "map.txt", "out", 16, 4, 50000, ['sp1', 'sp2'])
      res.p_value
      res.detailed
```

# Running GHDet from command line
------------------------------------

To facilitate analyses using the Python module, two scripts are provided to
conduct hybridization detection analyses directly from the command line:

- ``hdet_species.py``: runs a standard hybridization detection analysis on all the species in the data. The output is the p-value of the global test. The script also produce results on individual tests if the global test is significant. One can also provide a list of suspected hybrid, which only test if any of the suspected hybrid species is actually hybrid species or not. 
  evidence for hybridization.
- ``hdet_indiv.py``: tests if any of the individuals in the data is hybrid or not. It is possible to provide a list of suspected hybrid species, then the test will only detect individuals from the provided suspected hybrid species are hybrid or not.

# Examples

```bash
## species level

hdet_species.py -i data.txt -m map.txt -o out -n 16 -t 4 -s 50000 

## individual level

hdet_indiv.py -i data.txt -m map.txt -o out -n 16 -t 4 -s 50000 
```

# Performing Combination tests

The package can also be used to perform Combination tests. Three different combination tests have been implemented in the package. The tests are: [Cauchy Combination test (CCT)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7531765/), [MinP-CCT-MinP (MCM) test](https://www.nature.com/articles/s41598-022-07094-7), and [CCT-MinP-CCT-(CMC) test](https://www.nature.com/articles/s41598-022-07094-7). These tests can be performed as follows:

```python
import pyghdet as ghd

## CCT with different weights
ghd.cct([0.01,0.05,0.10, 0.53], [1,1,2,1])

## CCT with equal weights
ghd.cct([0.01, 0.10, 0.05, 0.54])

## CMC test
ghd.cmc([0.01,0.05,0.10,0.45])

## MCM test
ghd.mcm([0.01, 0.05, 0.40, 0.33])

```
