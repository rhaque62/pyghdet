from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys
missing_modules = []
INSTALL_ERROR   = False

# Try importing necessary modules to test
# and see if they are installed.
try:
    from Cython.Build import cythonize
except ImportError:
    missing_modules.append('cython')

try:
    import numpy
except ImportError:
    missing_modules.append('numpy')

try:
    import matplotlib
except ImportError:
    missing_modules.append('matplotlib')

try:
    import seaborn
except ImportError:
    missing_modules.append('seaborn')

try:
    import multiprocess
except ImportError:
    missing_modules.append('multiprocess')

try:
    import phyde
except ImportError:
    missing_modules.append('phyde')

try:
    import math
except ImportError:
    missing_modules.append('math')

try:
    import scipy
except ImportError:
    missing_modules.append('scipy')

try:
    import pandas
except ImportError:
    missing_modules.append('pandas')

try:
    import itertools
except ImportError:
    missing_modules.append('itertools')

try:
    import typing
except ImportError:
    missing_modules.append('typing')

if len(missing_modules) > 0:
    INSTALL_ERROR = True
    print("ERROR:")
    print("  You are missing the following required modules:")
    for m in missing_modules:
        print("  \t", m)
    print("\n")

if INSTALL_ERROR:
    print("ERROR:")
    print("  Unable to install phyde.")
    print("  Please see the documentation at https://github.com/rhaque62/pyghdet.git\n")
    sys.exit(-1)
else:
    setup(
        name="pyghdet",
        version="0.3.3",
        description="Global Hybridization detection using phylogenetic invariants",
        long_description="Global Hybridization detection using phylogenetic invariants",
        readme = "README.md",
        url="https://github.com/rhaque62/pyghdet",
        author="Rejuan Haque & Laura Kubatko",
        author_email="haque.62@osu.edu",
        packages=find_packages(),
        scripts=[
            'scripts/hdet_indiv.py',
            'scripts/hdet_species.py'
        ],
        license="GPLv3",
        classifiers=[
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Cython',
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
        ],
        zip_safe=False
    )