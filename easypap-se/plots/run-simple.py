#!/usr/bin/env python3

from graphTools import *
from expTools import *
import os

# options communes
options = {}
options["--label "] = ["final"]
options["-k "] = ["mandel"]
options["-i "] = [10]
options["-of "] = ["./plots/data/perf_data.csv"]
options["-s "] = [1024]

ompenv = {}
ompenv["OMP_NUM_THREADS="] = [1] + list(range(2, 13, 2))

nbrun = 15

# OMP_TILED GRAIN 8 16 32
options["-v "] = ["omp_tiled"]
options["-g "] = [8, 16, 32]
ompenv["OMP_SCHEDULE="] = ["dynamic"]

#execute('./run ', ompenv, options, nbrun, verbose=False, easyPath=".")

# OMP_LINE static
options["-v "] = ["omp_line"]
ompenv["OMP_SCHEDULE="] = ["static"]
options["-g "]=[1024]
execute('./run ', ompenv, options, nbrun, verbose=False, easyPath=".")

# OMP_LINE dynamic
ompenv["OMP_SCHEDULE="] = ["dynamic"]
options["-g "] = [1024]
execute('./run ', ompenv, options, nbrun, verbose=False, easyPath=".")
