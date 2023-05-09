#!/usr/bin/env python3
from graphTools import *
from expTools import *
import os

easypap_options = {}
easypap_options["--kernel "] = ["scrollup"]
easypap_options["--variant "] = ["seq", "ji"]
easypap_options["--output-file "] = ["scrollup.csv"]

omp_icv = {}  # OpenMP Internal Control Variables
omp_icv["OMP_NUM_THREADS="] = [1]
omp_icv["OMP_SCHEDULE="] = ["static"]

for size in [2**i for i in range(5, 13)]:
    easypap_options["--iterations "] = [4096*4096//(size*size)]
    easypap_options["--size "] = [size]
    execute('./run', omp_icv, easypap_options, nbrun=5)
