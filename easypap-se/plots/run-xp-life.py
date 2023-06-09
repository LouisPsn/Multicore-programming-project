#!/usr/bin/env python3

from graphTools import *
from expTools import *
import os

# Dictionnaire avec les options de compilations d'apres commande
options = {}
options["-k "] = ["life"]
options["-i "] = [40]
# options["-v "] = ["omp"]
# options["-wt "] = ["default", "AVX2"]
options["-s "] = [6208]
options["-a "] = ["meta3x3"]

# Pour renseigner l'option '-of' il faut donner le chemin depuis le fichier easypap
options["-of "] = ["./plots/data/perf_data.csv"]


# Dictionnaire avec les options OMP
ompenv = {}
# ompenv["OMP_NUM_THREADS="] = [1] + list(range(2, 9, 2))
# ompenv["OMP_PLACES="] = ["cores", "threads"]

ompenv["OMP_NUM_THREADS="] = [1]

nbrun = 1
# Lancement des experiences

# Lancement de la version seq avec le nombre de thread impose a 1
options["-v "] = ["seq"]
options["-wt "] = ["default"]

execute('./run', ompenv, options, nbrun, verbose=False, easyPath=".")

options ["-o "] = [""]
options["-v "] = ["ocl"]
execute('./run ', ompenv, options, nbrun, verbose=True, easyPath=".")