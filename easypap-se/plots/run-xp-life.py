#!/usr/bin/env python3

from graphTools import *
from expTools import *
import os

# Dictionnaire avec les options de compilations d'apres commande
options = {}
options["-k "] = ["life"]
options["-i "] = [100]
options["-v "] = ["seq", "omp"]
options["-wt "] = ["default", "sparse", "AVX2"]
options["-s "] = [6208]
options["-g "] = [4, 8, 16, 32]
options["-a "] = ["meta3x3"]

# Pour renseigner l'option '-of' il faut donner le chemin depuis le fichier easypap
options["-of "] = ["./plots/data/perf_data.csv"]


# Dictionnaire avec les options OMP
ompenv = {}
ompenv["OMP_NUM_THREADS="] = [1] + list(range(2, 9, 2))
ompenv["OMP_PLACES="] = ["cores", "threads"]

nbrun = 4
# Lancement des experiences
execute('./run ', ompenv, options, nbrun, verbose=True, easyPath=".")

# Lancement de la version seq avec le nombre de thread impose a 1
options["-v "] = ["seq"]
ompenv["OMP_NUM_THREADS="] = [1]
execute('../run', ompenv, options, nbrun, verbose=False, easyPath=".")
