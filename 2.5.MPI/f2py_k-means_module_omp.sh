#!/bin/bash

#f2py -m add add.f
#f2py3 -c -m calc_dist calc_dist.f90
#f2py3 -c --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -m k_means_mod k-means_mod.f90
#f2py -c --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -m k_means_mod k-means_mod.f90
f2py -c --fcompiler=intelem --f90exec=mpiifort --f90flags='-fopenmp' -lgomp -m k_means_mod k-means_mod.f90 #\
    #-I/cm/shared/apps/intel/mpi/5.0.3.048/intel64/include \
    #-I/cm/shared/apps/intel/mpi/5.0.3.048/intel64/include \
    #-m k_means_mod k-means_mod.f90
