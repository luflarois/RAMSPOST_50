# RAMSPOST root directory.

RAMSPOST=./src


# Source directories.

UTILS_LIB=./LIB
UTILS_INCS=./include

# MRC libraries.

LIBUTILS=$(UTILS_LIB)/libutils-ramspost.a


# Activate appropriate parts below, comment out others.
# Machine-dependent options.

CMACH=PC_LINUX1
#-----------------  INTEL COMPILER ---------------
#F_COMP=/opt/intel/fce/9.1.041/bin/ifort
#C_COMP=/opt/intel/cce/9.1.047/bin/icc
#LOADER=/opt/intel/fce/9.1.041/bin/ifort

#F_OPTS=-O0 -CB -traceback -FR #-fpe0 
#C_OPTS=-O0 -CB -traceback -FR #-fpe0 
#LOADER_OPTS=-O0 -CB -traceback -FR #-fpe0

#-----------------  LINUX Portland Group pgf90/pgcc ---------------
#F_COMP=pgf90
#F_OPTS=-O2 -i8
#C_COMP=pgcc
#C_OPTS=-O2 -DUNDERSCORE -DLITTLE
#LOADER=pgf90
#LOADER_OPTS=-v

#-----------------  GFORTRAN COMPILER ---------------
F_COMP=gfortran
F_OPTS=-O2 -g -fbacktrace -ffree-form -ffree-line-length-none -fno-range-check -fallow-argument-mismatch
C_COMP=gcc
C_OPTS=-O2 -g -fpermissive -std=gnu89
LOADER=gfortran
LOADER_OPTS=-O2 -g -fbacktrace

#-----------------------------------------------------
LIBS=
# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
