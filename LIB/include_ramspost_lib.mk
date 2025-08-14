#Makefile include include.mk
############################## Change Log ##################################
# 1.0.0.2
#
# 000908 MJB include.mk-mrc ##
#            Added MAKE environment varaible.
#            Added free format option to F_OPTS for some platforms. ##
# 000907 MJB include.mk-mrc ##
#            Changed the defualts to no NCAR Graphics and no parallel.
#            Also commented out the machine specifics to force the user to
#            select the appropriate machine for them. ##
# 000823 MJB include.mk-mrc ##
#            New - defines all make environment varaibles and is included
#            in all make files. ##
#
############################################################################

# Define make (gnu make works best).

MAKE=/usr/bin/make

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
#F_OPTS=-O2  -i8
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

# For IBM,HP,SGI,ALPHA use these:
#ARCHIVE=ar rs
# For SUN,CONVEX, LINUX
ARCHIVE=ar r
