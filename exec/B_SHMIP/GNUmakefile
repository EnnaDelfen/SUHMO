# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
include ../../mk/Make.defs
ebase := Suhmo

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2
USE_PETSC=FALSE
MPI=FALSE
USE_MF=TRUE

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 

base_dir = .

#BASE_DIR = ../../

#EXEC_DIR = .

SRC_DIR = ../../src  

UTIL_DIR = ../../util 

src_dirs = $(SRC_DIR) $(UTIL_DIR)  

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

VERBOSE=-v
#BISICLES_MK = ../mk/
