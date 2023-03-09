# -*- Mode: Makefile -*- 
SUHMO_HOME  ?= ../../
CHOMBO_HOME ?= $(SUHMO_HOME)/Submodules/Chombo_git/lib
AMR_HOME    ?= $(SUHMO_HOME)

### This makefile produces an executable for each 
### name in the `ebase' variable
include ../../mk/Make.defs
ebase := Suhmo

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2
USE_PETSC=FALSE
MPI=TRUE
DEBUG=FALSE
OPT=FALSE

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 

base_dir = .



SRC_DIR = ../../src  

BISICLES_SRC_DIR = ../../util/bisiclesSrc

UTIL_DIR = ../../util 

src_dirs = $(SRC_DIR) $(UTIL_DIR)  $(BISICLES_SRC_DIR)

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

VERBOSE=-v
