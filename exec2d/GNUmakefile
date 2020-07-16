# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
include ../mk/Make.defs
ebase := driver 

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2
USE_PETSC=FALSE
MPI=TRUE
USE_MF=TRUE
USE_HDF=FALSE

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTimeDependent AMRTools BoxTools BaseTools 

base_dir = .

BASE_DIR = ../

EXEC_DIR = .

SRC_DIR = ../src  

UTIL_DIR = ../util 

src_dirs = $(SRC_DIR) $(UTIL_DIR)  

include $(CHOMBO_HOME)/mk/Make.example

VERBOSE=-v
#BISICLES_MK = ../mk/

ifneq ($(PYTHON_INC),)
CPPFLAGS+=-I$(PYTHON_INC) -DHAVE_PYTHON	
LIBFLAGS+=$(PYTHON_LIBS)
endif

nolink: $(shell $(CHOMBO_HOME)/mk/reverse $(LibNames)) $(objects)  
	@echo BISICLES Objects built

doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""
