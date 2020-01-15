#Makefile to compile the flextopo program
#Usage (linux systems):
#Adjust the compiler to your own compiler and set the paths
#Open terminal; go to the directory and type "make"
#To remove old compiled code, type "make clean" 
#Written by: R.C. Nijzink june 2014
#References: 
#Adjustments:  

#The name of the exetuble (the final program)
PROGRAM  = flexsimple


FC      := /usr/bin/gfortran #compiler name and directory
#FC      := /usr/bin/x86_64-w64-mingw32-gfortran #64bit


SRC_PATH = ./src		#path of source files
OBJ_PATH = ./tmp		#path .o and .mod files

NFDIR := lib/

SRC_PATH    := $(abspath $(SRC_PATH:~%=${HOME}%))
OBJ_PATH   := $(abspath $(OBJ_PATH:~%=${HOME}%))


# compiler flags and options
FCFLAGS = 
FLFLAGS =

#source objects
SRC = src/mo_init_random_seed.f90 src/mo_sort.f90 src/mo_readdata.f90 src/mo_objectives.f90 src/mo_signatures.f90 src/mo_eval_signatures.f90 src/mo_model.f90 src/mo_window.f90 src/mo_multi_obj.f90 src/mo_sample.f90 src/mo_resample.f90 src/driver.f90


#mod files
MSRC = *.mod 
    

all: $(PROGRAM) 
		
$(PROGRAM): $(SRC) 
	$(FC) -o $@ $^ ;rm $(MSRC) #$(OBJ_PATH)
#	$(FC) -o $@ $^ -I${NFDIR}/include -L${NFDIR}/lib  -static ;mv $(MSRC) $(OBJ_PATH)
%.o: %.f
	$(FC) -o $@ $<




.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD $(OBJ_PATH)/*.o $(OBJ_PATH)/*.mod $(OBJ_PATH)*.MOD $(OBJ_PATH)*.exe $(PROGRAM)   
