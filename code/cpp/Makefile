all: bin/chen1995_sample.out

CC=g++
# The -ftrapv makes the program abort if there's and integer overflow.
# This makes easies to debug overflow errors, but only works if the overflow
# is from signed numbers.
SFLAGS=-std=c++11 -Ilib -Wall -Wpedantic -Wconversion -ftrapv -O3 #-O0 #
HEAD=$(shell git describe --abbrev=12 --dirty --always --tags)
CFLAGS=$(SFLAGS) -DHBM_GIT_HEAD_AT_COMPILATION=\"$(HEAD)\"

## CPLEX OPTIONS
# Options used by Leonardo Moura for compiling with CPlex
CPLEX_EFLAGS  = -m64 -fPIC -fno-strict-aliasing -fexceptions -DIL_STD #-DNDEBUG 

# Folders where CPlex is stored
SYSTEM        = x86-64_linux
LIBFORMAT     = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio128/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio128/concert

# linking cplex
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEX_LFLAGS  = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -ldl
# includes for cplex
CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEX_IFLAGS  = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

CPLEX_EFLAGS += $(CPLEX_LFLAGS) $(CPLEX_IFLAGS)

bin:
	mkdir bin

bin/chen1995_sample.out: main/chen1995_sample.cpp lib/chen1995_model.hpp Makefile bin
	$(CC) $(CFLAGS) $< -o $@ $(CPLEX_EFLAGS)

clean:
	rm -rf bin

bin/extremum_planes_tests.out: main/extremum_planes_tests.cpp lib/extremum_planes.hpp lib/catch.hpp Makefile bin
	$(CC) $(CFLAGS) $< -o $@

run_ep_tests: bin/extremum_planes_tests.out Makefile
	./bin/extremum_planes_tests.out

