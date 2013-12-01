#/*
#* Ising Model Monte Carlo Simulation via Cluster Algorithm
#* Copyright (C) 2013  Leopold Kellers
#*
#* This program is free software; you can redistribute it and/or
#* modify it under the terms of the GNU Lesser General Public
#* License as published by the Free Software Foundation; either
#* version 2.1 of the License, or (at your option) any later version.
#*
#* This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#* Lesser General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with this library; if not, write to the Free Software
#* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#*/

shell = /bin/sh
IDIR =.#include files are in current directory
GCC=gcc
ICC=icc
CC=$(GCC)
GCCFLAGS = -lm -pedantic -Wall -Wextra #flags specific to the gnu compiler
ICCFLAGS = -w2 #flags specific to the intel compiler
#general flags
CFLAGS=-I$(IDIR) -std=c99 -DD=$(D) -DRNG=$(RNG) $(RNGFLAGS) 
USEALGORITHMFLAGS = -DUSE_WOLFF #-DUSE_SWENDSEN_WANG # any combination of these
GCCOPTIMIZATIONFLAGS = -march=native -flto -ffast-math
ICCOPTIMIZATIONFLAGS = -ipo -xHost -static
COPTIMIZATIONFLAGS = -DNDEBUG -O3
CDEBUGFLAGS = -O0
GCCDEBUGFLAGS = -g
ICCDEBUGFLAGS = -debug
PROGNAME=ising_cluster_monte_carlo
TESTNAME=test_ising_cluster_monte_carlo

#Select the number of dimensions, the simulation should run in
#Must be a positive integer number obviously
D = 3

#add header files here:
_DEPS = monte_carlo.h data_analysis.h $(_RNG_DEPS) common.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
#the build time is only a fraction of a second, so we can save the hassle of
#building the object files separately
LIB_SRC = monte_carlo.c data_analysis.c common.c $(RNG_SRC)

#Choose the random number generator - can be 1 or 2 (beware of trailing whitespace!)
#1 for Double precision SIMD-oriented Fast Mersenne Twister
#2 for WELL Random number generator
RNG = 1
_RNG_DEPS_DSFMT = dsfmt/dSFMT.h
_RNG_SRC_DSFMT  = dsfmt/dSFMT.c
_RNG_DEPS_WELL  = well/WELL19937a.h
_RNG_SRC_WELL   = well/WELL19937a.c
ifeq ($(RNG),1)
	_RNG_DEPS=$(_RNG_DEPS_DSFMT)
	RNG_SRC=$(_RNG_SRC_DSFMT)
	RNGFLAGS=-DDSFMT_MEXP=19937
else 
	ifeq ($(RNG), 2)
		_RNG_DEPS = $(_RNG_DEPS_WELL)
		RNG_SRC  = $(_RNG_SRC_WELL)
		RNGFLAGS=-DTEMPERING
	else
		RNG_DEPS =
		RNG_SRC =
	endif
endif

ifeq ($(CC),$(GCC))
CFLAGS += $(GCCFLAGS)
COPTIMIZATIONFLAGS += $(GCCOPTIMIZATIONFLAGS)
CDEBUGFLAGS += $(GCCDEBUGFLAGS)
else
	ifeq ($(CC),$(ICC))
	CFLAGS += $(ICCFLAGS)
	COPTIMIZATIONFLAGS += $(ICCOPTIMIZATIONFLAGS)
	CDEBUGFLAGS += $(ICCDEBUGFLAGS)
	else
	#just hope that this unknown compiler supports all of the flags
	#which are still there as default flags... 
	endif
endif

default: all

all: clean
all: CFLAGS += $(COPTIMIZATIONFLAGS) $(USEALGORITHMFLAGS)
all: $(PROGNAME)

debug: clean
debug: CFLAGS += $(CDEBUGFLAGS)  $(USEALGORITHMFLAGS)
debug: $(PROGNAME)

profile: clean
profile: CFLAGS += -p $(COPTIMIZATIONFLAGS) $(USEALGORITHMFLAGS)
profile: $(PROGNAME)

#currently unused
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $(CFLAGS) $<

$(PROGNAME): MAINFILE = main.c
$(PROGNAME): $(MAINFILE) $(LIB_SRC)
	$(CC) -o $@ $(CFLAGS) $(MAINFILE) $(LIB_SRC)

$(TESTNAME): MAINFILE = test.c
$(TESTNAME): D = 3
$(TESTNAME): USEALGORITHMFLAGS = -DUSE_WOLFF -DUSE_SWENDSEN_WANG
$(TESTNAME): $(TESTFILE) $(LIB_SRC)
	$(CC) -o $@ $(CFLAGS) $(USEALGORITHMFLAGS) $(MAINFILE) $(LIB_SRC)


.PHONY: test

test: $(TESTNAME) 
	./$(TESTNAME) seed > /tmp/data
	gnuplot -persist -e "file=\"/tmp/data\"; L=32" fit_xi_comparison

.PHONY: run

run:
	./$(PROGNAME) -s seed

.PHONY: clean

clean:
		rm -f $(ODIR)/*.o *~ $(PROGNAME) $(IDIR)/*~ $(TESTNAME)

