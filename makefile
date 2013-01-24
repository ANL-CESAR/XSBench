program = XSBench

source = \
CalculateXS.c \
GridInit.c \
Main.c \
Materials.c \
XSutils.c

#===============================================================================
# User Options
#===============================================================================

COMPILER = gnu
OPTIMIZE = yes
DEBUG    = no
PROFILE  = no
PAPI     = no

#===============================================================================
# Compiler Flags and Options
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  GCC = gcc
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  GCC = /bgsys/drivers/toolchain/V1R1M2/gnu-linux/bin/powerpc64-bgq-linux-gcc
endif 

# Standard Flags
GCCFLAGS := -fopenmp -Wall -std=c99 -lm
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  GCCFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  GCCFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  GCCFLAGS += -O3
endif

# PAPI source 
ifeq ($(PAPI),yes)
  source += papi.c
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(source) do_flops.o XSbench_header.h
	$(GCC) $(GCCFLAGS) do_flops.o $(source) -o $@ $(LDFLAGS)
do_flops.o: do_flops.c
	$(GCC) -fopenmp -Wall -std=c99 -c do_flops.c -lm

clean:
	rm -rf XSBench XSBench.dSYM do_flops.o
edit:
	vim -p $(source) do_flops.c papi.c
