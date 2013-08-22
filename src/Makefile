#===============================================================================
# User Options
#===============================================================================

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
MPI       = no
PAPI      = yes
VEC_INFO  = no
VERIFY    = no
PAUSE     = no
BENCHMARK = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.c \
CalculateXS.c \
GridInit.c \
XSutils.c \
Materials.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  GCC = gcc
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  GCC = mpicc
endif

# Standard Flags
GCC_S_FLAGS := -fopenmp -std=gnu99

# Optional/Extra Flags
GCC_O_FLAGS =

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  GCC_O_FLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  GCC_O_FLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  GCC_O_FLAGS += -O3
endif

# Compiler Vectorization (needs -O3 flag) information
ifeq ($(VEC_INFO),yes)
  GCC_O_FLAGS += -ftree-vectorizer-verbose=6
endif

# PAPI source (you may need to provide -I and -L pointing
# to PAPI depending on your installation
ifeq ($(PAPI),yes)
  source += papi.c
  GCC_O_FLAGS += -DPAPI
  #GCC_S_FLAGS += -I/soft/apps/packages/papi/papi-5.1.1/include
  #LDFLAGS += -L/soft/apps/packages/papi/papi-5.1.1/lib -lpapi
  LDFLAGS += -lpapi
endif

# MPI
ifeq ($(MPI),yes)
  GCC = mpicc
  GCC_S_FLAGS += -DMPI
endif

# Verification of results mode
ifeq ($(VERIFY),yes)
  GCC_S_FLAGS += -DVERIFICATION
endif

# Injects artificial pause between XS lookups
# to represent time for particle tracking
ifeq ($(PAUSE),yes)
  GCC_S_FLAGS += -DPAUSE
endif

# Adds outer 'benchmarking' loop to do multiple trials for
# 1 < threads <= max_threads
ifeq ($(BENCHMARK),yes)
  GCC_S_FLAGS += -DBENCHMARK
endif
#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) do_flops.o XSbench_header.h
	$(GCC) $(GCC_S_FLAGS) $(GCC_O_FLAGS) do_flops.o $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(GCC) $(GCC_S_FLAGS) $(GCC_O_FLAGS) -c $< -o $@

do_flops.o: do_flops.c
	$(GCC) $(GCC_S_FLAGS) -c do_flops.c

clean:
	rm -rf XSBench XSBench.dSYM $(obj) do_flops.o

edit:
	vim -p $(source) do_flops.c papi.c XSbench_header.h

run:
	./XSBench

bgqrun:
	qsub -t 10 -n 1 -O test XSBench
