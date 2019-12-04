#===============================================================================
# User Options
#===============================================================================

# Compiler can be set below, or via environment variable
ifeq ($(CC),)
CC    = gcc
endif
OPTIMIZE    = yes
OPENMP      = yes
DEBUG       = no
PROFILE     = no
MPI         = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.c \
io.c \
Simulation.c \
GridInit.c \
XSutils.c \
Materials.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Make sure we didn't pick up lower case cc
ifeq (cc,$(CC))
CC    = gcc
endif

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifneq (,$(findstring gcc,$(CC)))
  CFLAGS += -flto
  ifeq ($(OPENMP),yes)
	CFLAGS += -fopenmp -DOPENMP
  endif
endif

# LLVM Compiler
ifneq (,$(findstring clang,$(CC)))
  CFLAGS += -flto
  ifeq ($(OPENMP),yes)
	CFLAGS += -fopenmp -DOPENMP
  endif
endif

# Intel Compiler
ifneq (,$(findstring intel,$(CC)))
  CFLAGS += -ipo
  ifeq ($(OPENMP),yes)
	CFLAGS += -fopenmp -DOPENMP
  endif
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# MPI
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) XSbench_header.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c XSbench_header.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) XSbench_header.h

run:
	./$(program)
