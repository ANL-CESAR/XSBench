#===============================================================================
# User Options
#===============================================================================

COMPILER    ?= pgi
OPTIMIZE    ?= yes
DEBUG       ?= yes
PROFILE     ?= no
MPI         ?= no
OPENACC     ?= yes
ALIGNED     ?= yes
CUDA_ARCH   ?= 80

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

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm

# PGI Compiler
ifeq ($(COMPILER),pgi)
  CC = pgcc
  CFLAGS = -c99 -Minform=warn -mp -Minfo=mp
  ifeq ($(OPTIMIZE),yes)
    CFLAGS += -fast
  endif
  ifeq ($(VEC_INFO),yes)
    CFLAGS += -Minfo=vect
  endif
  # OpenACC (assumes NVIDA GPU)
  ifeq ($(OPENACC),yes)
    CFLAGS += -acc -Minfo=accel -gpu=cc$(CUDA_ARCH)
  endif
endif

# GNU Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -std=gnu99 -Wall -fopenmp
  # Optimization
  ifeq ($(OPTIMIZE),yes)
    CFLAGS += -O3
    # Compiler Vectorization (needs -O3 flag) information
    ifeq ($(VEC_INFO),yes)
      CFLAGS += -ftree-vectorizer-verbose=6
    endif
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

# Align workload for comparison purposes
ifeq ($(ALIGNED),yes)
  CFLAGS += -DALIGNED_WORK
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
