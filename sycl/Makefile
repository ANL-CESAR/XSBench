#===============================================================================
# User Options
#===============================================================================

CC          ?= llvm
OPTIMIZE    ?= yes
DEBUG       ?= no
PROFILE     ?= no
MPI         ?= no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.cpp \
io.cpp \
Simulation.cpp \
GridInit.cpp \
XSutils.cpp \
Materials.cpp

obj = $(source:.cpp=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=c++17 -Wall

# Linker Flags
LDFLAGS = -lm

# Codeplay Compiler
ifeq ($(CC),codeplay)
  CC = compute++
  CFLAGS += -sycl -sycl-driver
  LDFLAGS += -lComputeCpp
endif

# LLVM Compiler
ifeq ($(CC),llvm)
  CC = clang++
  CFLAGS += -fsycl
  LDFLAGS += -lOpenCL
endif

# IntelOneAPI Compiler
ifeq ($(CC),oneapi)
  CC = icpx
  CFLAGS += -fsycl
endif

ifeq ($(CC),opensycl)
  CC = syclcc
  CFLAGS += -fsycl
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

%.o: %.cpp XSbench_header.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) XSbench_header.h

run:
	./$(program)
