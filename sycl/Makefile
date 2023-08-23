#===============================================================================
# User Options
#===============================================================================

COMPILER    ?= llvm
OPTIMIZE    ?= yes
DEBUG       ?= no
PROFILE     ?= no
MPI         ?= no
CHECK_MAX   ?= no
USE_CUDA    ?= yes
CUDA_ARCH   ?= sm_80

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
ifeq ($(COMPILER),codeplay)
  CC = compute++
  CFLAGS += -sycl -sycl-driver
  LDFLAGS += -lComputeCpp
endif

# LLVM Compiler
ifeq ($(COMPILER),llvm)
  CC = clang++
  CFLAGS += -fsycl
  ifeq ($(USE_CUDA),yes)
    CFLAGS += -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend --cuda-include-ptx=$(CUDA_ARCH)
  endif
  #LDFLAGS += -lOpenCL
endif

# IntelOneAPI Compiler
ifeq ($(COMPILER),oneapi)
  CC = icpx
  CFLAGS += -fsycl
endif

ifeq ($(COMPILER),opensycl)
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

# Disable or enable checking SYCL max mem alloc
ifeq ($(CHECK_MAX),yes)
	CFLAGS += -DCHECK_MAX_ALLOC
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
