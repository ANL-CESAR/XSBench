#===============================================================================
# User Options
#===============================================================================

KOKKOS_PATH = /gpfs/alpine/csc452/proj-shared/perf_port_builds/kokkos
CUDA_PATH   = /sw/summit/cuda/11.5.2
COMPILER    = nvidia
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
MPI         = no
KOKKOS_DEVICES = OpenMP,Cuda
KOKKOS_ARCH = Power9,Volta70
SM_VERSION = 70

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
CXXFLAGS := -Wall

# Linker Flags
LDFLAGS = -lm

# NVIDIA Compiler
ifeq ($(COMPILER),nvidia)
  CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
  CXXFLAGS += -Xcompiler -Wall -Xcompiler -O3 -arch=sm_$(SM_VERSION)
endif

# Clang Compiler
ifeq ($(COMPILER),llvm)
  CXX = clang++
  CXXFLAGS += -flto -fopenmp -DOPENMP
endif

# GCC Compiler
ifeq ($(COMPILER),gnu)
  CXX = g++
  CXXFLAGS += -flto -fopenmp -DOPENMP
endif


# Debug Flags
ifeq ($(DEBUG),yes)
  CXXFLAGS += -g
  LDFLAGS  += -g
ifeq ($(COMPILER),nvidia)
  CXXFLAGS += -G
  LDFLAGS  += -G
endif
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CXXFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CXXFLAGS += -O3
endif

# MPI
ifeq ($(MPI),yes)
  CXX = mpicxx
  CXXFLAGS += -DMPI
endif

#===============================================================================
# Targets to Build
#===============================================================================

default: $(program)

include $(KOKKOS_PATH)/Makefile.kokkos

$(program): $(obj) XSbench_header.hpp Makefile $(KOKKOS_LINK_DEPENDS)
	$(CXX) $(CXXFLAGS) $(KOKKOS_LDFLAGS) $(obj) $(KOKKOS_LIBS) -o $@ $(LDFLAGS)

%.o: %.cpp XSbench_header.hpp Makefile $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) XSbench_header.hpp

run:
	./$(program)
