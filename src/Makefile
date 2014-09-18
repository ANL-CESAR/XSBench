#===============================================================================
# User Options
#===============================================================================

OPENACC     = no

VERIFY      = no
BINARY_DUMP = no
BINARY_READ = no

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = yes
PROFILE     = no
MPI         = no
PAPI        = yes
VEC_INFO    = no
BENCHMARK   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.c \
io.c \
CalculateXS.c \
GridInit.c \
XSutils.c \
Materials.c

obj = $(source:.c=.o)
ifeq ($(OPENACC),yes)
  cpp = $(source:.c=.cpp)
  cetus_output = cetus_output
  cetus_output_cpp = $(addprefix $(cetus_output)/,$(cpp))
  target = ./bin
endif

#===============================================================================
# Sets Flags
#===============================================================================

space := 
space +=
comma := ,
macro_base =

ifeq ($(OPENACC),yes)

  OMP ?= 1
  OPENARCRTLIB = ${HOME}/openarc/openarcrt
  OPENARCLIB = ${HOME}/openarc/lib
  XSBENCH_DIR = ${HOME}/XSBench/src
  CUDA_ARCH = sm_20
  CUDA_INCLUDE = -I/usr/local/cuda-6.0/include
  CUDA_LIB = -lcuda
  GNVCC_FLAGS = -arch $(CUDA_ARCH) -O3

  $(info Ignoring the following user options: \
  COMPILER, OPTIMIZE, DEBUG, PROFILE, MPI, PAPI, VEC_INFO, BENCHMARK)

  CC = g++
  NVCC = nvcc
  CFLAGS = -DOPENACC -I$(OPENARCRTLIB) -I$(XSBENCH_DIR) $(CUDA_INCLUDE) $(GNVCC_FLAGS)

  ifeq ($(OMP),1)
    LDFLAGS += -lgomp
    ACCRT = openaccrtomp_cuda
    CFLAGS += -Xcompiler -fopenmp
  else
    ACCRT = openaccrt_cuda
  endif

  LDFLAGS = -lm -L$(OPENARCRTLIB)/ -l$(ACCRT) $(OPENCL_LIB) $(CUDA_LIB)

  macro_base = OPENACC

else

  CFLAGS := -std=gnu99
  LDFLAGS = -lm

  # Regular gcc Compiler
  ifeq ($(COMPILER),gnu)
    CC = gcc
    CFLAGS += -fopenmp
  endif

  # Intel Compiler
  ifeq ($(COMPILER),intel)
    CC = icc
    CFLAGS += -openmp 
  endif

  # BG/Q gcc Cross-Compiler
  ifeq ($(MACHINE),bluegene)
    CC = mpicc
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

  # Compiler Vectorization (needs -O3 flag) information
  ifeq ($(VEC_INFO),yes)
    CFLAGS += -ftree-vectorizer-verbose=6
  endif

  # PAPI source (you may need to provide -I and -L pointing
  # to PAPI depending on your installation
  ifeq ($(PAPI),yes)
    source += papi.c
    CFLAGS += -DPAPI
    #CFLAGS += -I/soft/apps/packages/papi/papi-5.1.1/include
    #LDFLAGS += -L/soft/apps/packages/papi/papi-5.1.1/lib -lpapi
    LDFLAGS += -lpapi
  endif

  # MPI
  ifeq ($(MPI),yes)
    CC = mpicc
    CFLAGS += -DMPI
  endif

  # Adds outer 'benchmarking' loop to do multiple trials for
  # 1 < threads <= max_threads
  ifeq ($(BENCHMARK),yes)
    CFLAGS += -DBENCHMARK
  endif

endif

# Verification of results mode
ifeq ($(VERIFY),yes)
  CFLAGS += -DVERIFICATION
  macro_base += VERIFICATION
endif

# Binary dump for file I/O based initialization
ifeq ($(BINARY_DUMP),yes)
  CFLAGS += -DBINARY_DUMP
  macro_base += BINARY_DUMP
endif

# Binary read for file I/O based initialization
ifeq ($(BINARY_READ),yes)
  CFLAGS += -DBINARY_READ
  macro_base += BINARY_READ
endif

macro = $(subst $(space),$(comma),$(macro_base))

#===============================================================================
# Targets to Build
#===============================================================================
ifeq ($(OPENACC),yes)

program: build_cetus_output $(cetus_output_cpp) | target_dir
	cd $(cetus_output); $(NVCC) $(CFLAGS) -o ../$(target)/$(program) $(cpp) $(LDFLAGS); cp openarc_kernel.cu ../$(target); cp $(OPENARCRTLIB)/binBuilder_cuda ../$(target); cp $(OPENARCRTLIB)/Timer ../$(target); cd ../

target_dir:
	mkdir -p $(target)

build_cetus_output:
	java -classpath $(OPENARCLIB)/cetus.jar:$(OPENARCLIB)/antlr.jar openacc.exec.ACC2GPUDriver -macro=$(macro) -verbosity=2 $(source)

clean:
	rm -f *.o;
	if [ -d ./$(cetus_output) ]; then cd $(cetus_output); rm -f *.o *~ memcheck memcheck.c; cd ../; fi
	rm -rf $(target)

else

$(program): $(obj) XSbench_header.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) papi.c XSbench_header.h

run:
	./$(program)

bgqrun:
	qsub -t 10 -n 1 -O test XSBench

endif
