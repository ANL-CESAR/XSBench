#Set USE_NVIDIA to 0 if OpenCL program runs on non-NVIDIA GPUs.
include make.header

#C compiler
CC = g++ 
NVCC = nvcc
CC_FLAGS_SERIAL = -O3 -Wall -ffloat-store 
CLIBS1 = -lm 
#TARGET is where the output binary is stored
TARGET = ./bin

#Set the output binary names for different compilers
BENCHMARK_SERIAL = $(TARGET)/xsbench_serial
BENCHMARK_OPENCL = $(TARGET)/xsbench_ARC_OPENCL
BENCHMARK_CUDA = $(TARGET)/xsbench_ARC_CUDA


#Choose the correct runtime library depending upon whether OMP is enabled or not.
ifeq ($(OMP),1)
CLIBS2_BASE = -lm -lgomp
xsbench_ARC_CUDA: CFLAGS2_BASE = -Xcompiler -fopenmp
xsbench_ARC_OPENCL: CFLAGS2_BASE = -fopenmp
ACCRT_BASE = openaccrtomp
else
CLIBS2_BASE = -lm 
CFLAGS2_BASE =
ACCRT_BASE = openaccrt
endif

#Check for the profiler mode
MODE ?= normal
ifeq ($(MODE), profile)
PREF = pf
else
PREF = 
endif

ifeq ($(ASPEN), aspenrt)
ASPENFILE = "aspenrt.cpp"
else 
ASPENFILE = 
endif

#Select underlying compiler-specific runtime library and compile flags
xsbench_ARC_OPENCL: ACCRT_BASE2 = $(ACCRT_BASE)_opencl$(PREF)
xsbench_ARC_OPENCL: CFLAGS2 = $(CFLAGS2_BASE) -I$(OPENARCLIB) $(OPENCL_INCLUDE) $(OPENCL_FLAGS)
xsbench_ARC_OPENCL: CLIBS2 = $(CLIBS2_BASE) -L$(OPENARCLIB)/ -l$(ACCRT_BASE2) $(OPENCL_LIB)

xsbench_ARC_CUDA: ACCRT_BASE2 = $(ACCRT_BASE)_cuda$(PREF)	
xsbench_ARC_CUDA: CFLAGS2 = $(CFLAGS2_BASE) -I$(OPENARCLIB) -I/home/alund/openarc/test/XSBench -I../../openarcrt/openacc.h $(CUDA_INCLUDE) $(GNVCC_FLAGS)
xsbench_ARC_CUDA: CLIBS2 = $(CLIBS2_BASE) -L$(OPENARCLIB)/ -l$(ACCRT_BASE2) $(OPENCL_LIB) $(CUDA_LIB)


ARC_CUDA: xsbench_ARC_CUDA
ARC_OPENCL: xsbench_ARC_OPENCL
all: all_CUDA
all_CUDA: xsbench_ARC_CUDA xsbench_serial
all_OPENCL: xsbench_ARC_OPENCL xsbench_serial
serial: xsbench_serial

#########
# Macro #
#########
#_N_ ?= 4096
OMP ?= 0
DEFSET_SERIAL ?= -D_N_=$(_N_) -DOMP=$(OMP)
DEFSET_ACC ?= -D_N_=$(_N_) -DOMP=$(OMP)

CETUS_OUTPUT_CPP = $(CETUS_OUTPUT)/CalculateXS.cpp $(CETUS_OUTPUT)/GridInit.cpp $(CETUS_OUTPUT)/io.cpp $(CETUS_OUTPUT)/Main.cpp $(CETUS_OUTPUT)/Materials.cpp $(CETUS_OUTPUT)/XSutils.cpp
CPP = CalculateXS.cpp GridInit.cpp io.cpp Main.cpp Materials.cpp XSutils.cpp

xsbench_serial: %.c XSbench_header.h
	$(CC) $(COMMONFLAGS) $(CC_FLAGS_SERIAL) $(DEFSET_SERIAL) -o $(BENCHMARK_SERIAL) %.c $(CLIBS1)

xsbench_ARC_CUDA: $(CETUS_OUTPUT_CPP)
	cd $(CETUS_OUTPUT); $(NVCC) $(CFLAGS2) $(DEFSET_ACC) -o ../$(BENCHMARK_CUDA) $(CPP) $(ASPENFILE) $(CLIBS2); cp openarc_kernel.cu ../$(TARGET); cp $(OPENARCLIB)/binBuilder_cuda ../$(TARGET); cp $(OPENARCLIB)/Timer ../$(TARGET); cd ../


xsbench_ARC_OPENCL: $(CETUS_OUTPUT_CPP)
	cd $(CETUS_OUTPUT); $(OPENCL_CC) $(CFLAGS2) $(DEFSET_ACC) -o ../$(BENCHMARK_OPENCL) $(CPP) $(ASPENFILE) $(CLIBS2); cp openarc_kernel.cl ../$(TARGET); cp $(OPENARCLIB)/binBuilder_opencl ../$(TARGET); cp $(OPENARCLIB)/Timer ../$(TARGET);  cd ../

clean:
	rm -f *.o;
	if [ -d ./$(CETUS_OUTPUT) ]; then cd $(CETUS_OUTPUT); rm -f *.o *~ memcheck aspenrt memcheck.c aspenrt.c; cd ..; fi
	if [ -d $(TARGET) ]; then rm -f $(BENCHMARK_SERIAL) $(BENCHMARK_CUDA) $(BENCHMARK_OPENCL) Timer; cd $(TARGET); rm -f openarc_kernel*.ptx openarc_kernel.cu openarc_kernel.cl binBuilder_cuda binBuilder_opencl Timer; fi

purge: clean
	if [ -d ./$(CETUS_OUTPUT) ]; then cd $(CETUS_OUTPUT); rm -f *.cpp *.h *.cl *.cu *.aspen *.c *.log aspenrt memcheck smallmodeltest; cd ..; fi
