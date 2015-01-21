include make.header

#C compiler
CC = g++ 
NVCC = nvcc
CC_FLAGS_SERIAL = -O3 -Wall -ffloat-store 
CLIBS1 = -lm 
TARGET = ./bin

#Set the output binary names for different compilers
BENCHMARK_SERIAL = $(TARGET)/xsbench_serial
BENCHMARK_OPENCL = $(TARGET)/xsbench_ARC_OPENCL
BENCHMARK_CUDA = $(TARGET)/xsbench_ARC_CUDA

#Choose the correct runtime library depending upon whether OMP is enabled or not.
OMP ?= 0
ifeq ($(OMP),1)
CLIBS2_BASE = -lm -lgomp
ACCRT_BASE = openaccrtomp
xsbench_ARC_CUDA: CFLAGS2_BASE = -Xcompiler -fopenmp
xsbench_ARC_OPENCL: CFLAGS2_BASE = -fopenmp
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

#Select underlying compiler-specific runtime library and compile flags
xsbench_ARC_OPENCL: ACCRT_BASE2 = $(ACCRT_BASE)_opencl$(PREF)
xsbench_ARC_OPENCL: CFLAGS2 = $(CFLAGS2_BASE) -I$(OPENARCLIB) -I$(XSBENCH_DIR) $(OPENCL_INCLUDE) $(OPENCL_FLAGS)
xsbench_ARC_OPENCL: CLIBS2 = $(CLIBS2_BASE) -L$(OPENARCLIB)/ -l$(ACCRT_BASE2) $(OPENCL_LIB)

xsbench_ARC_CUDA: ACCRT_BASE2 = $(ACCRT_BASE)_cuda$(PREF)	
xsbench_ARC_CUDA: CFLAGS2 = $(CFLAGS2_BASE) -I$(OPENARCLIB) -I$(XSBENCH_DIR) $(CUDA_INCLUDE) $(GNVCC_FLAGS)
xsbench_ARC_CUDA: CLIBS2 = $(CLIBS2_BASE) -L$(OPENARCLIB)/ -l$(ACCRT_BASE2) $(OPENCL_LIB) $(CUDA_LIB)


ARC_CUDA: xsbench_ARC_CUDA
ARC_OPENCL: xsbench_ARC_OPENCL
all: all_CUDA
all_CUDA: xsbench_ARC_CUDA xsbench_serial
all_OPENCL: xsbench_ARC_OPENCL xsbench_serial
serial: xsbench_serial

C = Main.c CalculateXS.c GridInit.c Materials.c XSutils.c io.c
CPP = $(C:.c=.cpp)
CETUS_OUTPUT_CPP = $(addprefix $(CETUS_OUTPUT)/, $(CPP))

xsbench_serial: $(C) | target_dir
	$(CC) $(COMMONFLAGS) $(CC_FLAGS_SERIAL) -o $(BENCHMARK_SERIAL) $(C) $(CLIBS1)

xsbench_ARC_CUDA: $(CETUS_OUTPUT_CPP) | target_dir
	cd $(CETUS_OUTPUT); $(NVCC) $(CFLAGS2) -o ../$(BENCHMARK_CUDA) $(CPP) $(CLIBS2); cp openarc_kernel.cu ../$(TARGET); cp $(OPENARCLIB)/binBuilder_cuda ../$(TARGET); cp $(OPENARCLIB)/Timer ../$(TARGET); cd ../

xsbench_ARC_OPENCL: $(CETUS_OUTPUT_CPP) | target_dir
	cd $(CETUS_OUTPUT); $(OPENCL_CC) $(CFLAGS2) -o ../$(BENCHMARK_OPENCL) $(CPP) $(CLIBS2); cp openarc_kernel.cl ../$(TARGET); cp $(OPENARCLIB)/binBuilder_opencl ../$(TARGET); cp $(OPENARCLIB)/Timer ../$(TARGET);  cd ../

target_dir:
	mkdir -p $(TARGET)

clean:
	rm -f *.o;
	if [ -d ./$(CETUS_OUTPUT) ]; then cd $(CETUS_OUTPUT); rm -f *.o *~ memcheck memcheck.c; cd ..; fi
	rm -rf $(TARGET)
