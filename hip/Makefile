#===============================================================================
# User Options
#===============================================================================

COMPILER    = amd
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

source = \
Main.hip \
io.hip \
Simulation.hip \
GridInit.hip \
XSutils.hip \
Materials.hip

obj = $(source:.hip=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS += -std=c++14

# Linker Flags
LDFLAGS = -lm

# AMD
ifeq ($(COMPILER),amd)
  CC = hipcc
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g -G
  LDFLAGS  += -g -G
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

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) XSbench_header.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.hip XSbench_header.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) XSbench_header.h

run:
	./$(program)
