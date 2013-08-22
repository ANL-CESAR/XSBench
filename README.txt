==============================================================================
                   __   __ ___________                 _                        
                   \ \ / //  ___| ___ \               | |                       
                    \ V / \ `--.| |_/ / ___ _ __   ___| |__                     
                    /   \  `--. \ ___ \/ _ \ '_ \ / __| '_ \                    
                   / /^\ \/\__/ / |_/ /  __/ | | | (__| | | |                   
                   \/   \/\____/\____/ \___|_| |_|\___|_| |_|                   

                                   Version 11

==============================================================================
Contact Information
==============================================================================

Organization:     Center for Exascale Simulation of Advanced Reactors (CESAR)
                  Argonne National Laboratory

Development Lead: John Tramm <jtramm@mcs.anl.gov>

==============================================================================
What is XSBench?
==============================================================================

XSBench is a mini-app representing a key computational kernel of the
Monte Carlo neutronics application OpenMC.

A full explanation of the theory and purpose of XSBench is provided in
docs/XSBench_Theory.pdf.

==============================================================================
Quick Start Guide
==============================================================================

Download----------------------------------------------------------------------

	For the most up-to-date version of XSBench, we recommend that you
	download from our git repository. This can be accomplished via
	cloning the repository from the command line, or by downloading a zip
	from our github page. Alternatively, you can download a tar file from
	the CESAR website directly.

	Git Repository Clone:
		
		Use the following command to clone XSBench to your machine:

		>$ git clone https://github.com/jtramm/XSBench.git

		Once cloned, you can update the code to the newest version
		using the following command (when in the XSBench directory):

		>$ git pull
	
	Git Zip Download:

		Simply use the "zip download" option on our webpage at:

		https://github.com/jtramm/XSBench
	
	CESAR Tar Download:

		A tar of the XSBench source code is available
		on the CESAR website at the following URL:
		
		https://cesar.mcs.anl.gov/content/software/neutronics

		Once downloaded, you can decompress XSBench using the following
		command on a linux or Mac OSX system:

		>$ tar -xvf XSBench-11.tar

		This will decompress the tar file into a directory called
		XSBench-11.

		To begin use of the XSBench code, you will have to navigate to
		the src directory:

		>$ cd XSBench-11/src

Compilation-------------------------------------------------------------------

	To compile XSBench with default settings, use the following
	command:

	>$ make

Running XSBench---------------------------------------------------------------

	To run XSBench with default settings, use the following command:

	>$ ./XSBench

	For non-default settings, XSBench supports the following command line
	options:	

	Usage: ./XSBench <options>
	Options include:
	  -t <threads>     Number of OpenMP threads to run
	  -s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)
	  -g <gridpoints>  Number of gridpoints per nuclide
	  -l <lookups>     Number of Cross-section (XS) lookups
	Default (no arguments given) is equivalent to: -s large -l 15000000

	-t <threads>

		Sets the number of OpenMP threads to run. By default, XSBench
		will run with 1 thread per hardware core. If the architecture
		supports hyperthreading, multiple threads will be run per
		core.

		If running in MPI mode, this will be the number of threads
		per MPI rank.

	-s <size>
		
		Sets the size of the Hoogenboom-Martin reactor model. There
		are four options: 'small', 'large', 'XL', and 'XXL'. By default,
		the 'large' option is selected. 

		The H-M size corresponds to the number of nuclides present
		in the fuel region.  The small version has 34 fuel nuclides,
		whereas the large version has 321 fuel nuclides. This
		significantly slows down the runtime of the program as the
		data structures are much larger, and more lookups are required
		whenever a lookup occurs in a fuel material.  Note that the
		program defaults to "Large" if no specification is made.

		The additional size options, "XL" and "XXL", do not directly correspond
		to any particular physical model. They are similar to the H-M
		"large" option, except the number of gridpoints per nuclide
		has been increased greatly. This creates an extremely
		large energy grid data structure (XL: 120GB, XXL: 252GB), which is
		unlikely to fit on a single node, but is useful for experimentation
		purposes on novel architectures.

	-g <gridpoints>

		Sets the number of gridpoints per nuclide. By default, this
		value is set to 11,303. This corresponds to the average number
		of actual gridpoints per nuclide in the H-M Large model as run
		by OpenMC with the actual ACE ENDF cross-section data. 

		Note that this option will override the number of default grid
		-points as set by the '-s' option.

	-l <lookups>
		
		Sets the number of cross-section (XS) lookups to perform. By
		default, this value is set to 15,000,000. Users may want to
		increase this value if they wish to extend the runtime of
		XSBench, perhaps to produce more reliable performance counter
		data - as extending the run will decrease the percentage of
		runtime spent on initialization.

==============================================================================
Debugging, Optimization & Profiling
==============================================================================

There are also a number of switches that can be set in the makefile.

Here is a sample of the control panel at the top of the makefile:

COMPILER  = gnu
OPTIMIZE  = yes
DEBUG     = no
PROFILE   = no
MPI       = no
PAPI      = no
VEC_INFO  = no
VERIFY    = no
PAUSE     = no
BENCHMARK = no

-> Optimization enables the -O3 optimization flag.

-> Debugging enables the -g flag.

-> Profiling enables the -pg flag.

-> MPI enables MPI support in the code.

-> The PAPI flag is explained below.

-> VEC_INFO enables some additional information regarding the success or
failure of the compiler's use of vectorization techniques during
compilation.

-> VERIFY enables a verification mode, the details of which are explained below.

-> Pause adds an arbitrary pause in between macro-XS lookups, in order to
represent particle tracking / communication etc that occurs in the real
OpenMC. Note that the amount of time spent waiting

-> Benchmark mode tests all possible thread configurations on the given
computer. I.e., if your computer supports 16 threads, XSBench will
automatically do 1 <= nthreads <= 16 lookup loops

==============================================================================
MPI Support
==============================================================================

While XSBench is primarily used to investigate "on node parallelism" issues,
some systems provide power & performance statistics batched in multi-node
configurations. To accommodate this, XSBench provides an MPI mode which
runs the code on all MPI ranks simultaneously. There is no decomposition
across ranks of any kind, and all ranks accomplish the same work. There is
only one point of MPI communication (a reduce) at the end, which aggregates
the timing statistics and averages them across MPI ranks before printing them
out.

MPI support can be enabled with the makefile flag "MPI". If you are not using
the mpicc wrapper on your system, you may need to alter the makefile to
make use of your desired compiler.


==============================================================================
Verification Support
==============================================================================

XSBench has the ability to verify that consistent and correct results are
achieved. This mode is enabled by altering the "VERIFY" setting to 'yes' in
the makefile, i.e.:

VERIFY = yes

Once enabled, the code will generate a hash of the results and display it
with the other data once the code has completed executing. This hash can
then be verified against hashes that other versions or configurations of
the code generate. For instance, running XSBench with 4 threads vs 8 threads
(on a machine that supports that configuration) should generate the
same hash number. Changing the model / run parameters should NOT generate
the same hash number (i.e., increasing the number of lookups, number
of gridpoints, etc, will result in different hashes). 

Verification mode uses a RNG with a static seed. The randomized lookup
parameters are generated within a critical region. This ensures that the
same set of lookups are performed regardless of the number of threads
used. Then, after each lookup is completed, the lookup parameters and
the cross section vector are hashed together. This local hash is then
atomically added to a global running hash.

Note that the verification mode runs much slower, due to the use of
atomics within the threading loop. 

Below are the expected checksums for default runs of each size (-s):

small : 74966788162
large : 74994938929 

==============================================================================
PAPI Performance Counters
==============================================================================

PAPI performance counters is a performance counting library that can
offer information regarding the frequency of specific events (such as
memory loads, cache misses, branch prediction failures, etc) that occur
when the code is executed. XSBench supports use of these performance
counters, although it is left to the user to select the particular
performance counters and locations to instrument.

By default, PAPI is disabled.

To enable PAPI, set in the makefile:

PAPI = yes

Note that you may need to change the relevant library paths for papi to
work (as these are dependent on your machine).  The library path can be
specified in the makefile, and the header path is specified in the
XSBench_header.h file.

To select the performance counters you are interested in, open
the file papi.c and alter the events[] array to the events
you would like to count. 

==============================================================================
Running on ANL BlueGene/Q (Vesta & Mira)
==============================================================================

Compilation is done using the included makefile, as follows:

>$ make MACHINE=bluegene

Note that the INFO macro in the XSbench_header.h file should be set to
0 when running on BG/Q to remove the run status portions of the output,
which cuts down on unnecessary file I/O, i.e.:

#define INFO 0

Also, note that you may need to add the following line to your .soft
file in order to use the mpicc compiler wrapper:

+mpiwrapper-gcc

Then, be sure to use the "resoft" command to update your software, i.e.,:

>$ resoft

When running in c16 mode, the maximum number of gridpoints per nuclide
is 900 (when running in "Large" mode). More points will cause the 1GB
memory limit to be broken.

A basic test run on 1 node can be achieved (assuming you have an allocation)
using the makefile and the following command:

>$ make bgqrun

Further information on queuing can be found at:

https://www.alcf.anl.gov/resource-guides/vesta-queuing

==============================================================================
Adding Extra Flops
==============================================================================

One of the areas we're investigating is the effect of adding additional
flops per each load from the nuclide xs arrays. Adding flops has so far
shown to increase scaling, indicating that there is in fact a bottleneck
being caused by the memory loads.

To enable this feature, go to the XSBench_header.h file and uncomment
out the "#define ADD_EXTRAS" line.

By default, there are no "extra" flops in the benchmark. To add these
in, open the XSbench_header.h file and change the EXTRA_FLOPS definition
to the desired number of additional flops. i.e.,

"#define EXTRA_FLOPS 10"

Note that even just 5-10 extra flops will greatly increase the running
time of the program.

==============================================================================
Adding Extra Memory Loads
==============================================================================

One of the areas we're investigating is the effect of adding additional
memory loads per each required load from the nuclide xs arrays. Adding
loads has so far been rather inconclusive, as randomizing the loads
requires a number of flops be performed to generate the random index
number to load from.

To enable this feature, go to the XSBench_header.h file and uncomment
out the "#define ADD_EXTRAS" line.

By default, there are no "extra" loads in the benchmark. To add these
in, open the XSbench_header.h file and change the EXTRA_LOADS definition
to the desired number of additional loads. i.e.,

"#define EXTRA_LOADS 10"

Note that even just a few extra loads will greatly increase the running
time of the program.

==============================================================================
