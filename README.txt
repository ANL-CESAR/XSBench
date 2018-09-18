==============================================================================
                   __   __ ___________                 _                        
                   \ \ / //  ___| ___ \               | |                       
                    \ V / \ `--.| |_/ / ___ _ __   ___| |__                     
                    /   \  `--. \ ___ \/ _ \ '_ \ / __| '_ \                    
                   / /^\ \/\__/ / |_/ /  __/ | | | (__| | | |                   
                   \/   \/\____/\____/ \___|_| |_|\___|_| |_|                   

                                   Version 18

==============================================================================
Contact Information
==============================================================================

Organization:     Center for Exascale Simulation of Advanced Reactors (CESAR)
                  Argonne National Laboratory

Development Lead: John Tramm   <jtramm@anl.gov>
                  Ron Rahaman  <rahaman@anl.gov>
                  Amanda Lund  <alund@anl.gov>

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
	from our github page.

	Git Repository Clone:
		
		Use the following command to clone XSBench to your machine:

		>$ git clone https://github.com/jtramm/XSBench.git

		Once cloned, you can update the code to the newest version
		using the following command (when in the XSBench directory):

		>$ git pull
	
	Git Zip Download:

		Simply use the "zip download" option on our webpage at:

		https://github.com/jtramm/XSBench

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
	  -m <simulation method>   Simulation method (history, event)
	  -t <threads>     Number of OpenMP threads to run
	  -s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)
	  -g <gridpoints>  Number of gridpoints per nuclide (overrides -s defaults)
	  -G <grid type>   Grid search type (unionized, nuclide, hash). Defaults to unionized.
	  -p <particles>   Number of particle histories
	  -l <lookups>     Number of Cross-section (XS) lookups per particle history
	  -h <hash bins>   Number of hash bins (only relevant when used with "-G hash")
	Default is equivalent to: -s large -l 34 -p 500000 -G unionized

	-m <simulation method>

		Sets the simulation method, either "history" or "event". These
		options represent the history based or event based algorithms
		respectively. The default is the history based method. These two
		methods represent different methods of parallelizing the Monte
		Carlo transport method. In the history based method, the central
		mode of parallelism is expressed over particles, which each require
		some number of macroscopic cross sections to be executed in series
		and in a dependent order. The event based method expresses its
		parallelism over a large pool of independent macroscopic cross
		section lookups that can be executed in any order without dependence.
		They key difference between the two methods is the dependence/independence
		of the macroscopic cross section loop.

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
	
	-G <grid type>

		Sets the grid search type (unionized, nuclide, hash). Defaults to unionized.
		The unionized grid is what is typically used in Monte Carlo codes, as
		it offers the fastest speed. However, the increase in speed comes in
		a significant increase in memory usage as a union of all the separate
		nuclide grids must be formed and stored in memory. The "nuclide" mode
		uses only the basic nuclide grid data, with no unionization. This is
		slower as a binary search must be performed on every nuclide for each
		macroscopic XS lookup, rather than only once when using the unionized
		grid.

		Finally, the "hash" mode is a newer algorithm now used by
		many full Monte Carlo codes which offers speed nearly equivalent
		to the unionized energy grid method, but with only a small fraction
		of the memory overhead. More details on the hash lookup algorithm
		can be found in the CHANGES file for version 16, or in the following
		publication:
  
		Forrest B Brown. New hash-based energy lookup algorithm for monte
		carlo codes. Trans. Am. Nucl. Soc., 111:659–662, 2014.
		http://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-14-27037

	-p <particles>
		
		Sets the number of particle histories to simulate.
		By default, this value is set to 500,000. Users may want to
		increase this value if they wish to extend the runtime of
		XSBench, perhaps to produce more reliable performance counter
		data - as extending the run will decrease the percentage of
		runtime spent on initialization. Real MC simulations in a full
		application may use up to several billion particles per generation,
		so there is great flexibility in this variable.

	-l <lookups>
		
		Sets the number of cross-section (XS) lookups to perform per particle.
		By default, this value is set to 34, which represents the average
		number of XS lookups per particle over the course of its lifetime in
		a light water reactor problem. Users should only alter this value if
		they are trying to capture the behavior of a different type of reactor
		(e.g., one with a fast spectrum), where the number of lookups per
		history may be different.

	-h <hash bins>

		Sets the number of hash bins (only relevant when using the hash
		lookup algorithm, as selected with "-G hash"). Default is 10,000.

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
BINARY_DUMP = no
BINARY_READ = no

-> Optimization enables the -O3 optimization flag.

-> Debugging enables the -g flag.

-> Profiling enables the -pg flag. When profiling the code, you may
   wish to significantly increase the number of lookups (with the -l
   flag) in order to wash out the initialization phase of the code.

-> MPI enables MPI support in the code.

-> The PAPI flag is explained below.

-> VEC_INFO enables some additional information regarding the success or
   failure of the compiler's use of vectorization techniques during
   compilation.

-> VERIFY enables a verification mode, the details of which are explained below.

-> Binary dump mode writes a binary file containing a randomized data set
   of cross sections. This can be used in tandem with the binary read mode
   to skip generation of cross section data every time the program is run.
   Note that if you create the grid when specifying the -G flag as
   "nuclide", data for the unionized energy grid will not be written, and
   therefore any subsequent runs using that file in binary read mode must
   also use the -G nuclide option. Files generated for the unionized grid
   can also be used when running in the nuclide grid mode.

-> Binary read mode reads the binary file created by the binary dump mode
   as a (usually) much faster substitution for randomly generating XS
   data on-the-fly. This mode is particularly useful if running on
   simulators where walltime minimization is extremely critical for
   logistical reasons.

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
the same hash number (i.e., increasing the number of particles, number
of gridpoints, etc, will result in different hashes). 

Note that the verification mode runs a little slower, due to need to hash
each macroscopic cross section result. Therefore, performance measurements
should generally not be made when verification mode is on. Rather,
verification mode should be used to ensure any changes to the code have not
broken it, and then be disabled before performance metrics are recorded.

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
Binary File Support
==============================================================================

The flags:

BINARY_DUMP = no
BINARY_READ = no

Can be set to yes in order to write or read a binary file containing
a randomized XS data set (both nuclide grids and unionized grids). This
feature may be extremely useful for users running on simulators where
walltime minimization is critical for logistical purposes, or for users
who are doing many sequential runs.

Note that identical input parameters (problem size, etc) must be used
when reading and writing a binary file. No runtime checks are made
to validate that the file correctly corresponds to the selected input
parameters.

Also note that if you create the grid when specifying the -G flag as
"nuclide", data for the unionized energy grid will not be written, and
therefore any subsequent runs using that file in binary read mode must
also use the "-G nuclide" option. Files generated for the full unionized grid
can also be used when running in the nuclide grid mode.

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
Citing XSBench
==============================================================================

Papers citing the XSBench program in general should refer to:

J. R. Tramm, A. R. Siegel, T. Islam, and M. Schulz, “XSBench - The
Development and Verification of a Performance Abstraction for Monte
Carlo Reactor Analysis,” presented at PHYSOR 2014 - The Role
of Reactor Physics toward a Sustainable Future, Kyoto.

A PDF of this paper can be accessed directly at this link:

http://www.mcs.anl.gov/papers/P5064-0114.pdf

Bibtex Entry:

@inproceedings{Tramm:wy,
author = {Tramm, John R and Siegel, Andrew R and Islam, Tanzima and Schulz, Martin},
title = {{XSBench - The Development and Verification of a Performance Abstraction for Monte Carlo Reactor Analysis}},
booktitle = {PHYSOR 2014 - The Role of Reactor Physics toward a Sustainable Future},
address = {Kyoto}
}

==============================================================================
