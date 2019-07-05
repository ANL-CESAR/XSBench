==============================================================================
                   __   __ ___________                 _                        
                   \ \ / //  ___| ___ \               | |                       
                    \ V / \ `--.| |_/ / ___ _ __   ___| |__                     
                    /   \  `--. \ ___ \/ _ \ '_ \ / __| '_ \                    
                   / /^\ \/\__/ / |_/ /  __/ | | | (__| | | |                   
                   \/   \/\____/\____/ \___|_| |_|\___|_| |_|                   

                                   Version 19

==============================================================================
Contact Information
==============================================================================

Organization:     Computational Science Division
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
docs/XSBench_Theory.pdf. More information is also available in the publications
listed at the bottom of this document.

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

Selecting A Source Version----------------------------------------------------

	XSBench has been implemented in multiple different languages to target
	a variety of computational architectures and accelerators. The available
	implementations can be found in the "XSBench/src" directory:

	1. XSBench/src/openmp-threading

	This is the "default" version of XSBench that is appropriate for serial
	and multicore CPU architectures. The method of parallelism is via the
	OpenMP threading model.

	2. XSBench/src/openmp-offload

	This method of parallelism uses OpenMP 4.5 (or newer) to map program
	data to a remote accelerator memory space and run targeted kernels on
	the accelerator. This method of parallelism could be used for a wide
	variety of architectures (besides CPUs) that support OpenMP 4.5 targeting.

	NOTE: The Makefile will likely not work by default and will need to be
	adjusted to utilize your OpenMP accelerator compiler.

	3. XSBench/src/CUDA

	This version of XSBench is written in CUDA for use with NVIDIA GPU
	architectures.

	NOTE: You will likely want to specify in the makefile the SM version
	for the card you are running on.

	4. XSBench/src/OpenCL

	This version of XSBench is written in OpenCL, and can be used for CPU,
	GPU, FPGA, or other architecture that supports OpenCL. It was written
	with GPUs in mind, so if running on other architectures you may need to
	heavily re-optimize the code. You will also likely need to edit the
	makefile to supply the path to your OpenCL compiler.
	

Compilation-------------------------------------------------------------------

	To compile XSBench with default settings, navigate to your selected src
	directory and use the following command:

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
	  -b <binary mode> Read or write all data structures to file. If reading, this will skip initialization phase. (read, write)
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

	-b <binary mode>

		This optional mode can read or write the simulation data structures
		to disk. Options are ("read" or "write"). This may be useful if
		it is necessary to minimize the initialization phase of the program,
		which has a non-trivial runtime. The generated file is named
		"XS_data.dat" and will be located in the current working directory.
		The same file name and location will be used when reading. Note that
		as the file is binary, it may not be portable between compilers and
		computer systems. NOTE: When running in the "read" mode, you must
		be running with an identical program configuration as when the file
		was generated. E.g., if the file was generated with the "-G nuclide"
		argument, subsequent runs reading from that file must use the same
		configuration flags.

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

-> Optimization enables the -O3 optimization flag.

-> Debugging enables the -g flag.

-> Profiling enables the -pg flag. When profiling the code, you may
   wish to significantly increase the number of lookups (with the -l
   flag) in order to wash out the initialization phase of the code.

-> MPI enables MPI support in the code.

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

Legacy versions of XSBench had a special "Veriication" compiler flag option
to enable verification of the results. However, a much more performant and
portable verification scheme was developed and is now used for all
configurations -- therefore, it is not necessary to compile with or without
the verification mode as it is always enabled by default.

XSBench generates a hash of the results at the end of the simulation and displays
it with the other data once the code has completed executing. This hash can
then be verified against hashes that other versions or configurations of
the code generate. For instance, running XSBench with 4 threads vs 8 threads
(on a machine that supports that configuration) should generate the
same hash number. Changing the model / run parameters should NOT generate
the same hash number (i.e., increasing the number of particles, number
of gridpoints, etc, will result in different hashes). However, changing
the type of lookup performed (e.g., nuclide, unionized, or hash) should result
in the same hash being generated. Changing the simulation mode (history or
event) will generate different hashes. 

==============================================================================
Binary File Support
==============================================================================

Instead of initializing the randomized synthetic cross section data structres
in XSBench everytime it is run, you may optionally have XSBench generate
a data set and write it to file. It can then be read on subsequent runs to
speed up initialization. This process is controlled with the 
"-b (read, write)" command line argument.

Can be set to yes in order to write or read a binary file containing
a randomized XS data set (both nuclide grids, hash grids, and unionized grids).
This feature may be extremely useful for users running on simulators where
walltime minimization is critical for logistical purposes, or for users
who are doing many sequential runs.

Note that identical input parameters (problem size, etc) must be used
when reading and writing a binary file. No runtime checks are made
to validate that the file correctly corresponds to the selected input
parameters.

Also note that if you create the grid when specifying the -G flag as
"nuclide", data for the unionized energy grid will not be written, and
therefore any subsequent runs using that file in binary read mode must
also use the "-G nuclide" option.

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
