XSBench
=======

Background
------------------------------------------------------

XSBench is a simple application that executes only the most
computationally expensive steps of Monte Carlo particle transport
\- the calculation of macroscopic cross sections \- in an effort to
expose bottlenecks within multi-core, shared memory architectures.

In a particle transport simulation, every time a particle changes
energy or crosses a material boundary, a new macroscopic cross
section must be calculated. The time spent looking up and calculating
the required cross section information often accounts for well over
half of the total active runtime of the simulation. XSBench uses a
unionized energy grid to facilitate cross section lookups for
randomized particle energies. There are a similar number of energy
grid points, material types, and nuclides as are used in the
Hoogenboom-Martin benchmark. The probability of particles residing
in any given material is weighted based on that material's commonality
inside the Hoogenboom-Martin benchmark geometry.

The end result is that the essential computational conditions and
tasks of fully featured Monte Carlo transport codes are retained
in XSBench, without the additional complexity and overhead inherent
in a fully featured code. This provides a much simpler and clearer
platform for stressing different architectures, and ultimately for
making determinations as to where hardware bottlenecks occur as
cores are added.


User Guide
------------------------------------------------------

By default, XSBench will with 1 thread per hardware core. If the
architecture supports hyperthreading, two threads will be run per
core.

To compile and run XSBench with default settings, use the following
commands:

>$ make

>$ make run

To compile with profiling (-pg) enabled:

>$ make profile

To alter the number of threads used, run the code with the desired
number of threads as an argument:

>$ ./XSBench 4

PAPI Performance Counters
------------------------------------------------------

By default, PAPI is disabled.

To enable PAPI, open the XSBench_header.h file and add (or uncomment)
the following definition to the file:

> #define __PAPI

Then, compile the code with the following command:

>$ make papi

Note that you may need to change the relevant library paths for papi
to work. The library path can be specified in the makefile, and the
header path is specified in the XSBench_header.h file.
