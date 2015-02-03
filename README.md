# OCCA-XSBench

## Introduction

OCCA-XSBench is a port of XSBench, a neutronics proxy-app. It uses the OCCA
multithreading API to target both multicore CPUs and many-core hardware
accellerators.  OCCA-XSBench features a unified, highly-optimized kernel that
may be compiled for either OpenMP and CUDA.  

XSBench and OCCA were developed by research groups affiliated with the [Center
for Exascale Simulation of Advanced Reactors](https://cesar.mcs.anl.gov/), a
DOE Office of Science Co-design Center.  

[XSBench][xsbench] is a mini-app representing a key computational kernel of the
Monte Carlo neutronics application [OpenMC][openmc].  It was originally
developed by John Tramm, Andrew Siegel, et al. in the Mathematics and Computer
Science Division at Argonne National Laboratory.  A full explanation of the
theory and purpose of XSBench is provided in `docs/XSBench_Theory.pdf`.  

[OCCA][occa] is an extensible multi-threading programming API that allows a
unified source code to target many multi-threading languages, including OpenMP,
CUDA, and OpenCL.  It was originally developed by David Medina, Tim Warburton,
et. al in the Department of Computational and Applied Mathematics at Rice
University.  Documentation and downloads for OCCA can be found at
[libocca.org][occa].

For more information about OCCA-XSBench, please contact Ron
Rahaman \(rahaman@anl.gov \) and Amanda Lund \(alund@anl.gov \) at Argonne National
Laboratory.  

[xsbench]: https://github.com/ANL-CESAR/XSBench "XSBench"
[openmc]:  https://mit-crpg.github.io/openmc/ "OpenMC"
[occa]:  http://libocca.org "OCCA"


## Installation

#### Dependencies

The OCCA library must first be installed and configured.  OCCA-XSBench
currently supports OpenMP and CUDA modes; hence, OCCA must be compiled with
OpenMP and CUDA modes enabled. 

#### Compilation

The XSBench driver may be compiled with `make`.  The makefile
(`src/Makefile`) has the following settings:

* `COMPILER = ( gnu | intel ):` Specifies compiler for driver.
* `OPTIMIZE = ( yes | no ):` Sets the -O3 flag.  
* `DEBUG  = ( yes | no ):` Sets the -g flag.
* `PROFILE = ( yes | no );` Sets the -pg flag.
* `VERIFY = ( yes | no ):` Enables verification mode.  See [Verification](#verification), below.

The XSBench kernel is compiled at runtime using the environment variabls:
`$OCCA_OPENMP_COMPILER`, `$OCCA_OPENMP_COMPILER_FLAGS`, `$OCCA_CUDA_COMPILER`,
and `$OCCA_CUDA_COMPILER_FLAGS`.  

Currently, the multithreading mode \(OpenMP or CUDA\) is specified by a macro
in `src/Main.c`.  Future versions will allow the mode to be set at runtime.  

##Usage

#### Synopsis

`./XSBench [-t <threads>] [-s 'small'|'large'|'XL'|'XXL'] [-g <gridpoints>] [-l <lookups>]`

#### Options

- `-t <threads>`

  Sets the number of OpenMP threads to run. Defaults to the value of
  `$OMP_NUM_THREADS`, if it is set. Otherwise, defaults to
  `omp_get_num_procs()`.
 

- `-s 'small'|'large'|'XL'|'XXL'`

  Sets the size of the Hoogenboom-Martin reactor model. Defaults to `'large'`.
  
  The H-M size corresponds to the number of nuclides present
  in the fuel region.  The `'small'` version has 34 fuel nuclides,
  whereas the `'large'` version has 321 fuel nuclides. This
  significantly slows down the runtime of the program as the
  data structures are much larger, and more lookups are required
  whenever a lookup occurs in a fuel material.    

  The additional size options, `'XL'` and `'XXL'`, do not directly correspond
  to any particular physical model. They are similar to the H-M "large" option,
  except the number of gridpoints per nuclide has been increased greatly. This
  creates an extremely large energy grid data structure \(XL: 120GB, XXL:
  252GB\), which is unlikely to fit on a single node, but is useful for
  experimentation purposes on novel architectures.

- `-g <gridpoints>`

  Sets the number of gridpoints per nuclide. Defaults to 11,303. 
  
  This corresponds to the average number of actual gridpoints per nuclide in
  the H-M Large model as run by OpenMC with the actual ACE ENDF cross-section
  data. 

- `-l <lookups>`

  Sets the number of cross-section \(XS\) lookups to perform. Defaults to
  15,000,000. 
  
  Users may want to increase this value if they wish to extend the runtime of
  XSBench, perhaps to produce more reliable performance counter data - as
  extending the run will decrease the percentage of runtime spent on
  initialization.

#### Verification

  If compiled in verification mode \(see [Compilation](#compilation)\), a weak test for
  consistency is performed.  The program will accumulate cross-section values
  and report an non-deteriministic expectation value \(approximately 23.5 +/-
  0.5\ for the 'small' benchmark and 73.2 +/- 0.2 for the 'large' benchmark).  

  The verification does not produce a deterministic value because the random
  number stream is not reproducible and because the floating-point arithmatic
  is non-associative.  Future versions will produce a deterministic integer
  checksum.  

## References

- J. R. Tramm, A. R. Siegel, T. Islam, and M. Schulz. "XSBench - The
  Development and Verification of a Performance Abstraction for Monte
  Carlo Reactor Analysis." *PHYSOR 2014 - The Role
  of Reactor Physics toward a Sustainable Future*, Kyoto, 2014. [pdf](http://www.mcs.anl.gov/papers/P5064-0114.pdf)
- D. S. Medina, M. St-Cyr, T. Warburton.  "OCCA: A unifiied approach to
  multi-threading languages." *CoRR* abs/1403.0968, 2014. [pdf](http://arxiv.org/abs/1403.0968) 
