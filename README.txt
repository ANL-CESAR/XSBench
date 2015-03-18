==============================================================================
                   __   __ ___________                 _                        
                   \ \ / //  ___| ___ \               | |                       
                    \ V / \ `--.| |_/ / ___ _ __   ___| |__                     
                    /   \  `--. \ ___ \/ _ \ '_ \ / __| '_ \                    
                   / /^\ \/\__/ / |_/ /  __/ | | | (__| | | |                   
                   \/   \/\____/\____/ \___|_| |_|\___|_| |_|                   

                                   Version 13

                                  OpenACC Ports

==============================================================================
Contact Information
==============================================================================

Organizations:    

    Center for Exascale Simulation of Advanced Reactors (CESAR), 
        Argonne National Laboratory

    Future Technologies Group,
        Oak Ridge National Laboratory

Developers: 

    John Tramm   <jtramm@anl.gov>
    Ron Rahaman  <rahaman@anl.gov>
    Amanda Lund  <alund@anl.gov>
    Seyong Lee   <ornl.gov>

==============================================================================
What is XSBench?
==============================================================================

XSBench is a mini-app representing a key computational kernel of the Monte
Carlo neutronics application OpenMC.  A full explanation of the theory and
purpose of XSBench is provided in docs/XSBench_Theory.pdf.

This branch contains OpenACC ports of XSBench.  Currently, multiple OpenACC
ports are maintained for multiple compilers.  Please refer the README in each
subdirectory for more information on the respective port.  

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
