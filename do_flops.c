#include "XSbench_header.h"

// Note that this file needs to be compiled separately with the -O0 flag
// to prevent the compiler from optimizing away the flops. The included
// makefile takes this into account.
void do_flops( void )
{
	double a = 1.33;
	double b = 2.34;
	int i;
	
	for( i = 0; i < EXTRA_FLOPS; i++ )
	{
		a = a * b;
	}
}

// Due to the way I did my logic, there's no way of knowing where
// "high" is exactly in its nuclide_grids array. It could be at
// the end, or at the beginning, we have no idea. 
void do_loads( NuclideGridPoint ** restrict nuclide_grids,
          NuclideGridPoint * high, int n_gridpoints )
{
	int i;
	double load;
	for( i = 0; i < EXTRA_LOADS; i++ )
	{
		;	
	}
}
