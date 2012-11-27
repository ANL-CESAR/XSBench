#include "XSbench_header.h"
/*
To Do: I want to add a variable number of memory loads here.
How do I do this without adding way more flops than loads?
IDEAS:
1. Just do N extra loads from random locations. (Requires use of RNG).
   Obviously, the RNG results in a non-trivial amount of extra flops.
2. Load from some modular location ahead in memory, repeat pattern N
   times. This avoids the RNG, but still requires a few flops per load.
3. Just load the next N elements off the grid. There will be huge
   cacheing efficiency here, but just the same if N is high enough
   it should put a lot of stress on memory badwidth, which is the whole
   idea of adding this lever. Problem here is that we have no idea
   where the original element was on the grid.....

I suppose as of right now, I kind of like 3 best, since it doesn't
add many flops, but does add a ton of loads (even if the effect will
be proportionally less ( N % cache ). 
*/

void calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids,
                           int idx, double * restrict xs_vector ){
	// pull ptr from energy grid
	NuclideGridPoint * high = energy_grid[idx].xs_ptrs[nuc];
	NuclideGridPoint * low = high - 1;
	
	// assignments
	double e_h, e_l, xs_h, xs_l;
	e_h= high->energy;
	e_l = low->energy;

	// Total XS
	xs_h = high->total_xs;
	xs_l = low->total_xs;
	xs_vector[0] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	do_flops();
	do_loads( nuclide_grids, high, n_gridpoints );
	// Elastic XS
	xs_h = high->elastic_xs;
	xs_l = low->elastic_xs;
	xs_vector[1] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	do_flops();
	// Absorbtion XS
	xs_h = high->absorbtion_xs;
	xs_l = low->absorbtion_xs;
	xs_vector[2] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	do_flops();
	// Fission XS
	xs_h = high->fission_xs;
	xs_l = low->fission_xs;
	xs_vector[3] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	do_flops();
	// Nu Fission XS
	xs_h = high->nu_fission_xs;
	xs_l = low->nu_fission_xs;
	xs_vector[4] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	do_flops();
}

/*
As we make the jump from 1 -> cores, the amount of time spent in macro_xs
jumps significantly (from 36 to 170 cpu seconds). Nearly a factor of 5x. This
jump occurs while practically everything else stays the same (as expected).

So, the question is: What's going on in macro_xs that causes poor scaling?
It's not grid search (that's measured to stay the same). What else could
it be?

My initial thoughts was that the xs_vector allocation was expensive, but
now I'm thinking that's not a big deal. If that were the problem, we'd be
going slow, but scaling would be just fine.

*/
void calculate_macro_xs( double p_energy, int mat, int n_isotopes,
                           int n_gridpoints, int * restrict num_nucs,
                           double ** restrict concs,
						   GridPoint * restrict energy_grid,
                           NuclideGridPoint ** restrict nuclide_grids,
						   int ** restrict mats,
                           double * restrict macro_xs_vector ){
	//double * xs_vector = (double *) malloc( 5 * sizeof(double) );
	double xs_vector[5];
	int p_nuc;
	int idx = 0;	
	double conc;

	// cleans out macro_xs_vector
	for( int k = 0; k < 5; k++ )
		macro_xs_vector[k] = 0;

	// CONVERT TO BINARY SEARCH
	/*
	for( int i = 0; i < n_isotopes * n_gridpoints; i++ )
		if( energy_grid[i].energy <= p_energy )
			idx = i;
	*/

	// binary search
	idx = grid_search( n_isotopes * n_gridpoints, p_energy,
	                   energy_grid);	

	// I think our contention problem may be here.
	// Think about it. This is a pretty small set of information compared
	// to the giant energy grid arrays.
	for( int j = 0; j < num_nucs[mat]; j++ )
	{
		p_nuc = mats[mat][j];
		conc = concs[mat][j];
		calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		                    n_gridpoints, energy_grid,
		                    nuclide_grids, idx, xs_vector );
		for( int k = 0; k < 5; k++ )
			macro_xs_vector[k] += xs_vector[k] * conc;
	}
	//free(xs_vector);
	
	for( int k = 0; k < 5; k++ )
		macro_xs_vector[k] = macro_xs_vector[k] / num_nucs[mat];
}

// binary search for energy on unionized energy grid
int grid_search( int n, double quarry, GridPoint * A)
{
	int min = 0;
	int max = n-1;
	int mid = 0;
	
	while( max >= min )
	{
		mid = min + floor( (max-min) / 2.0);
		if( A[mid].energy < quarry )
			min = mid+1;
		else if( A[mid].energy > quarry )
			max = mid-1;
		else
			return mid;
	}
	return mid;
}
