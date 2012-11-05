#include "XSbench_header.h"

double calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids,
                           int idx ){
	// pull ptr from energy grid
	NuclideGridPoint * high = energy_grid[idx].xs_ptrs[nuc];
	NuclideGridPoint * low = high - 1;
	
	// assignments
	double e_h, e_l, xs_h, xs_l, xs;
	e_h= high->energy;
	e_l = low->energy;
	xs_h = high->micro_xs;
	xs_l = low->micro_xs;

	// interpolate
	xs = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);

	return xs;
}

double calculate_macro_xs( double p_energy, int mat, int n_isotopes,
                           int n_gridpoints, int * num_nucs,
                           double ** concs, GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids,
													 int ** mats ){
	int p_nuc;
	int idx = 0;	
	double conc;
	double macro_xs = 0;

	// CONVERT TO BINARY SEARCH
	/*
	for( int i = 0; i < n_isotopes * n_gridpoints; i++ )
		if( energy_grid[i].energy <= p_energy )
			idx = i;
	*/

	// binary search
	idx = grid_search( n_isotopes * n_gridpoints, p_energy,
	                   energy_grid);	

	for( int j = 0; j < num_nucs[mat]; j++ )
	{
		p_nuc = mats[mat][j];
		conc = concs[mat][j];
		macro_xs += calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		                                n_gridpoints, energy_grid,
		                                nuclide_grids, idx ) * conc;
	}
	
	macro_xs = macro_xs / num_nucs[mat];
	
	return macro_xs;
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
