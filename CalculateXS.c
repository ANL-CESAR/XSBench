#include "XSbench_header.h"

void calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids,
                           int idx, double * xs_vector ){
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
	// Elastic XS
	xs_h = high->elastic_xs;
	xs_l = low->elastic_xs;
	xs_vector[1] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	// Absorbtion XS
	xs_h = high->absorbtion_xs;
	xs_l = low->absorbtion_xs;
	xs_vector[2] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	// Fission XS
	xs_h = high->fission_xs;
	xs_l = low->fission_xs;
	xs_vector[3] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	// Nu Fission XS
	xs_h = high->nu_fission_xs;
	xs_l = low->nu_fission_xs;
	xs_vector[4] = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
}

void calculate_macro_xs( double p_energy, int mat, int n_isotopes,
                           int n_gridpoints, int * num_nucs,
                           double ** concs, GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids, int ** mats,
                           double * macro_xs_vector ){
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
