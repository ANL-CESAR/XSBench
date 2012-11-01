#include "XSbench_header.h"

double calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids ){
	// find ptr in energy grid
	NuclideGridPoint * high = energy_grid[1].xs_ptrs[1];
	NuclideGridPoint * low;
	for( int i = 0; i < n_isotopes * n_gridpoints; i++ )
		if( energy_grid[i].energy <= p_energy )
			high = energy_grid[i].xs_ptrs[nuc];
	
	low = high - 1;
	
	double e_h, e_l, xs_h, xs_l, xs;
	e_h= high->energy;
	e_l = low->energy;
	xs_h = high->micro_xs;
	xs_l = low->micro_xs;

	// interpolate
	xs = xs_h - (e_h - p_energy) * (xs_h - xs_l) / (e_h - e_l);
	return xs;
}
