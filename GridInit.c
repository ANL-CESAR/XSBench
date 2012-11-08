#include "XSbench_header.h"

void generate_grids( NuclideGridPoint ** nuclide_grids,
                     int n_isotopes, int n_gridpoints ) {
	for( int i = 0; i < n_isotopes; i++ )
		for( int j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy =((double)rand() / (double)RAND_MAX);
			nuclide_grids[i][j].micro_xs =((double)rand() /(double)RAND_MAX );
		}
}

void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, int n_isotopes)
{
	if( DEBUG ) printf("Sorting Nuclide Energy Grids...\n");
	
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	for( int i = 0; i < n_isotopes; i++ )
		qsort( nuclide_grids[i], n_isotopes, sizeof(NuclideGridPoint),
		       cmp );
}

GridPoint * generate_energy_grid( int n_isotopes, int n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids) {
	if( DEBUG ) printf("Generating Unionized Energy Grid...\n");
	
	int n_unionized_grid_points = n_isotopes*n_gridpoints;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	GridPoint * energy_grid = (GridPoint *)malloc( n_unionized_grid_points
	                                               * sizeof( GridPoint ) );
	
	if( DEBUG ) printf("Copying and Sorting all nuclide grids...\n");
	
	NuclideGridPoint ** n_grid_sorted = gpmatrix( n_isotopes, n_gridpoints );
	
	memcpy( n_grid_sorted[0], nuclide_grids[0], n_isotopes*n_gridpoints*
	                                      sizeof( NuclideGridPoint ) );
	
	qsort( &n_grid_sorted[0][0], n_unionized_grid_points,
	       sizeof(NuclideGridPoint), cmp);

	if( DEBUG ) printf("Assigning energies to unionized grid...\n");
	
	for( int i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].energy = ( n_grid_sorted[0] + i)->energy;
	
	gpmatrix_free(n_grid_sorted);
	
	// also, need to allocate pointer arrays

	for( int i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].xs_ptrs = (NuclideGridPoint **)
		                         malloc( n_isotopes*sizeof(NuclideGridPoint*));

	return energy_grid;
}


void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
                    int n_isotopes, int n_gridpoints ){
	
	if( DEBUG ) printf("Assigning pointers to Unionized Energy Grid...\n");
	
	for( int i = 0; i < n_isotopes * n_gridpoints ; i++ )
	{
		// for each energy grid level, we need to look at all the nuclide
		// grids. For each nuclide grid, we need to scan through, and
		// find where our energy matches -OR- the first energy level that
		// is GREATER than (?) the quarry.
		double quarry = energy_grid[i].energy;
		if( DEBUG && i % 100 == 0 )
			printf("\rAligning Unionized Grid...(%.1lf%% complete)",
			       100.0 * (double) i / (n_isotopes*n_gridpoints) );
		for( int j = 0; j < n_isotopes; j++ )
		{
			int nuc_id = j;
			// Now scan through the relevant nuclide grid
			// This would be much better done as a binary search
			
			// log n binary search
			energy_grid[i].xs_ptrs[nuc_id] = 
				binary_search( nuclide_grids[nuc_id], quarry, n_gridpoints);

			// n search
			/*
			for( int k = 0; k < n_gridpoints; k++ )
			{
				if( k == n_gridpoints - 1 )
					energy_grid[i].xs_ptrs[nuc_id] = &nuclide_grids[nuc_id][k];
				else if( quarry <= nuclide_grids[nuc_id][k].energy )
					energy_grid[i].xs_ptrs[nuc_id] = &nuclide_grids[nuc_id][k];
			}
			*/
		}
	}
	if( DEBUG ) printf("\n");
}
