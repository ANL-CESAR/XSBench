#include "XSbench_header.h"

void generate_grids( NuclideGridPoint ** nuclide_grids,
                     int n_isotopes, int n_gridpoints ) {
	for( int i = 0; i < n_isotopes; i++ )
		for( int j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy =((double)rand() / (double)RAND_MAX);
			nuclide_grids[i][j].total_xs =((double)rand() /(double)RAND_MAX );
			nuclide_grids[i][j].elastic_xs =((double)rand() /(double)RAND_MAX );
			nuclide_grids[i][j].absorbtion_xs=((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].fission_xs =((double)rand() /(double)RAND_MAX );
			nuclide_grids[i][j].nu_fission_xs=((double)rand()/(double)RAND_MAX);
		}
}

void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, int n_isotopes)
{
	printf("Sorting Nuclide Energy Grids...\n");
	
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	for( int i = 0; i < n_isotopes; i++ )
		qsort( nuclide_grids[i], n_isotopes, sizeof(NuclideGridPoint),
		       cmp );
}

GridPoint * generate_energy_grid( int n_isotopes, int n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids) {
	printf("Generating Unionized Energy Grid...\n");
	
	int n_unionized_grid_points = n_isotopes*n_gridpoints;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	GridPoint * energy_grid = (GridPoint *)malloc( n_unionized_grid_points
	                                               * sizeof( GridPoint ) );
	
	printf("Copying and Sorting all nuclide grids...\n");
	
	NuclideGridPoint ** n_grid_sorted = gpmatrix( n_isotopes, n_gridpoints );
	
	memcpy( n_grid_sorted[0], nuclide_grids[0], n_isotopes*n_gridpoints*
	                                      sizeof( NuclideGridPoint ) );
	
	qsort( &n_grid_sorted[0][0], n_unionized_grid_points,
	       sizeof(NuclideGridPoint), cmp);

	printf("Assigning energies to unionized grid...\n");
	
	for( int i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].energy = ( n_grid_sorted[0] + i)->energy;
	
	gpmatrix_free(n_grid_sorted);
	
	for( int i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].xs_ptrs = (NuclideGridPoint **)
		                         malloc( n_isotopes*sizeof(NuclideGridPoint*));

	return energy_grid;
}


void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
                    int n_isotopes, int n_gridpoints ){
	
	printf("Assigning pointers to Unionized Energy Grid...\n");
	
	for( int i = 0; i < n_isotopes * n_gridpoints ; i++ )
	{
		double quarry = energy_grid[i].energy;
		if( INFO && i % 500 == 0 )
			printf("\rAligning Unionized Grid...(%.0lf%% complete)",
			       100.0 * (double) i / (n_isotopes*n_gridpoints) );
		for( int j = 0; j < n_isotopes; j++ )
		{
			int nuc_id = j;
			// log n binary search
			energy_grid[i].xs_ptrs[nuc_id] = 
				binary_search( nuclide_grids[nuc_id], quarry, n_gridpoints);
		}
	}
	printf("\n");
}
