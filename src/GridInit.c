#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

GridPoint * generate_hash_table( NuclideGridPoint ** nuclide_grids,
                          long n_isotopes, long n_gridpoints, long hash_bins )
{
	printf("Generating Hash Grid...\n");

	GridPoint * energy_grid = (GridPoint *)malloc( hash_bins * sizeof( GridPoint ) );
	int * full = (int *) malloc( n_isotopes * hash_bins * sizeof(int) );

	for( long i = 0; i < hash_bins; i++ )
		energy_grid[i].xs_ptrs = &full[n_isotopes * i];

	double du = 1.0 / hash_bins;

	// For each energy level in the hash table
	#pragma omp parallel for
	for( long e = 0; e < hash_bins; e++ )
	{
		double energy = e * du;

		// We need to determine the bounding energy levels for all isotopes
		for( long i = 0; i < n_isotopes; i++ )
		{
			energy_grid[e].xs_ptrs[i] = grid_search_nuclide( n_gridpoints, energy, nuclide_grids[i], 0, n_gridpoints-1);
		}
	}

	return energy_grid;
}

// Generates randomized energy grid for each nuclide
// Note that this is done as part of initialization (serial), so
// rand() is used.
void generate_grids( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints ) {
	for( long i = 0; i < n_isotopes; i++ )
		for( long j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy       =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].total_xs     =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].elastic_xs   =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].absorbtion_xs=((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].fission_xs   =((double)rand()/(double)RAND_MAX);
			nuclide_grids[i][j].nu_fission_xs=((double)rand()/(double)RAND_MAX);
		}
}

// Verification version of this function (tighter control over RNG)
void generate_grids_v( NuclideGridPoint ** nuclide_grids,
                     long n_isotopes, long n_gridpoints ) {
	for( long i = 0; i < n_isotopes; i++ )
		for( long j = 0; j < n_gridpoints; j++ )
		{
			nuclide_grids[i][j].energy       = rn_v();
			nuclide_grids[i][j].total_xs     = rn_v();
			nuclide_grids[i][j].elastic_xs   = rn_v();
			nuclide_grids[i][j].absorbtion_xs= rn_v();
			nuclide_grids[i][j].fission_xs   = rn_v();
			nuclide_grids[i][j].nu_fission_xs= rn_v();
		}
}

// Sorts the nuclide grids by energy (lowest -> highest)
void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, long n_isotopes,
                         long n_gridpoints )
{
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	for( long i = 0; i < n_isotopes; i++ )
		qsort( nuclide_grids[i], n_gridpoints, sizeof(NuclideGridPoint),
		       cmp );
	
	// error debug check
	/*
	for( int i = 0; i < n_isotopes; i++ )
	{
		printf("NUCLIDE %d ==============================\n", i);
		for( int j = 0; j < n_gridpoints; j++ )
			printf("E%d = %lf\n", j, nuclide_grids[i][j].energy);
	}
	*/
}

// Allocates unionized energy grid, and assigns union of energy levels
// from nuclide grids to it.
GridPoint * generate_energy_grid( long n_isotopes, long n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids) {
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	if( mype == 0 ) printf("Generating Unionized Energy Grid...\n");
	
	long n_unionized_grid_points = n_isotopes*n_gridpoints;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	GridPoint * energy_grid = (GridPoint *)malloc( n_unionized_grid_points
	                                               * sizeof( GridPoint ) );
	if( mype == 0 ) printf("Copying and Sorting all nuclide grids...\n");
	
	NuclideGridPoint ** n_grid_sorted = gpmatrix( n_isotopes, n_gridpoints );
	
	  	
	memcpy( n_grid_sorted[0], nuclide_grids[0], n_isotopes*n_gridpoints*
	                                      sizeof( NuclideGridPoint ) );
	
	qsort( &n_grid_sorted[0][0], n_unionized_grid_points,
	       sizeof(NuclideGridPoint), cmp);
	
	if( mype == 0 ) printf("Assigning energies to unionized grid...\n");
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].energy = n_grid_sorted[0][i].energy;
	

	gpmatrix_free(n_grid_sorted);
	
	int * full = (int *) malloc( n_isotopes * n_unionized_grid_points
	                             * sizeof(int) );
	if( full == NULL )
	{
		fprintf(stderr,"ERROR - Out Of Memory!\n");
		exit(1);
	}
	
	for( long i = 0; i < n_unionized_grid_points; i++ )
		energy_grid[i].xs_ptrs = &full[n_isotopes * i];
	
	// debug error checking
	/*
	for( int i = 0; i < n_unionized_grid_points; i++ )
		printf("E%d = %lf\n", i, energy_grid[i].energy);
	*/

	return energy_grid;
}

// Initializes the unionized energy grid, by locating the appropriate
// location in the nulicde grid for each entry on the unionized grid.
// This function should not be profiling when doing performance analysis or tuning.
// This is because this function does not really represent
// a real function in OpenMC, it is only necessary to initialize the UEG
// in this mini-app for later use in the simulation region of the code. 
// Analysts who find that this function is eating up a large portion of the runtime
// should exclude it from their analysis, or otherwise wash it out by increasing
// the amount of time spent in the simulation portion of the application by running
// more lookups with the "-l" command line argument.
// This function is not parallelized due to false sharing issues that arrise when writing
// to the UEG.
void initialization_do_not_profile_set_grid_ptrs( GridPoint * restrict energy_grid, NuclideGridPoint ** restrict nuclide_grids,
                    long n_isotopes, long n_gridpoints )
{
	int mype = 0;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if( mype == 0 ) printf("Assigning pointers to Unionized Energy Grid...\n");

	int * idx_low = (int *) calloc( n_isotopes, sizeof(int));
	double * energy_high = (double *) malloc( n_isotopes * sizeof(double));

	for( int i = 0; i < n_isotopes; i++ )
		energy_high[i] = nuclide_grids[i][1].energy;

	for( long e = 0; e < n_isotopes * n_gridpoints; e++ )
	{
		double unionized_energy = energy_grid[e].energy;
		for( long i = 0; i < n_isotopes; i++ )
		{
			if( unionized_energy < energy_high[i]  )
				energy_grid[e].xs_ptrs[i] = idx_low[i];
			else if( idx_low[i] == n_gridpoints - 2 )
				energy_grid[e].xs_ptrs[i] = idx_low[i];
			else
			{
				idx_low[i]++;
				energy_grid[e].xs_ptrs[i] = idx_low[i];
				energy_high[i] = nuclide_grids[i][idx_low[i]+1].energy;
			}

		}
	}

	free(idx_low);
	free(energy_high);

	/* Alternative method with high stride access, less efficient
	//#pragma omp parallel for schedule(dynamic,1)
	for( long i = 0; i < n_isotopes; i++ )
	{
		int nuclide_grid_idx_low = 0;
		double nuclide_energy_high = nuclide_grids[i][nuclide_grid_idx_low+1].energy;

		// Loop over entries in the UEG
		for( long e = 0; e < n_isotopes * n_gridpoints; e++ )
		{
			double unionized_energy = energy_grid[e].energy;

			if( unionized_energy < nuclide_energy_high  )
				energy_grid[e].xs_ptrs[i] = nuclide_grid_idx_low;
			else if( nuclide_grid_idx_low == n_gridpoints - 2 )
				energy_grid[e].xs_ptrs[i] = nuclide_grid_idx_low;
			else
			{
				nuclide_grid_idx_low++;
				energy_grid[e].xs_ptrs[i] = nuclide_grid_idx_low;
				nuclide_energy_high = nuclide_grids[i][nuclide_grid_idx_low+1].energy;
			}
		}
	}
	*/
}
