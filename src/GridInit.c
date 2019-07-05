#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int double_compare(const void * a, const void * b)
{
	double A = *((double *) a);
	double B = *((double *) b);

	if( A > B )
		return 1;
	else if( A < B )
		return -1;
	else
		return 0;
}

int NGP_compare(const void * a, const void * b)
{
	NuclideGridPoint A = *((NuclideGridPoint *) a);
	NuclideGridPoint B = *((NuclideGridPoint *) b);

	if( A.energy > B.energy )
		return 1;
	else if( A.energy < B.energy )
		return -1;
	else
		return 0;
}

SimulationData flat_grid_init( Inputs in )
{
	// Structure to hold all allocated simuluation data arrays
	SimulationData SD;

	// Keep track of how much data we're allocating
	size_t nbytes = 0;

	////////////////////////////////////////////////////////////////////
	// Initialize Nuclide Grids
	////////////////////////////////////////////////////////////////////

	// First, we need to initialize our nuclide grid. This comes in the form
	// of a flattened 2D array that hold all the information we need to define
	// the cross sections for all isotopes in the simulation. 
	// The grid is composed of "NuclideGridPoint" structures, which hold the
	// energy level of the grid point and all associated XS data at that level.
	// An array of structures (AOS) is used instead of
	// a structure of arrays, as the grid points themselves are accessed in 
	// a random order, but all cross section interaction channels and the
	// energy level are read whenever the gridpoint is accessed, meaning the
	// AOS is more cache efficient.
	
	// Initialize Nuclide Grid
	SD.length_nuclide_grid = in.n_isotopes * in.n_gridpoints;
	SD.nuclide_grid     = (NuclideGridPoint *) malloc( SD.length_nuclide_grid * sizeof(NuclideGridPoint));
	assert(SD.nuclide_grid != NULL);
	nbytes += SD.length_nuclide_grid * sizeof(NuclideGridPoint);
	for( int i = 0; i < SD.length_nuclide_grid; i++ )
	{
		SD.nuclide_grid[i].energy        = rn_v();
		SD.nuclide_grid[i].total_xs      = rn_v();
		SD.nuclide_grid[i].elastic_xs    = rn_v();
		SD.nuclide_grid[i].absorbtion_xs = rn_v();
		SD.nuclide_grid[i].fission_xs    = rn_v();
		SD.nuclide_grid[i].nu_fission_xs = rn_v();
	}

	// Sort so that each nuclide has data stored in ascending energy order.
	for( int i = 0; i < in.n_isotopes; i++ )
		qsort( &SD.nuclide_grid[i*in.n_gridpoints], in.n_gridpoints, sizeof(NuclideGridPoint), NGP_compare);
	
	// error debug check
	/*
	for( int i = 0; i < in.n_isotopes; i++ )
	{
		printf("NUCLIDE %d ==============================\n", i);
		for( int j = 0; j < in.n_gridpoints; j++ )
			printf("E%d = %lf\n", j, SD.nuclide_grid[i * in.n_gridpoints + j].energy);
	}
	*/
	

	////////////////////////////////////////////////////////////////////
	// Initialize Acceleration Structure
	////////////////////////////////////////////////////////////////////
	
	if( in.grid_type == NUCLIDE )
	{
		SD.length_unionized_energy_array = 0;
		SD.length_index_grid = 0;
	}
	
	if( in.grid_type == UNIONIZED )
	{
		// Allocate space to hold the union of all nuclide energy data
		SD.length_unionized_energy_array = in.n_isotopes * in.n_gridpoints;
		SD.unionized_energy_array = (double *) malloc( SD.length_unionized_energy_array * sizeof(double));
		assert(SD.unionized_energy_array != NULL );
		nbytes += SD.length_unionized_energy_array * sizeof(double);

		// Copy energy data over from the nuclide energy grid
		for( int i = 0; i < SD.length_unionized_energy_array; i++ )
			SD.unionized_energy_array[i] = SD.nuclide_grid[i].energy;

		// Sort unionized energy array
		qsort( SD.unionized_energy_array, SD.length_unionized_energy_array, sizeof(double), double_compare);

		// Allocate space to hold the acceleration grid indices
		SD.length_index_grid = SD.length_unionized_energy_array * in.n_isotopes;
		SD.index_grid = (int *) malloc( SD.length_index_grid * sizeof(int));
		assert(SD.index_grid != NULL);
		nbytes += SD.length_index_grid * sizeof(int);

		// Generates the double indexing grid
		int * idx_low = (int *) calloc( in.n_isotopes, sizeof(int));
		assert(idx_low != NULL );
		double * energy_high = (double *) malloc( in.n_isotopes * sizeof(double));
		assert(energy_high != NULL );

		for( int i = 0; i < in.n_isotopes; i++ )
			energy_high[i] = SD.nuclide_grid[i * in.n_gridpoints + 1].energy;

		for( long e = 0; e < SD.length_unionized_energy_array; e++ )
		{
			double unionized_energy = SD.unionized_energy_array[e];
			for( long i = 0; i < in.n_isotopes; i++ )
			{
				if( unionized_energy < energy_high[i]  )
					SD.index_grid[e * in.n_isotopes + i] = idx_low[i];
				else if( idx_low[i] == in.n_gridpoints - 2 )
					SD.index_grid[e * in.n_isotopes + i] = idx_low[i];
				else
				{
					idx_low[i]++;
					SD.index_grid[e * in.n_isotopes + i] = idx_low[i];
					energy_high[i] = SD.nuclide_grid[i * in.n_gridpoints + idx_low[i] + 1].energy;	
				}
			}
		}

		free(idx_low);
		free(energy_high);
	}

	if( in.grid_type == HASH )
	{
		SD.length_unionized_energy_array = 0;
		SD.length_index_grid  = in.hash_bins * in.n_isotopes;
		SD.index_grid = (int *) malloc( SD.length_index_grid * sizeof(int)); 
		assert(SD.index_grid != NULL);
		nbytes += SD.length_index_grid * sizeof(int);

		double du = 1.0 / in.hash_bins;

		// For each energy level in the hash table
		#pragma omp parallel for
		for( long e = 0; e < in.hash_bins; e++ )
		{
			double energy = e * du;

			// We need to determine the bounding energy levels for all isotopes
			for( long i = 0; i < in.n_isotopes; i++ )
			{
				SD.index_grid[e * in.n_isotopes + i] = grid_search_nuclide( in.n_gridpoints, energy, SD.nuclide_grid + i * in.n_gridpoints, 0, in.n_gridpoints-1);
			}
		}
	}

	////////////////////////////////////////////////////////////////////
	// Initialize Materials and Concentrations
	////////////////////////////////////////////////////////////////////
	
	// Set the number of nuclides in each material
	SD.num_nucs  = load_num_nucs(in.n_isotopes);
	SD.length_num_nucs = 12; // There are always 12 materials in XSBench

	// Intialize the flattened 2D grid of material data. The grid holds
	// a list of nuclide indices for each of the 12 material types. The
	// grid is allocated as a full square grid, even though not all
	// materials have the same number of nuclides.
	SD.mats = load_mats(SD.num_nucs, in.n_isotopes, &SD.max_num_nucs);

	// Intialize the flattened 2D grid of nuclide concentration data. The grid holds
	// a list of nuclide concentrations for each of the 12 material types. The
	// grid is allocated as a full square grid, even though not all
	// materials have the same number of nuclides.
	SD.concs = load_concs_v(SD.num_nucs, SD.max_num_nucs);


	return SD;

}

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
