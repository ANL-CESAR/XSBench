#include "XSbench_header.h"

SimulationData grid_init_do_not_profile( Inputs in, int mype )
{
	// Structure to hold all allocated simuluation data arrays
	SimulationData SD;

	// Keep track of how much data we're allocating
	size_t nbytes = 0;

	// Set the initial seed value
	uint64_t seed = 42;	

	////////////////////////////////////////////////////////////////////
	// Initialize Nuclide Grids
	////////////////////////////////////////////////////////////////////
	
	if(mype == 0) printf("Intializing nuclide grids...\n");

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
		SD.nuclide_grid[i].energy        = LCG_random_double(&seed);
		SD.nuclide_grid[i].total_xs      = LCG_random_double(&seed);
		SD.nuclide_grid[i].elastic_xs    = LCG_random_double(&seed);
		SD.nuclide_grid[i].absorbtion_xs = LCG_random_double(&seed);
		SD.nuclide_grid[i].fission_xs    = LCG_random_double(&seed);
		SD.nuclide_grid[i].nu_fission_xs = LCG_random_double(&seed);
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
		if(mype == 0) printf("Intializing unionized grid...\n");

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
		if(mype == 0) printf("Intializing hash grid...\n");
		SD.length_unionized_energy_array = 0;
		SD.length_index_grid  = in.hash_bins * in.n_isotopes;
		SD.index_grid = (int *) malloc( SD.length_index_grid * sizeof(int)); 
		assert(SD.index_grid != NULL);
		nbytes += SD.length_index_grid * sizeof(int);

		double du = 1.0 / in.hash_bins;

		// For each energy level in the hash table
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
	if(mype == 0) printf("Intializing material data...\n");
	
	// Set the number of nuclides in each material
	SD.num_nucs  = load_num_nucs(in.n_isotopes);
	SD.length_num_nucs = 12; // There are always 12 materials in XSBench

	// Intialize the flattened 2D grid of material data. The grid holds
	// a list of nuclide indices for each of the 12 material types. The
	// grid is allocated as a full square grid, even though not all
	// materials have the same number of nuclides.
	SD.mats = load_mats(SD.num_nucs, in.n_isotopes, &SD.max_num_nucs);
	SD.length_mats = SD.length_num_nucs * SD.max_num_nucs;

	// Intialize the flattened 2D grid of nuclide concentration data. The grid holds
	// a list of nuclide concentrations for each of the 12 material types. The
	// grid is allocated as a full square grid, even though not all
	// materials have the same number of nuclides.
	SD.concs = load_concs(SD.num_nucs, SD.max_num_nucs);
	SD.length_concs = SD.length_mats;

	if(mype == 0) printf("Intialization complete. Allocated %.0lf MB of data.\n", nbytes/1024.0/1024.0 );

	return SD;

}
