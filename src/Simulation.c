#include "XSbench_header.h"

void run_event_based_simulation(Inputs in, GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids, int * num_nucs, int ** mats, double ** concs, int mype, unsigned long long * vhash_result)
{
	if( mype == 0)	
		printf("Beginning event based simulation...\n");
	
	////////////////////////////////////////////////////////////////////////////////
	// First, we convert our various 2D arrays (some of which are jagged)
	// to 1D for easier target transfers
	////////////////////////////////////////////////////////////////////////////////
	
	// Figure out some sizing parameters, based on type of lookup method used
	int length_energy_grid;
	int length_index_grid;
	if( in.grid_type == UNIONIZED )
	{
		length_energy_grid = in.n_isotopes * in.n_gridpoints;
		length_index_grid  = length_energy_grid * in.n_isotopes;
	}
	else if( in.grid_type == HASH )
	{
		length_energy_grid = in.hash_bins;
		length_index_grid  = length_energy_grid * in.n_isotopes;
	}
	else
	{
		length_energy_grid = 0;
		length_index_grid  = 0;
	}

	// Convert Energy Grid Compount structure to separate egrid and index 1D arrays
	// (Not necessary if only using the nuclide grid)
	double * egrid = NULL;
	int * index_data = NULL;

	if( in.grid_type == UNIONIZED || in.grid_type == HASH )
	{
		egrid = (double *) malloc( length_energy_grid * sizeof(double));
		index_data =  (int *) malloc( length_index_grid  * sizeof(int));
		for( int gp = 0; gp < length_energy_grid; gp++ )
		{
			egrid[gp] = energy_grid[gp].energy;
			for( int n = 0; n < in.n_isotopes; n++ )
			{
				index_data[gp * in.n_isotopes + n] = energy_grid[gp].xs_ptrs[n];
			}
		}
	}

	// Convert 2D nuclide grid array to 1D. This was a contiguously
	// allocated 2D array so just need to grab the data pointer.
	NuclideGridPoint * nuclide_grids_t = nuclide_grids[0];
	int length_nuclide_grids   = in.n_isotopes * in.n_gridpoints;

	// Convert the materials and concentrations (jagged) 2D arrays to 1D.
	// Since they're jagged, need to figure out maximum element and allocate
	// full size 1D array with extra delements.
	
	// Determine maximum number of nuclides for all materials
	int num_mats = 12;
	int max_num_nucs = 0;
	for( int m = 0; m < num_mats; m++ )
	{
		if( num_nucs[m] > max_num_nucs )
			max_num_nucs = num_nucs[m];
	}
	
	int length_mats  = num_mats * max_num_nucs;
	int length_concs = length_mats;

	// Copy materials and concs jagged 2D vectors over to 1D
	int * mats_t     = (int *)    calloc( length_mats,  sizeof(int)); 
	double * concs_t = (double *) calloc( length_concs, sizeof(double));

	for( int m = 0; m < num_mats; m++ )
	{
		for( int n = 0; n < num_nucs[m]; n++ ) 
		{
			mats_t[m*max_num_nucs + n] = mats[m][n];
			concs_t[m*max_num_nucs + n] = concs[m][n];
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////
	// SUMMARY: Data Structure Manifest
	// Heap arrays (and lengths) that would need to be offloaded to an accelerator
	////////////////////////////////////////////////////////////////////////////////
	// int * num_nucs   :::: Length = 12
	// double * concs_t :::: Length = length_concs
	// int * mats_t     :::: Length = length_mats
	// double * egrid   :::: Length = length_energy_grid
	// int * index_data :::: Length = length_index_grid
	// NuclideGridPoint * nuclide_grids_t :::: Length = length_nuclide_grids
	// 
	// Note: "egrid" and "index_data" could be of length 0, if nuclide grid only
	// method was selected by user, i.e., if in.grid_type == NUCLIDE
	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	// Begin Actual Simulation Loop 
	////////////////////////////////////////////////////////////////////////////////
	unsigned long long verification = 0;
	#pragma omp parallel for schedule(guided) reduction(+:verification)
	for( int i = 0; i < in.lookups; i++ )
	{
		// Particles are seeded by their particle ID
		unsigned long seed = ((unsigned long) i+ (unsigned long)1)* (unsigned long) 13371337;

		// Randomly pick an energy and material for the particle
		double p_energy = rn(&seed);
		int mat      = pick_mat(&seed); 

		// debugging
		//printf("E = %lf mat = %d\n", p_energy, mat);

		double macro_xs_vector[5] = {0};

		// This returns the macro_xs_vector, but we're not going
		// to do anything with it in this program, so return value
		// is written over.
	
		calculate_macro_xs( p_energy, mat, in.n_isotopes,
			in.n_gridpoints, num_nucs, concs_t,
			egrid, index_data, nuclide_grids_t, mats_t,
			macro_xs_vector, in.grid_type, in.hash_bins, max_num_nucs );

		// For verification, and to prevent the compiler from optimizing
		// all work out, we interrogate the returned macro_xs_vector array
		// to find its maximum value index, then increment the verification
		// value by that index. In this implementation, we prevent thread
		// contention by using an OMP reduction on it. For other accelerators,
		// a different approach might be required (e.g., atomics, reduction
		// of thread-specific values in large array via CUDA thrust, etc)
		double max = -1.0;
		int max_idx = 0;
		for(int i = 0; i < 5; i++ )
		{
			if( macro_xs_vector[i] > max )
			{
				max = macro_xs_vector[i];
				max_idx = i;
			}
		}
		verification += max_idx;
	}
	*vhash_result = verification;
}
