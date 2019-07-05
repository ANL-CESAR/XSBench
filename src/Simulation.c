#include "XSbench_header.h"

unsigned long long run_event_based_simulation(Inputs in, SimulationData SD, int mype)
{
	if( mype == 0)	
		printf("Beginning event based simulation...\n");
	
	////////////////////////////////////////////////////////////////////////////////
	// SUMMARY: Simulation Data Structure Manifest of "SD" Object
	// Heap arrays (and lengths) that would need to be offloaded to an accelerator
	////////////////////////////////////////////////////////////////////////////////
	// int * num_nucs;                     // Length = length_num_nucs;
	// double * concs;                     // Length = length_concs
	// int * mats;                         // Length = length_mats
	// double * unionized_energy_array;    // Length = length_unionized_energy_array
	// int * index_grid;                   // Length = length_index_grid
	// NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
	 
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

		// Perform macroscopic Cross Section Lookup
		calculate_macro_xs(
				p_energy,        // Sampled neutron energy (in lethargy)
				mat,             // Sampled material type neutron is in
				in.n_isotopes,   // Total number of isotopes in simulation
				in.n_gridpoints, // Number of gridpoints per isotope in simulation
				SD.num_nucs,     // 1-D array with number of nuclides per material
				SD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
				SD.unionized_energy_array, // 1-D Unionized energy array
				SD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
				SD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
				SD.mats,         // Flattened 2-D array with nuclide indices for each type of material
				macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
				in.grid_type,    // Lookup type (nuclide, hash, or unionized)
				in.hash_bins,    // Number of hash bins used (if using hash lookups)
				SD.max_num_nucs  // Maximum number of nuclides present in any material
				);

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

	return verification;
}
