#include "XSbench_header.h"

unsigned long long run_event_based_simulation(Inputs in, SimulationData SD, int mype)
{
	if( mype == 0)	
		printf("Beginning event based simulation...\n");
	
	////////////////////////////////////////////////////////////////////////////////
	// SUMMARY: Simulation Data Structure Manifest for "SD" Object
	// Here we list all heap arrays (and lengths) in SD that would need to be
	// offloaded manually if using an accelerator with a seperate memory space
	////////////////////////////////////////////////////////////////////////////////
	// int * num_nucs;                     // Length = length_num_nucs;
	// double * concs;                     // Length = length_concs
	// int * mats;                         // Length = length_mats
	// double * unionized_energy_array;    // Length = length_unionized_energy_array
	// int * index_grid;                   // Length = length_index_grid
	// NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
	// 
	// Note: "unionized_energy_array" and "index_grid" can be of zero length
	//        depending on lookup method.
	//
	// Note: "Lengths" are given as the number of objects in the array, not the
	//       number of bytes.
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
				mat,             // Sampled material type index neutron is in
				in.n_isotopes,   // Total number of isotopes in simulation
				in.n_gridpoints, // Number of gridpoints per isotope in simulation
				SD.num_nucs,     // 1-D array with number of nuclides per material
				SD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
				SD.unionized_energy_array, // 1-D Unionized energy array
				SD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
				SD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
				SD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
				macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
				in.grid_type,    // Lookup type (nuclide, hash, or unionized)
				in.hash_bins,    // Number of hash bins used (if using hash lookup type)
				SD.max_num_nucs  // Maximum number of nuclides present in any material
				);

		// For verification, and to prevent the compiler from optimizing
		// all work out, we interrogate the returned macro_xs_vector array
		// to find its maximum value index, then increment the verification
		// value by that index. In this implementation, we prevent thread
		// contention by using an OMP reduction on the verification value.
		// For accelerators, a different approach might be required
		// (e.g., atomics, reduction of thread-specific values in large
		// array via CUDA thrust, etc).
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

unsigned long long run_history_based_simulation(Inputs in, SimulationData SD, int mype)
{
	if( mype == 0)	
		printf("Beginning history based simulation...\n");

	
	////////////////////////////////////////////////////////////////////////////////
	// SUMMARY: Simulation Data Structure Manifest for "SD" Object
	// Here we list all heap arrays (and lengths) in SD that would need to be
	// offloaded manually if using an accelerator with a seperate memory space
	////////////////////////////////////////////////////////////////////////////////
	// int * num_nucs;                     // Length = length_num_nucs;
	// double * concs;                     // Length = length_concs
	// int * mats;                         // Length = length_mats
	// double * unionized_energy_array;    // Length = length_unionized_energy_array
	// int * index_grid;                   // Length = length_index_grid
	// NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
	// 
	// Note: "unionized_energy_array" and "index_grid" can be of zero length
	//        depending on lookup method.
	//
	// Note: "Lengths" are given as the number of objects in the array, not the
	//       number of bytes.
	////////////////////////////////////////////////////////////////////////////////

	unsigned long long verification = 0;

	// Begin outer lookup loop over particles. This loop is independent.
	#pragma omp parallel for schedule(guided) reduction(+:verification)
	for( int p = 0; p < in.particles; p++ )
	{
		// Particles are seeded by their particle ID
		unsigned long seed = ((unsigned long) p+ (unsigned long)1)* (unsigned long) 13371337;

		// Randomly pick an energy and material for the particle
		double p_energy = rn(&seed);
		int mat      = pick_mat(&seed); 

		// Inner XS Lookup Loop
		// This loop is dependent!
		// i.e., Next iteration uses data computed in previous iter.
		for( int i = 0; i < in.lookups; i++ )
		{
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
			for(int j = 0; j < 5; j++ )
			{
				if( macro_xs_vector[j] > max )
				{
					max = macro_xs_vector[j];
					max_idx = j;
				}
			}
			verification += max_idx;

			// Randomly pick next energy and material for the particle
			// Also incorporates results from macro_xs lookup to
			// enforce loop dependency.
			// In a real MC app, this dependency is expressed in terms
			// of branching physics sampling, whereas here we are just
			// artificially enforcing this dependence based on altering
			// the seed
			for( int j = 0; j < 5; j++ )
				seed += macro_xs_vector[j] * (j+1)*1337*1337;

			p_energy = rn(&seed);
			mat      = pick_mat(&seed); 
		}

	}
	return verification;
}
