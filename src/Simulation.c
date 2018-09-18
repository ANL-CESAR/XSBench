#include "XSbench_header.h"
void run_event_based_simulation(Inputs in, GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids, int * num_nucs, int ** mats, double ** concs, int mype, unsigned long long * vhash_result)
{
	if( mype == 0)	
		printf("Beginning event based simulation...\n");

	unsigned long long vhash = 0;
	// OpenMP compiler directives - declaring variables as shared or private
	// The reduction is only needed when in verification mode.
	#pragma omp parallel default(none) \
	shared( in, energy_grid, nuclide_grids, \
			mats, concs, num_nucs, mype) \
	reduction(+:vhash)
	{	
		// Initialize parallel PAPI counters
		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif

		double * xs = (double *) calloc(5, sizeof(double));

		// Initialize RNG seeds for threads
		int thread = omp_get_thread_num();

		// XS Lookup Loop
		// This loop is independent. Represents lookup events for many particles executed independently in one loop.
		//     i.e., All iterations can be processed in any order and are not related
		#pragma omp for schedule(guided)
		for( int i = 0; i < in.lookups; i++ )
		{
			// Status text
			if( INFO && mype == 0 && thread == 0 && i % 2000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double) in.lookups / (double) in.nthreads ))
						/ (double) in.nthreads * 100.0);
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
					in.n_gridpoints, num_nucs, concs,
					energy_grid, nuclide_grids, mats,
					macro_xs_vector, in.grid_type, in.hash_bins );

			// Copy results from above function call onto heap
			// so that compiler cannot optimize function out
			// (only occurs if -flto flag is used)
			// This operation is only done to avoid optimizing out
			// calculate_macro_xs -- we do not care about what is
			// in the "xs" array
			memcpy(xs, macro_xs_vector, 5*sizeof(double));

			// Verification hash calculation
			// This method provides a consistent hash accross
			// architectures and compilers.
			#ifdef VERIFICATION
			char line[256];
			sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
					p_energy, mat,
					macro_xs_vector[0],
					macro_xs_vector[1],
					macro_xs_vector[2],
					macro_xs_vector[3],
					macro_xs_vector[4]);
			unsigned long long vhash_local = hash(line, 10000);

			vhash += vhash_local;
			#endif
		}

		// Prints out thread local PAPI counters
		#ifdef PAPI
		if( mype == 0 && thread == 0 )
		{
			printf("\n");
			border_print();
			center_print("PAPI COUNTER RESULTS", 79);
			border_print();
			printf("Count          \tSmybol      \tDescription\n");
		}
		{
			#pragma omp barrier
		}
		counter_stop(&eventset, num_papi_events);
		#endif

	}
	*vhash_result = vhash;
}

void run_history_based_simulation(Inputs in, GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids, int * num_nucs, int ** mats, double ** concs, int mype, unsigned long long * vhash_result)
{
	if( mype == 0)	
		printf("Beginning history based simulation...\n");
	unsigned long long vhash = 0;
	// OpenMP compiler directives - declaring variables as shared or private
	// The reduction is only needed when in verification mode.
	#pragma omp parallel default(none) \
	shared( in, energy_grid, nuclide_grids, \
			mats, concs, num_nucs, mype) \
	reduction(+:vhash)
	{	
		// Initialize parallel PAPI counters
		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif

		double * xs = (double *) calloc(5, sizeof(double));

		// Initialize RNG seeds for threads
		int thread = omp_get_thread_num();

		// Particle loop 
		// (independent - can be processed in any order and in parallel)
		// Only present in History based method (default)
		#pragma omp for schedule(guided)
		for( int p = 0; p < in.particles; p++ )
		{
			// Particles are seeded by their particle ID
			unsigned long seed = ((unsigned long) p+ (unsigned long)1)* (unsigned long) 13371337;

			// Randomly pick an energy and material for the particle
			double p_energy = rn(&seed);
			int mat      = pick_mat(&seed); 

			// Status text
			if( INFO && mype == 0 && thread == 0 && p % 100 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(p / ( (double)in.particles / (double) in.nthreads ))
						/ (double) in.nthreads * 100.0);

			// XS Lookup Loop
			// This loop is dependent!
			// i.e., Next iteration uses data computed in previous iter.
			for( int i = 0; i < in.lookups; i++ )
			{
				// debugging
				//printf("E = %lf mat = %d\n", p_energy, mat);

				double macro_xs_vector[5] = {0};

				// This returns the macro_xs_vector, but we're not going
				// to do anything with it in this program, so return value
				// is written over.
				calculate_macro_xs( p_energy, mat, in.n_isotopes,
						in.n_gridpoints, num_nucs, concs,
						energy_grid, nuclide_grids, mats,
						macro_xs_vector, in.grid_type, in.hash_bins );

				// Copy results from above function call onto heap
				// so that compiler cannot optimize function out
				// (only occurs if -flto flag is used)
				// This operation is only done to avoid optimizing out
				// calculate_macro_xs -- we do not care about what is
				// in the "xs" array
				memcpy(xs, macro_xs_vector, 5*sizeof(double));

				// Verification hash calculation
				// This method provides a consistent hash accross
				// architectures and compilers.
				#ifdef VERIFICATION
				char line[256];
				sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
						p_energy, mat,
						macro_xs_vector[0],
						macro_xs_vector[1],
						macro_xs_vector[2],
						macro_xs_vector[3],
						macro_xs_vector[4]);
				unsigned long long vhash_local = hash(line, 10000);

				vhash += vhash_local;
				#endif

				// Randomly pick next energy and material for the particle
				// Also incorporates results from macro_xs lookup to
				// enforce loop dependency.
				// In a real MC app, this dependency is expressed in terms
				// of branching physics sampling, whereas here we are just
				// artificially enforcing this dependence based on altering
				// the seed
				for( int x = 0; x < 5; x++ )
					seed += macro_xs_vector[x] * (x+1)*1337*1337;

				p_energy = rn(&seed);
				mat      = pick_mat(&seed); 
			}
		}

		// Prints out thread local PAPI counters
		#ifdef PAPI
		if( mype == 0 && thread == 0 )
		{
			printf("\n");
			border_print();
			center_print("PAPI COUNTER RESULTS", 79);
			border_print();
			printf("Count          \tSmybol      \tDescription\n");
		}
		{
			#pragma omp barrier
		}
		counter_stop(&eventset, num_papi_events);
		#endif

	}
	*vhash_result = vhash;
}
