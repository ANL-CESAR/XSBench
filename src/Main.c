#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	int version = 17;
	int mype = 0;
	int thread;
	double omp_start, omp_end;
	unsigned long long vhash = 0;
	int nprocs = 1;

	#ifdef MPI
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	// rand() is only used in the serial initialization stages.
	// A custom RNG is used in parallel portions.
	#ifdef VERIFICATION
	srand(26);
	#else
	srand(time(NULL));
	#endif

	// Process CLI Fields -- store in "Inputs" structure
	Inputs in = read_CLI( argc, argv );
	
	// Set number of OpenMP Threads
	omp_set_num_threads(in.nthreads); 

	// Print-out of Input Summary
	if( mype == 0 )
		print_inputs( in, nprocs, version );

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	#ifndef BINARY_READ
	if( mype == 0) printf("Generating Nuclide Energy Grids...\n");
	#endif
	
	NuclideGridPoint ** nuclide_grids = gpmatrix(in.n_isotopes,in.n_gridpoints);
	
	#ifdef VERIFICATION
	generate_grids_v( nuclide_grids, in.n_isotopes, in.n_gridpoints );	
	#else
	generate_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );	
	#endif
	
	// Sort grids by energy
	#ifndef BINARY_READ
	if( mype == 0) printf("Sorting Nuclide Energy Grids...\n");
	sort_nuclide_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );
	#endif

	// If using a unionized grid search, initialize the energy grid
	// Otherwise, leave these as null
	GridPoint * energy_grid = NULL;
	int * index_data = NULL;

	if( in.grid_type == UNIONIZED )
	{
		// Prepare Unionized Energy Grid Framework
		#ifndef BINARY_READ
		energy_grid = generate_energy_grid( in.n_isotopes,
											in.n_gridpoints, nuclide_grids ); 	
		#else
		energy_grid = (GridPoint *)malloc( in.n_isotopes *
										   in.n_gridpoints * sizeof( GridPoint ) );
		index_data = (int *) malloc( in.n_isotopes * in.n_gridpoints
										   * in.n_isotopes * sizeof(int));
		for( int i = 0; i < in.n_isotopes*in.n_gridpoints; i++ )
			energy_grid[i].xs_ptrs = &index_data[i*in.n_isotopes];
		#endif

		// Double Indexing. Filling in energy_grid with pointers to the
		// nuclide_energy_grids.
		#ifndef BINARY_READ
		initialization_do_not_profile_set_grid_ptrs( energy_grid, nuclide_grids, in.n_isotopes, in.n_gridpoints );
		#endif
	}
	else if( in.grid_type == HASH )
	{
		energy_grid = generate_hash_table( nuclide_grids, in.n_isotopes, in.n_gridpoints, in.hash_bins );
	}

	#ifdef BINARY_READ
	if( mype == 0 ) printf("Reading data from \"XS_data.dat\" file...\n");
	binary_read(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, in.grid_type);
	#endif
	
	// Get material data
	if( mype == 0 )
		printf("Loading Mats...\n");
	int *num_nucs  = load_num_nucs(in.n_isotopes);
	int **mats     = load_mats(num_nucs, in.n_isotopes);

	#ifdef VERIFICATION
	double **concs = load_concs_v(num_nucs);
	#else
	double **concs = load_concs(num_nucs);
	#endif

	#ifdef BINARY_DUMP
	if( mype == 0 ) printf("Dumping data to binary file...\n");
	binary_dump(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, in.grid_type);
	if( mype == 0 ) printf("Binary file \"XS_data.dat\" written! Exiting...\n");
	return 0;
	#endif

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================

	// Outer benchmark loop can loop through all possible # of threads
	#ifdef BENCHMARK
	for( int bench_n = 1; bench_n <=omp_get_num_procs(); bench_n++ )
	{
		in.nthreads = bench_n;
		omp_set_num_threads(in.nthreads);
 	#endif

	if( mype == 0 )
	{
		printf("\n");
		border_print();
		center_print("SIMULATION", 79);
		border_print();
	}

	omp_start = omp_get_wtime();
  
	//initialize papi with one thread (master) here
	#ifdef PAPI
	if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	#endif	

	// OpenMP compiler directives - declaring variables as shared or private
	#pragma omp parallel default(none) \
	private(thread) \
	shared( in, energy_grid, nuclide_grids, \
	        mats, concs, num_nucs, mype, vhash) 
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

		// Initialize RNG seeds for threads
		thread = omp_get_thread_num();

		// Particle loop
		// (independent - can be processed in any order and in parallel)
		#pragma omp for schedule(dynamic)
		for( int p = 0; p < in.particles; p++ )
		{
			// Particles are seeded by their particle ID
			unsigned long seed = (p+1)*13371337;

			// Randomly pick an energy and material for the particle
			double p_energy = rn(&seed);
			int mat      = pick_mat(&seed); 

			// Status text
			if( INFO && mype == 0 && thread == 0 && p % 100 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(p / ( (double)in.particles / (double) in.nthreads ))
						/ (double) in.nthreads * 100.0);

			// XS Lookup Loop
			// (dependent! Next iteration uses data computed in previous iter.)
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
				#pragma omp atomic
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

	#ifndef PAPI
	if( mype == 0)	
	{	
		printf("\n" );
		printf("Simulation complete.\n" );
	}
	#endif

	omp_end = omp_get_wtime();

	// Final Hash Step
	vhash = vhash % 1000000;
	
	// Print / Save Results and Exit
	print_results( in, mype, omp_end-omp_start, nprocs, vhash );

	#ifdef BENCHMARK
	}
	#endif

	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
