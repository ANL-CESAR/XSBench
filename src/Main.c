#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	
	int version = 11;
	int mype = 0;
	int max_procs = omp_get_num_procs();
	int n_isotopes, n_gridpoints, lookups, i, thread, nthreads, mat;
	unsigned long seed;
	double omp_start, omp_end, p_energy;
	char * HM;
	unsigned long long vhash = 0;

	#ifdef MPI
	int nprocs;
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

	// Process CLI Fields
	Inputs input = read_CLI( argc, argv );
	
	// Set CLI variables
	nthreads =     input.nthreads;
	n_isotopes =   input.n_isotopes;
	n_gridpoints = input.n_gridpoints;
	lookups =      input.lookups;
	HM =           input.HM;

	// Set number of OpenMP Threads
	omp_set_num_threads(nthreads); 
		
	// =====================================================================
	// Calculate Estimate of Memory Usage
	// =====================================================================

	size_t single_nuclide_grid = n_gridpoints * sizeof( NuclideGridPoint );
	size_t all_nuclide_grids   = n_isotopes * single_nuclide_grid;
	size_t size_GridPoint      = sizeof(GridPoint) + n_isotopes*sizeof(int);
	size_t size_UEG            = n_isotopes * n_gridpoints * size_GridPoint;
	size_t memtotal;
	int mem_tot;

	memtotal          = all_nuclide_grids + size_UEG;
	all_nuclide_grids = all_nuclide_grids / 1048576;
	size_UEG          = size_UEG / 1048576;
	memtotal          = memtotal / 1048576;
	mem_tot           = memtotal;

	// =====================================================================
	// Print-out of Input Summary
	// =====================================================================
	
	if( mype == 0 )
	{
		logo(version);
		center_print("INPUT SUMMARY", 79);
		border_print();
		#ifdef VERIFICATION
		printf("Verification Mode:            on\n");
		#endif
		printf("Materials:                    %d\n", 12);
		printf("H-M Benchmark Size:           %s\n", HM);
		printf("Total Nuclides:               %d\n", n_isotopes);
		printf("Gridpoints (per Nuclide):     ");
		fancy_int(n_gridpoints);
		printf("Unionized Energy Gridpoints:  ");
		fancy_int(n_isotopes*n_gridpoints);
		printf("XS Lookups:                   "); fancy_int(lookups);
		#ifdef MPI
		printf("MPI Ranks:                    %d\n", nprocs);
		printf("OMP Threads per MPI Rank:     %d\n", nthreads);
		printf("Mem Usage per MPI Rank (MB):  "); fancy_int(mem_tot);
		#else
		printf("Threads:                      %d\n", nthreads);
		printf("Est. Memory Usage (MB):       "); fancy_int(mem_tot);
		#endif
		if( EXTRA_FLOPS > 0 )
			printf("Extra Flops:                  %d\n", EXTRA_FLOPS);
		if( EXTRA_LOADS > 0 )
			printf("Extra Loads:                  %d\n", EXTRA_LOADS);
		border_print();
		center_print("INITIALIZATION", 79);
		border_print();
	}

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	if( mype == 0) printf("Generating Nuclide Energy Grids...\n");
	
	NuclideGridPoint ** nuclide_grids = gpmatrix( n_isotopes, n_gridpoints );
	
	#ifdef VERIFICATION
	generate_grids_v( nuclide_grids, n_isotopes, n_gridpoints );	
	#else
	generate_grids( nuclide_grids, n_isotopes, n_gridpoints );	
	#endif
	
	// Sort grids by energy
	if( mype == 0) printf("Sorting Nuclide Energy Grids...\n");
	sort_nuclide_grids( nuclide_grids, n_isotopes, n_gridpoints );

	// Prepare Unionized Energy Grid Framework
	GridPoint * energy_grid = generate_energy_grid( n_isotopes, n_gridpoints,
	                                                nuclide_grids ); 	

	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	set_grid_ptrs( energy_grid, nuclide_grids, n_isotopes, n_gridpoints );
	
	// Get material data
	if( mype == 0 ) printf("Loading Mats...\n");
	int *num_nucs  = load_num_nucs(n_isotopes);
	int **mats     = load_mats(num_nucs, n_isotopes);
	#ifdef VERIFICATION
	double **concs = load_concs_v(num_nucs);
	#else
	double **concs = load_concs(num_nucs);
	#endif

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	#ifdef BENCHMARK
	for( int jrt = 1; jrt <=omp_get_num_procs(); jrt++ )
	{
		nthreads = jrt;
		omp_set_num_threads(nthreads);
 	#endif

	if( mype == 0 )
	{
		printf("\n");
		border_print();
		center_print("SIMULATION", 79);
		border_print();
	}

	omp_start = omp_get_wtime();
  
	#ifdef PAPI
	/* initialize papi with one thread here  */
	if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	#endif	

	// OpenMP compiler directives - declaring variables as shared or private
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat, seed) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs, mype, vhash) 
	{	
		#ifdef PAPI
		int eventset = PAPI_NULL; 
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events);
		}
		#endif

		double macro_xs_vector[5];
		thread = omp_get_thread_num();
		seed   = (thread+1)*19+17;

		#pragma omp for schedule(dynamic)
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( INFO && mype == 0 && thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)lookups / (double) nthreads )) / (double) nthreads * 100.0);
			// Randomly pick an energy and material for the particle
			#ifdef VERIFICATION
			#pragma omp critical
			{
				p_energy = rn_v();
				mat      = pick_mat(&seed); 
			}
			#else
			p_energy = rn(&seed);
			mat      = pick_mat(&seed); 
			#endif
			
			// debugging
			//printf("E = %lf mat = %d\n", p_energy, mat);
				
			// This returns the macro_xs_vector, but we're not going
			// to do anything with it in this program, so return value
			// is written over.
			calculate_macro_xs( p_energy, mat, n_isotopes,
			                    n_gridpoints, num_nucs, concs,
			                    energy_grid, nuclide_grids, mats,
                                macro_xs_vector );

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

			// Artificial pause injected to represent particle
			// tracking calculation time. Similar to adding dummy flops.
			#ifdef PAUSE
			struct timespec ts, rts;
			ts.tv_sec = 0;
			//ts.tv_nsec = 100000; // .1 ms
			ts.tv_nsec = 1000000; // 1 ms
			nanosleep(&ts, &rts);
			#endif
		}

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
	
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	
	// Calculate Lookups per sec
	int lookups_per_sec = (int) ((double) lookups / (omp_end-omp_start));
	
	// If running in MPI, reduce timing statistics and calculate average
	#ifdef MPI
	int total_lookups = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&lookups_per_sec, &total_lookups, 1, MPI_INT,
	           MPI_SUM, 0, MPI_COMM_WORLD);
	//total_lookups = total_lookups / nprocs;
	#endif
	
	// Print output
	if( mype == 0 )
	{
		border_print();
		center_print("RESULTS", 79);
		border_print();

		// Print the results
		printf("Threads:     %d\n", nthreads);
		#ifdef MPI
		printf("MPI ranks:   %d\n", nprocs);
		#endif
		if( EXTRA_FLOPS > 0 )
		printf("Extra Flops: %d\n", EXTRA_FLOPS);
		if( EXTRA_LOADS > 0 )
		printf("Extra Loads: %d\n", EXTRA_LOADS);
		#ifdef MPI
		printf("Total Lookups/s:            ");
		fancy_int(total_lookups);
		printf("Avg Lookups/s per MPI rank: ");
		fancy_int(total_lookups / nprocs);
		#else
		printf("Runtime:     %.3lf seconds\n", omp_end-omp_start);
		printf("Lookups:     "); fancy_int(lookups);
		printf("Lookups/s:   ");
		fancy_int(lookups_per_sec);
		#endif
		#ifdef VERIFICATION
		printf("Verification checksum: %llu\n", vhash);
		#endif
		border_print();

		// For bechmarking, output lookup/s data to file
		if( SAVE )
		{
			FILE * out = fopen( "results.txt", "a" );
			fprintf(out, "%d\t%d\n", nthreads, lookups_per_sec);
			fclose(out);
		}
	}	
	#ifdef BENCHMARK
	}
	#endif

	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
