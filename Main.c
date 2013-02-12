#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	
	unsigned long seed;
	size_t memtotal;
	int n_isotopes; // H-M Large is 355, H-M Small is 68
	int n_gridpoints = 900;
	int lookups = 15000000;
	int i, thread, nthreads, mat;
	double omp_start, omp_end, p_energy;
	int max_procs = omp_get_num_procs();
	char * HM;
	int bgq_mode = 0;
	

	// rand() is only used in the serial initialization stages.
	// A custom RNG is used in parallel portions.
	srand(time(NULL));

	// Process CLI Fields
	// Useage:   ./XSBench <# threads> <H-M Size ("Small or "Large")>
	if( argc == 2 )
	{
		nthreads = atoi(argv[1]);	// first arg sets # of threads
		n_isotopes = 355;			// defaults to H-M Large
	}
	else if( argc == 3 )
	{
		nthreads = atoi(argv[1]);	// first arg sets # of threads
		// second arg species small or large H-M benchmark
		if( strcmp( argv[2], "small") == 0 || strcmp( argv[2], "Small" ) == 0)
			n_isotopes = 68;
		else
			n_isotopes = 355;
	}
	else if( argc == 4 )
	{
		bgq_mode = atoi(argv[3]);  // BG/Q mode (16,8,4,2,1) 
		nthreads = atoi(argv[1]);	// first arg sets # of threads
		// second arg species small or large H-M benchmark
		if( strcmp( argv[2], "small") == 0 || strcmp( argv[2], "Small" ) == 0)
			n_isotopes = 68;
		else
			n_isotopes = 355;
	}
	else
	{
		nthreads = max_procs;		// defaults to full CPU usage
		n_isotopes = 355;			// defaults to H-M Large
	}

	// Sets H-M size name
	if( n_isotopes == 68 )
		HM = "Small";
	else
		HM = "Large";

/*	
	// Deals with BG/Q mode setting.
	// c16 = 16 ranks - we only want core  0
	// c8  = 8  ranks - we only want cores 0-1
	// c4 =  4  ranks - we only want cores 0-3
	// c2 =  2  ranks - we only want cores 0-7
	// c1 =  1  rank  - we want all cores  0-15
	#ifdef __bgq__
	switch(bgq_mode)
	{	
		case 16:
			if( get_bgq_core() != 0 )
			{
				printf("My CPU # is %d... sleeping...\n", get_bgq_core());
				return 0;
			}
			break;
		case 8:
			if( get_bgq_core() > 1 )
			{
				printf("My CPU # is %d... sleeping...\n", get_bgq_core());
				return 0;
			}
			break;
		case 4:
			if( get_bgq_core() > 3 )
			{
				printf("My CPU # is %d... sleeping...\n", get_bgq_core());
				return 0;
			}
			break;
		case 2:
			if( get_bgq_core() > 7 )
			{
				printf("My CPU # is %d... sleeping...\n", get_bgq_core());
				return 0;
			}
			break;
		default:
			break;
	}
	printf("My CPU # is %d... Running test...\n", get_bgq_core());
	#endif
*/		
	// Set number of OpenMP Threads
	omp_set_num_threads(nthreads); 
		
	// calculate estimate for memory useage
	size_t single_nuclide_grid = n_gridpoints * sizeof( NuclideGridPoint );
	size_t all_nuclide_grids = n_isotopes * single_nuclide_grid;
	size_t size_GridPoint =sizeof(GridPoint)+n_isotopes*sizeof(NuclideGridPoint *);
	size_t size_UEG = n_isotopes*n_gridpoints * size_GridPoint;
	memtotal = all_nuclide_grids + size_UEG;
	
	all_nuclide_grids = all_nuclide_grids  / 1048576;
	size_UEG = size_UEG / 1048576;
	memtotal = memtotal / 1048576;

	// =====================================================================
	// Print-out of Input Summary
	// =====================================================================
	
	logo();
	center_print("INPUT SUMMARY", 79);
	printf(
	"###################################################################"
	"#############\n"
	);
	printf("Materials:                    %d\n", 12);
	printf("H-M Benchmark Size:           %s\n", HM);
	printf("Total Isotopes:               %d\n", n_isotopes);
	printf("Gridpoints (per Nuclide):     %d\n", n_gridpoints);
	printf("Unionized Energy Gridpoints:  %d\n", n_isotopes*n_gridpoints);
	printf("XS Lookups:                   %d\n", lookups);
	printf("Threads:                      %d\n", nthreads);
	printf("Nuclide Grid Memory:          %zu\n", all_nuclide_grids);
	printf("UEG Memory:                   %zu\n", size_UEG);
	printf("Est. Memory Usage (MB):       %zu\n", memtotal);
	if( EXTRA_FLOPS > 0 )
		printf("Extra Flops:                  %d\n", EXTRA_FLOPS);
	if( EXTRA_LOADS > 0 )
		printf("Extra Loads:                  %d\n", EXTRA_LOADS);
	printf(
	"###################################################################"
	"#############\n");
	center_print("INITIALIZATION", 79);
	printf(
	"###################################################################"
	"#############\n"
	);
	
	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	printf("Generating Nuclide Energy Grids...\n");
	
	NuclideGridPoint ** nuclide_grids = gpmatrix( n_isotopes, n_gridpoints );
	
	generate_grids( nuclide_grids, n_isotopes, n_gridpoints );	
	
	// Sort grids by energy
	sort_nuclide_grids( nuclide_grids, n_isotopes, n_gridpoints );

	// Prepare Unionized Energy Grid Framework
	GridPoint * energy_grid = generate_energy_grid( n_isotopes, n_gridpoints,
	                                                nuclide_grids ); 	

	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	set_grid_ptrs( energy_grid, nuclide_grids, n_isotopes, n_gridpoints );
	
	// Get material data
	printf("Loading Mats...\n");
	int *num_nucs = load_num_nucs(n_isotopes);
	int **mats = load_mats(num_nucs, n_isotopes);
	double **concs = load_concs(num_nucs);

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================
	
	printf(
	"###################################################################"
	"#############\n");
	center_print("SIMULATION", 79);
	printf(
	"###################################################################"
	"#############\n"
	);

	omp_start = omp_get_wtime();
	
	#ifdef __PAPI
	int eventset = PAPI_NULL; 
	int num_papi_events;
	counter_init(&eventset, &num_papi_events);
	#endif
	
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat, seed) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs)
	{	
		double macro_xs_vector[5];
		thread = omp_get_thread_num();
		seed = (thread+1)*19+17;
		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( INFO && thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);
			
			// Randomly pick an energy and material for the particle
			p_energy = rn(&seed);
			mat = pick_mat(&seed); 
		
			// This returns the macro_xs_vector, but we're not going
			// to do anything with it in this program, so return value
			// is written over.
			calculate_macro_xs( p_energy, mat, n_isotopes,
			                    n_gridpoints, num_nucs, concs,
			                    energy_grid, nuclide_grids, mats,
                                macro_xs_vector );
		}
	}
	
	printf("\n" );
	printf("Simulation complete.\n" );

	omp_end = omp_get_wtime();
	
	// =====================================================================
	// Print / Save Results and Exit
	// =====================================================================
	
	printf(
	"###################################################################"
	"#############\n");
	center_print("RESULTS", 79);
	printf(
	"###################################################################"
	"#############\n"
	);

	// Print the results
	printf("Threads:     %d\n", nthreads);
	if( EXTRA_FLOPS > 0 )
	printf("Extra Flops: %d\n", EXTRA_FLOPS);
	if( EXTRA_LOADS > 0 )
	printf("Extra Loads: %d\n", EXTRA_LOADS);
	printf("Runtime:     %.3lf seconds\n", omp_end-omp_start);
	printf("Lookups:     %d\n", lookups);
	printf("Lookups/s:   %.0lf\n", (double) lookups / (omp_end-omp_start));
	printf(
	"###################################################################"
	"#############\n"
	);

	// For bechmarking, output lookup/s data to file
	if( SAVE )
	{
		FILE * out = fopen( "results.txt", "a" );
		fprintf(out, "c%d\t%d\t%.0lf\n", bgq_mode, nthreads, (double) lookups / (omp_end-omp_start));
		fclose(out);
	}
	
	#ifdef __PAPI
	counter_stop(&eventset, num_papi_events);
	#endif

	return 0;
}
