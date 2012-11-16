#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	unsigned long seed;
	int n_isotopes = 68;
	int n_gridpoints = 10000;
	int lookups = 100000000;
	int i, thread, nthreads, mat;
	double omp_start, omp_end, p_energy;
	int max_procs = omp_get_num_procs();

	#ifdef __PAPI
	int eventset = PAPI_NULL; 
	counter_init(&eventset);
	#endif
	
	// rand() is only used in the serial initialization stages.
	// A custom RNG is used in parallel portions.
	srand(time(NULL));

	if( argc == 2 )
		nthreads = atoi(argv[1]);
	else
		nthreads = max_procs;
	
	omp_set_num_threads(nthreads); 

	logo();
	
	// Print out input summary
	center_print("INPUT SUMMARY", 79);
	printf(
	"###################################################################"
	"#############\n"
	);
	printf("Materials:                    %d\n", 12);
	printf("Total Isotopes:               %d\n", n_isotopes);
	printf("Gridpoints (per Nuclide):     %d\n", n_gridpoints);
	printf("Unionized Energy Gridpoints:  %d\n", n_isotopes*n_gridpoints);
	printf("XS Lookups:                   %d\n", lookups);
	printf("Threads:                      %d\n", nthreads);
	printf(
	"###################################################################"
	"#############\n");
	center_print("INITIALIZATION", 79);
	printf(
	"###################################################################"
	"#############\n"
	);

	// Allocate & fill energy grids
	printf("Generating Nuclide Energy Grids...\n");
	
	NuclideGridPoint ** nuclide_grids = gpmatrix( n_isotopes, n_gridpoints );
	generate_grids( nuclide_grids, n_isotopes, n_gridpoints );	
	
	// Sort grids by energy
	sort_nuclide_grids( nuclide_grids, n_isotopes );

	// Prepare Unionized Energy Grid Framework
	GridPoint * energy_grid = generate_energy_grid( n_isotopes, n_gridpoints,
	                                                nuclide_grids ); 	

	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	set_grid_ptrs( energy_grid, nuclide_grids, n_isotopes, n_gridpoints );
	
	// Get material data
	printf("Loading Mats...\n");
	int *num_nucs = load_num_nucs();
	int **mats = load_mats(num_nucs);
	double **concs = load_concs(num_nucs);

	printf(
	"###################################################################"
	"#############\n");
	center_print("SIMULATION", 79);
	printf(
	"###################################################################"
	"#############\n"
	);

	omp_start = omp_get_wtime();
	
	// Energy grid built. Now to enter parallel region
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat, seed) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs)
	{	
		double macro_xs_vector[5];
		thread = omp_get_thread_num();
		seed = thread+1;
		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( INFO && thread == 0 && i % 10000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);
			// Randomly pick an energy and material for the particle
			//p_energy = (double) rand() / (double) RAND_MAX;
			p_energy = rn(&seed);
			mat = pick_mat(&seed); 
		
			// This returns the macro_xs_vector, but we're not going
			//to do anything with it in this program, so return value
			//is written over stored.
			calculate_macro_xs( p_energy, mat, n_isotopes,
			                    n_gridpoints, num_nucs, concs,
			                    energy_grid, nuclide_grids, mats,
                                macro_xs_vector );
		}
		//free(macro_xs_vector);	
	}
	printf("\n" );
	printf("Simulation complete.\n" );

	omp_end = omp_get_wtime();
	
	printf(
	"###################################################################"
	"#############\n");
	center_print("RESULTS", 79);
	printf(
	"###################################################################"
	"#############\n"
	);

	// Print the results
	printf("Threads:   %d\n", nthreads);
	printf("Runtime:   %.3lf seconds\n", omp_end-omp_start);
	printf("Lookups:   %d\n", lookups);
	printf("Lookups/s: %.0lf\n", (double) lookups / (omp_end-omp_start));
	printf(
	"###################################################################"
	"#############\n"
	);

	// For bechmarking, output lookup/s data to file
	if( SAVE )
	{
		FILE * out = fopen( "results.txt", "a" );
		fprintf(out, "%.0lf\n", (double) lookups / (omp_end-omp_start));
		fclose(out);
	}
	
	#ifdef __PAPI
	counter_stop(&eventset);
	#endif

	return 0;
}
