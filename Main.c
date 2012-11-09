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
	printf(
	"\tINPUT SUMMARY\n"
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
	"#############\n"
	"\tINITIALIZATION\n"
	"###################################################################"
	"#############\n"
	);

	// Allocate & fill energy grids
	if( DEBUG ) printf("Generating Nuclide Energy Grids...\n");
	
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
	if( INFO ) printf("Loading Mats...\n");
	int *num_nucs = load_num_nucs();
	int **mats = load_mats(num_nucs);
	double **concs = load_concs(num_nucs);

	if( INFO ) printf("Using %d threads.\n", nthreads);
	
	printf(
	"###################################################################"
	"#############\n"
	"\tSIMULATION\n"
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
		double * macro_xs_vector = (double *) malloc( 5 * sizeof(double));
		thread = omp_get_thread_num();
		seed = thread+1;
		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			// Status text
			if( DEBUG && thread == 0 && i % 10000 == 0 )
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
		free(macro_xs_vector);	
	}
	if( DEBUG ) printf("\n" );
	printf("Simulation complete.\n" );

	omp_end = omp_get_wtime();
	
	printf(
	"###################################################################"
	"#############\n"
	"\tRESULTS\n"
	"###################################################################"
	"#############\n"
	);

	// Print the results
	if( INFO ) printf("Runtime:   %.3lf seconds\n", omp_end-omp_start);
	if( INFO ) printf("Lookups:   %d\n", lookups);
	if( INFO ) printf("Lookups/s: %.0lf\n",
		               (double) lookups / (omp_end-omp_start));
	printf(
	"###################################################################"
	"#############\n"
	);

	return 0;
}
