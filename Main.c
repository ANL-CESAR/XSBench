#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	srand(time(NULL));
	int n_isotopes = 68;
	int n_gridpoints = 200;
	int lookups = 1000000;
	int max_procs = omp_get_num_procs();
	int i, thread, nthreads;
	double omp_start, omp_end;

	if( argc == 2 )
		nthreads = atoi(argv[1]);
	else
		nthreads = max_procs;
	
	omp_set_num_threads(nthreads); 

	logo();
	
	// Allocate & fill energy grids
	if( DEBUG ) printf("Generating Nuclide Energy Grids...\n");
	
	NuclideGridPoint ** nuclide_grids = gpmatrix( n_isotopes, n_gridpoints );
	generate_grids( nuclide_grids, n_isotopes, n_gridpoints );	
	
	// Sort grids
	sort_nuclide_grids( nuclide_grids, n_isotopes );

	// Build Unionized Grid Framework
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

	omp_start = omp_get_wtime();
	
	// variables we'll need
	double p_energy;
	double macro_xs;
	int mat;

	// Energy grid built. Now to make a loop.
	#pragma omp parallel default(none) \
	private(i, thread, p_energy, macro_xs, mat) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs)
	{	
		thread = omp_get_thread_num();

		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			if( DEBUG && thread == 0 && i % 100 == 0 )
				printf("\rRunning Sim... Calculating XS's... (%.1lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);

			p_energy = (double) rand() / (double) RAND_MAX;
		
			mat = pick_mat(); 
		
			macro_xs = calculate_macro_xs( p_energy, mat, n_isotopes,
			                               n_gridpoints, num_nucs, concs,
			                               energy_grid, nuclide_grids, mats );
		}	
	}
	if( DEBUG ) printf("\n" );

	omp_end = omp_get_wtime();

	if( INFO ) printf("Runtime:   %.3lf seconds\n", omp_end-omp_start);
	if( INFO ) printf("Lookups:   %d\n", lookups);
	if( INFO ) printf("Lookups/s: %.0lf\n",
		                (double) lookups / (omp_end-omp_start));

	return 0;
}
