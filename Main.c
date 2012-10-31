#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	srand(time(NULL));
	int n_isotopes = 300;
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
	
	if( INFO ) printf("Using %d threads.\n", nthreads);

	omp_start = omp_get_wtime();
	
	// Energy grid built. Now to make a loop.
	#pragma omp parallel default(none) \
	private(i, thread) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads )
	{	
		thread = omp_get_thread_num();

		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			if( DEBUG && thread == 0 && i % 100 == 0 )
				printf("\rRunning Sim... Calculating XS's... (%.1lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);

			double p_energy = (double) rand() / (double) RAND_MAX;
			int p_nuc = rand() % n_isotopes;

			calculate_micro_xs( p_energy, p_nuc, n_isotopes,
					n_gridpoints, energy_grid, nuclide_grids );
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
