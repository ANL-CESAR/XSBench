#include "XSbench_header.h"

int main( int argc, char* argv[] )
{
	srand(time(NULL));
	int n_isotopes = 300;
	int n_gridpoints = 100;
	int lookups = 100000;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	int max_procs = omp_get_num_procs();
	omp_set_num_threads( max_procs );

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

	double omp_start = omp_get_wtime();
	
	// Energy grid built. Now to make a loop.
	int i;
	#pragma omp parallel for default(none) \
	private( i ) \
	shared( n_isotopes, n_gridpoints, energy_grid, nuclide_grids, lookups )
	for( i = 0; i < lookups; i++ )
	{
		if( DEBUG ) printf("\rRunning MC... (%d of %d completed)",i+1,lookups );

		double p_energy = (double) rand() / (double) RAND_MAX;
		int p_nuc = rand() % n_isotopes;

		calculate_micro_xs( p_energy, p_nuc, n_isotopes,
                        n_gridpoints, energy_grid, nuclide_grids );
	}	
	if( DEBUG ) printf("\n" );

	double omp_end = omp_get_wtime();

	if( INFO ) printf("Runtime:   %.3lf seconds\n", omp_end-omp_start);
	if( INFO ) printf("Lookups:   %d\n", lookups);
	if( INFO ) printf("Lookups/s: %lf\n", (double) lookups / (omp_end-omp_start));

	return 0;
}
