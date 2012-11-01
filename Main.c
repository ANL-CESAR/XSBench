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
	int **mats;
	double **concs;
	int *num_nucs;
	load_mats( mats, concs, num_nucs );

	for( int k = 0; k < 12; k++ )
	{		
		printf("mat = %d\n", k);
		for( int l = 0; l < num_nucs[k]; l++ )
		{
			printf("nuc = %d. conc = %lf\n", l, concs[k][l]);
		}
	}

	if( INFO ) printf("Using %d threads.\n", nthreads);

	omp_start = omp_get_wtime();
	
	// variables we'll need
	

	// Energy grid built. Now to make a loop.
	#pragma omp parallel default(none) \
	private(i, thread) \
	shared( max_procs, n_isotopes, n_gridpoints, \
	energy_grid, nuclide_grids, lookups, nthreads, \
	mats, concs, num_nucs)
	{	
		thread = omp_get_thread_num();

		if( INFO ) printf("entering parallel region...\n");
		#pragma omp for
		for( i = 0; i < lookups; i++ )
		{
			if( INFO ) printf("i = %d\n", i);
			if( DEBUG && thread == 0 && i % 100 == 0 )
				printf("\rRunning Sim... Calculating XS's... (%.1lf%% completed)",
						i / ( lookups / (double) nthreads ) * 100.0);

			double p_energy = (double) rand() / (double) RAND_MAX;
		
			printf("picking mats");	
			int mat = pick_mat(); 
			
			double macro_xs = 0;
			int p_nuc;
			double conc;
			for( int j = 0; j < num_nucs[mat]; j++ )
			{
				printf("mat = %d. num_nucs[mat] = %d. mats[mat][j] = %d. "
				       "concs[mat][j] = %lf\n", mat, num_nucs[mat],
							 mats[mat][j], concs[mat][j]);
				p_nuc = mats[mat][j];
				conc = concs[mat][j];
				macro_xs += calculate_micro_xs( p_energy, p_nuc, n_isotopes,
					                  n_gridpoints, energy_grid, nuclide_grids )
				            * conc;
			}
			macro_xs = macro_xs / num_nucs[mat];
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
