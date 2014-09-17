#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

// Generates randomized energy grid for each nuclide
// Note that this is done as part of initialization (serial), so rand() is used.
void generate_grids(double *** nuclide_grids, long n_isotopes, long n_gridpoints)
{
	long i, j;
	for(i=0; i<n_isotopes; i++){
		for(j=0; j<n_gridpoints; j++){
			nuclide_grids[i][j][0] = ((double)rand()/(double)RAND_MAX); //energy
			nuclide_grids[i][j][1] = ((double)rand()/(double)RAND_MAX); //total xs
			nuclide_grids[i][j][2] = ((double)rand()/(double)RAND_MAX); //elastic xs
			nuclide_grids[i][j][3] = ((double)rand()/(double)RAND_MAX); //absorption xs
			nuclide_grids[i][j][4] = ((double)rand()/(double)RAND_MAX); //fission xs
			nuclide_grids[i][j][5] = ((double)rand()/(double)RAND_MAX); //nu fission xs
		}
	}
}

// Verification version of this function (tighter control over RNG)
void generate_grids_v(double *** nuclide_grids, long n_isotopes, long n_gridpoints)
{
	long i, j;
	for(i=0; i<n_isotopes; i++){
		for(j=0; j<n_gridpoints; j++){
			nuclide_grids[i][j][0] = rn_v(); //energy
			nuclide_grids[i][j][1] = rn_v(); //total xs
			nuclide_grids[i][j][2] = rn_v(); //elastic xs
			nuclide_grids[i][j][3] = rn_v(); //absorption xs
			nuclide_grids[i][j][4] = rn_v(); //fission xs
			nuclide_grids[i][j][5] = rn_v(); //nu fission xs
		}
	}
}

// Sorts the nuclide grids by energy (lowest -> highest)
void sort_nuclide_grids(double *** nuclide_grids, long n_isotopes, long n_gridpoints)
{
	long i, j;
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;

	for(i=0; i<n_isotopes; i++){
		qsort(&nuclide_grids[i][0], n_gridpoints, sizeof(double*), cmp);
	}
}

// Allocates unionized energy grid, and assigns union of energy levels
// from nuclide grids to it.
double * generate_energy_grid(long n_isotopes, long n_gridpoints, double *** nuclide_grids)
{
	long i, j, n_unionized_grid_points;
	int mype = 0;
	double * energy_grid;
	int (*cmp) (const void *, const void *);
	cmp = d_compare;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if(mype == 0) printf("Generating Unionized Energy Grid...\n");

	n_unionized_grid_points = n_isotopes*n_gridpoints;
	energy_grid = (double *)malloc(n_unionized_grid_points*sizeof(double));

	if(mype == 0) printf("Copying and sorting nuclide grid energies...\n");

	for(i=0; i<n_isotopes; i++){
		for(j=0; j<n_gridpoints; j++){
			energy_grid[i*n_gridpoints + j] = nuclide_grids[i][j][0];
		}
	}
	qsort(energy_grid, n_unionized_grid_points, sizeof(double), cmp);

	return energy_grid;
}

// Searches each nuclide grid for the closest energy level and assigns
// pointer from unionized grid to the correct spot in the nuclide grid.
// This process is time consuming, as the number of binary searches
// required is:  binary searches = n_gridpoints * n_isotopes^2
int * generate_grid_ptrs(long n_isotopes, long n_gridpoints, double *** nuclide_grids, double * energy_grid)
{
	long i, j;
	double quarry;
	int mype = 0;
	int * grid_ptrs;

	#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	grid_ptrs = (int *) malloc(n_isotopes*n_gridpoints*n_isotopes*sizeof(int));

	if(mype == 0) printf("Assigning pointers to Unionized Energy Grid...\n");

	#pragma omp parallel for default(none) \
	shared( energy_grid, grid_ptrs, nuclide_grids, n_isotopes, n_gridpoints, mype ) \
	private( quarry, i, j )
	for(i=0; i<n_isotopes*n_gridpoints; i++){
		quarry = energy_grid[i];
		if(INFO && mype == 0 && omp_get_thread_num() == 0 && i % 200 == 0)
			printf("\rAligning Unionized Grid...(%.0lf%% complete)",
			       100.0 * (double) i / (n_isotopes*n_gridpoints /
				                         omp_get_num_threads())     );
		for(j=0; j<n_isotopes; j++){
			// j is the nuclide i.d.
			// log n binary search
			grid_ptrs[n_isotopes*i + j] = binary_search(nuclide_grids[j], quarry, n_gridpoints);
		}
	}
	if(mype == 0) printf("\n");

	return grid_ptrs;
}
