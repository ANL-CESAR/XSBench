#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

// Generates randomized energy grid for each nuclide
// Note that this is done as part of initialization (serial), so
// rand() is used.
void generate_grids( Nuclide * nuclides, int n_isotopes, int n_gridpoints ) {
	for( int i = 0; i < n_isotopes; i++ )
	{
		nuclides[i].energy = (double *) malloc( sizeof(double) * n_gridpoints );
		nuclides[i].XS = (double *) malloc( sizeof(NuclideGridPoint) * n_gridpoints ); 
		for( int j = 0; j < n_gridpoints; j++ )
		{
			nuclides[i].energy[j]          =((double)rand()/(double)RAND_MAX);
			nuclides[i].XS[j].total_xs     =((double)rand()/(double)RAND_MAX);
			nuclides[i].XS[j].elastic_xs   =((double)rand()/(double)RAND_MAX);
			nuclides[i].XS[j].absorbtion_xs=((double)rand()/(double)RAND_MAX);
			nuclides[i].XS[j].fission_xs   =((double)rand()/(double)RAND_MAX);
			nuclides[i].XS[j].nu_fission_xs=((double)rand()/(double)RAND_MAX);
		}
	}
}

// Verification version of this function (tighter control over RNG)
void generate_grids_v( Nuclide * nuclides,
                     int n_isotopes, int n_gridpoints ) {
	for( int i = 0; i < n_isotopes; i++ )
	{
		nuclides[i].energy = (double *) malloc( sizeof(double) * n_gridpoints );
		nuclides[i].XS = (double *) malloc( sizeof(NuclideGridPoint) * n_gridpoints ); 
		for( int j = 0; j < n_gridpoints; j++ )
		{
			nuclides[i].energy[j]          = rn_v();
			nuclides[i].XS[j].total_xs     = rn_v();
			nuclides[i].XS[j].elastic_xs   = rn_v();
			nuclides[i].XS[j].absorbtion_xs= rn_v();
			nuclides[i].XS[j].fission_xs   = rn_v();
			nuclides[i].XS[j].nu_fission_xs= rn_v();
		}
	}
}

// Sorts the nuclide grids by energy (lowest -> highest)
// This requires a manual implementation of qsort, as Qsort can't sort by
// things it doesn't have access to
void sort_nuclide_grids( Nuclide * nuclides, int n_isotopes, int n_gridpoints )
{
	int (*cmp) (const void *, const void *);
	cmp = NGP_compare;
	
	for( int i = 0; i < n_isotopes; i++ )
	{
		Nuclide * A = nuclides[i];
		nuclide_qsort( A.energy, A.XS, 0, n_gridpoints-1 );
	}
	
	// error debug check
	/*
	for( int i = 0; i < n_isotopes; i++ )
	{
		printf("NUCLIDE %d ==============================\n", i);
		for( int j = 0; j < n_gridpoints; j++ )
			printf("E%d = %lf\n", j, nuclide_grids[i][j].energy);
	}
	*/
}

void nuclide_qsort(double * energy, NuclideGridPoint * XS, int low,int high)
{
	int pivot,j,i;
	double dtemp;
	NuclideGridPoint XStemp;

	if(low<high)
	{
		pivot = low;
		i = low;
		j = high;

		while(i<j)
		{
			while((energy[i]<=energy[pivot])&&(i<high))
				i++;

			while(energy[j]>energy[pivot])
				j--;

			if(i<j)
			{ 
				dtemp=energy[i];
				XStemp=XS[i];

				energy[i]=energy[j];
				XS[i]=XS[j];

				energy[j]=dtemp;
				XS[j] = XStemp;
			}
		}

		dtemp=energy[pivot];
		XStemp=XS[pivot]

		energy[pivot]=energy[j];
		XS[pivot] = XS[j];

		energy[j]=dtemp;
		XS[j]=XStemp;

		nuclide_qsort(energy,XS,low,j-1);
		nuclide_qsort(energy,XS,j+1,high);
	}
}

// Allocates unionized energy grid, and assigns union of energy levels
// from nuclide grids to it.
Grid generate_energy_grid( int n_isotopes, int n_gridpoints,
		Nuclide * nuclides) {
	int mype = 0;

#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#endif

	if( mype == 0 ) printf("Generating Unionized Energy Grid...\n");

	Grid grid;
	int pts = n_isotopes * n_gridpoints;
	grid.energy = (double *) malloc( sizeof(double *) * pts );
	grid.nuclides = (int *) malloc( sizeof(GridPoint) * pts );

	if( mype == 0 ) printf("Copying and Sorting all nuclide grids...\n");
	
	for( int i = 0; i < n_isotopes; i++ )
		memcpy( &grid.energy[i*n_gridpoints], nuclides[i].energy, n_gridpoints*sizeof(double));
	
	qsort( grid.energy, sizeof(double), pts, dbl_cmp );

	return grid;
}

dbl_cmp(const void * a, const void * b)
{
	double * A = (double *) a;
	double * B = (double *) b;
	if( A > B )
		return 1;
	else if ( A < B )
		return -1;
	else
		return 0;
}


//JRT JRT JRT JRT -- THIS IS WHERE YOU NEED TO PICKUP WORK!!!!

// Searches each nuclide grid for the closest energy level and assigns
// pointer from unionized grid to the correct spot in the nuclide grid.
// This process is time consuming, as the number of binary searches
// required is:  binary searches = n_gridpoints * n_isotopes^2
void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
		int n_isotopes, int n_gridpoints )
{
	int mype = 0;

#ifdef MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#endif

	if( mype == 0 ) printf("Assigning pointers to Unionized Energy Grid...\n");
#pragma omp parallel for default(none) \
	shared( energy_grid, nuclide_grids, n_isotopes, n_gridpoints, mype )\
	num_threads(32)
	for( int i = 0; i < n_isotopes * n_gridpoints ; i++ )
	{
		double quarry = energy_grid[i].energy;
		if( INFO && mype == 0 && omp_get_thread_num() == 0 && i % 200 == 0 )
			printf("\rAligning Unionized Grid...(%.0lf%% complete)",
					100.0 * (double) i / (n_isotopes*n_gridpoints /
						omp_get_num_threads())     );
		for( int j = 0; j < n_isotopes; j++ )
		{
			// j is the nuclide i.d.
			// log n binary search
			energy_grid[i].xs_ptrs[j] = 
				binary_search( nuclide_grids[j], quarry, n_gridpoints);
		}
	}
	if( mype == 0 ) printf("\n");

	//test
	/*
	   for( int i=0; i < n_isotopes * n_gridpoints; i++ )
	   for( int j = 0; j < n_isotopes; j++ )
	   printf("E = %.4lf\tNuclide %d->%p->%.4lf\n",
	   energy_grid[i].energy,
	   j,
	   energy_grid[i].xs_ptrs[j],
	   (energy_grid[i].xs_ptrs[j])->energy
	   );
	 */
}
