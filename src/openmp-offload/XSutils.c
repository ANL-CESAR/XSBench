#include "XSbench_header.h"

int double_compare(const void * a, const void * b)
{
	double A = *((double *) a);
	double B = *((double *) b);

	if( A > B )
		return 1;
	else if( A < B )
		return -1;
	else
		return 0;
}

int NGP_compare(const void * a, const void * b)
{
	NuclideGridPoint A = *((NuclideGridPoint *) a);
	NuclideGridPoint B = *((NuclideGridPoint *) b);

	if( A.energy > B.energy )
		return 1;
	else if( A.energy < B.energy )
		return -1;
	else
		return 0;
}


size_t estimate_mem_usage( Inputs in )
{
	size_t single_nuclide_grid = in.n_gridpoints * sizeof( NuclideGridPoint );
	size_t all_nuclide_grids   = in.n_isotopes * single_nuclide_grid;
	size_t size_UEG            = in.n_isotopes*in.n_gridpoints*sizeof(double) + in.n_isotopes*in.n_gridpoints*in.n_isotopes*sizeof(int);
	size_t size_hash_grid      = in.hash_bins * in.n_isotopes * sizeof(int);
	size_t memtotal;

	if( in.grid_type == UNIONIZED )
		memtotal          = all_nuclide_grids + size_UEG;
	else if( in.grid_type == NUCLIDE )
		memtotal          = all_nuclide_grids;
	else
		memtotal          = all_nuclide_grids + size_hash_grid;

	memtotal          = ceil(memtotal / (1024.0*1024.0));
	return memtotal;
}

