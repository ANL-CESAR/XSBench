#include "XSbench_header.h"

int d_compare(const void * a, const void * b)
{
	double i = *(double *) a;
	double j = *(double *) b;

	if(i > j) return 1;
	else if(i < j) return -1;
	return 0;
}

// Binary Search function for nuclide grid
// Returns ptr to energy less than the quarry that is closest to the quarry
int binary_search(double * A, double quarry, int n)
{
	int min = 0;
	int max = n-1;
	int mid;
	
	// checks to ensure we're not reading off the end of the grid
	if(A[0*6 + 0] > quarry) return 0;
	else if(A[(n-1)*6 + 0] < quarry) return n-2;
	
	// Begins binary search	
	while(max >= min)
	{
		mid = min + floor((max-min) / 2.0);
		if(A[mid*6 + 0] < quarry) min = mid+1;
		else if(A[mid*6 + 0] > quarry) max = mid-1;
		else return mid;
	}
	return max;
}

// Park & Miller Multiplicative Conguential Algorithm
// From "Numerical Recipes" Second Edition
double rn(unsigned long * seed)
{
	double ret;
	unsigned long n1;
	unsigned long a = 16807;
	unsigned long m = 2147483647;
	n1 = ( a * (*seed) ) % m;
	*seed = n1;
	ret = (double) n1 / m;
	return ret;
}

// RNG Used for Verification Option.
// This one has a static seed (must be set manually in source).
// Park & Miller Multiplicative Conguential Algorithm
// From "Numerical Recipes" Second Edition
double rn_v(void)
{
	static unsigned long seed = 1337;
	double ret;
	unsigned long n1;
	unsigned long a = 16807;
	unsigned long m = 2147483647;
	n1 = ( a * (seed) ) % m;
	seed = n1;
	ret = (double) n1 / m;
	return ret;
}

unsigned int hash(unsigned char *str, int nbins)
{
	unsigned int hash = 5381;
	int c;
	while (c = *str++) hash = ((hash << 5) + hash) + c;
	return hash % nbins;
}

size_t estimate_mem_usage(long n_isotopes, long n_gridpoints)
{
	size_t single_nuclide_grid = n_gridpoints * 6 * sizeof(double);
	size_t all_nuclide_grids   = n_isotopes * single_nuclide_grid;
	size_t size_GridPoint      = sizeof(double) + n_isotopes*sizeof(int);
	size_t size_UEG            = n_isotopes*n_gridpoints * size_GridPoint;
	size_t memtotal;

	memtotal          = all_nuclide_grids + size_UEG;
	all_nuclide_grids = all_nuclide_grids / 1048576;
	size_UEG          = size_UEG / 1048576;
	memtotal          = memtotal / 1048576;
	return memtotal;
}

void binary_dump(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid, int * grid_ptrs)
{
	FILE * fp = fopen("XS_data.dat", "wb");

	fwrite(nuclide_grids, sizeof(double), n_isotopes*n_gridpoints*6, fp);
	fwrite(energy_grid, sizeof(double), n_isotopes*n_gridpoints, fp);
	fwrite(grid_ptrs, sizeof(int), n_isotopes*n_gridpoints*n_isotopes, fp);

	fclose(fp);
}

void binary_read(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid, int * grid_ptrs)
{
	int stat;
	FILE * fp = fopen("XS_data.dat", "rb");

	stat = fread(nuclide_grids, sizeof(double), n_isotopes*n_gridpoints*6, fp);
	stat = fread(energy_grid, sizeof(double), n_isotopes*n_gridpoints, fp);
	stat = fread(grid_ptrs, sizeof(int), n_isotopes*n_gridpoints*n_isotopes, fp);

	fclose(fp);
}

double timer()
{
	struct timeval time;
	gettimeofday(&time, 0);
	return time.tv_sec + time.tv_usec / 1000000.0;
}
