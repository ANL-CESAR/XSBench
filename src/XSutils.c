#include "XSbench_header.h"

// Allocates energy grid pointer matrix
int ** pmatrix(size_t m, size_t n)
{
	int i;
	int * full = (int *)malloc(m*n*sizeof(int));
	int ** M   = (int **)malloc(m*sizeof(int*));

	for(i=0; i<m; i++){
		M[i] = full;
		full += n;
	}
	return M;
}

// Allocates nuclide matrix
double *** gpmatrix(size_t m, size_t n)
{
	int i, j;
	int l = 6;
	double * full = (double *)malloc(m*n*l*sizeof(double));
	double ** N   = (double **)malloc(m*n*sizeof(double*));
	double *** M  = (double ***)malloc(m*sizeof(double**));

	for(i=0; i<m; i++){
		M[i] = N;
		N += n;
		for(j=0; j<n; j++){
			M[i][j] = full;
			full += l;
		}
	}
	return M;
}

// Frees nuclide matrix
void gpmatrix_free(double *** M)
{
	free(**M);
	free(*M);
	free(M);
}

void pmatrix_free(double ** M)
{
	free(*M);
	free(M);
}

// Compare function for two grid points. Used for sorting during init
int NGP_compare(const void * a, const void * b)
{
	double *i = *(double **) a;
	double *j = *(double **) b;
	
	if(i[0] > j[0]) return 1;
	else if(i[0] < j[0]) return -1;
	else return 0;
}

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

size_t estimate_mem_usage(Inputs in)
{
	size_t single_nuclide_grid = in.n_gridpoints * 6 * sizeof(double);
	size_t all_nuclide_grids   = in.n_isotopes * single_nuclide_grid;
	size_t size_GridPoint      = sizeof(double) + in.n_isotopes*sizeof(int);
	size_t size_UEG            = in.n_isotopes*in.n_gridpoints * size_GridPoint;
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
	long i, j;
	// Dump Nuclide Grid Data
	for(i=0; i<n_isotopes; i++){
		//fwrite(nuclide_grids[i], 6*sizeof(double), n_gridpoints, fp);
		for(j=0; j<n_gridpoints; j++){
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 0], sizeof(double), 1, fp);
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 1], sizeof(double), 1, fp);
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 2], sizeof(double), 1, fp);
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 3], sizeof(double), 1, fp);
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 4], sizeof(double), 1, fp);
			fwrite(&nuclide_grids[i*n_gridpoints*6 + j*6 + 5], sizeof(double), 1, fp);
		}
	}

	//Dump UEG Data
	fwrite(energy_grid, sizeof(double), n_isotopes*n_gridpoints, fp);
	fwrite(grid_ptrs, sizeof(int), n_isotopes*n_gridpoints*n_isotopes, fp);

	fclose(fp);
}

void binary_read(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid, int * grid_ptrs)
{
	int stat;
	long i, j;
	FILE * fp = fopen("XS_data.dat", "rb");
	// Read Nuclide Grid Data
	for(i=0; i<n_isotopes; i++){
		//stat = fread(nuclide_grids[i], 6*sizeof(double), n_gridpoints, fp);
		for(j=0; j<n_gridpoints; j++){
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 0], sizeof(double), 1, fp);
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 1], sizeof(double), 1, fp);
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 2], sizeof(double), 1, fp);
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 3], sizeof(double), 1, fp);
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 4], sizeof(double), 1, fp);
			stat = fread(&nuclide_grids[i*n_gridpoints*6 + j*6 + 5], sizeof(double), 1, fp);
		}
	}
	//Read UEG Data
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
