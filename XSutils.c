#include "XSbench_header.h"

/*
#ifdef __bgq__
#include </bgsys/drivers/ppcfloor/spi/include/kernel/location.h>
#endif
*/

// Allocates nuclide matrix
NuclideGridPoint ** gpmatrix(size_t m, size_t n)
{
	int i,j;
	NuclideGridPoint * full = (NuclideGridPoint *) malloc( m * n *
	                          sizeof( NuclideGridPoint ) );
	NuclideGridPoint ** M = (NuclideGridPoint **) malloc( m *
	                          sizeof(NuclideGridPoint *) );

	for( i = 0, j=0; i < m*n; i++ )
		if( i % n == 0 )
			M[j++] = &full[i];

	return M;
}

// Frees nuclide matrix
void gpmatrix_free( NuclideGridPoint ** M )
{
	free( *M );
	free( M );
}

// Compare function for two grid points. Used for sorting during init
int NGP_compare( const void * a, const void * b )
{
	NuclideGridPoint *i, *j;

	i = (NuclideGridPoint *) a;
	j = (NuclideGridPoint *) b;

	if( i->energy > j->energy )
		return 1;
	else if ( i->energy < j->energy)
		return -1;
	else
		return 0;
}

// Prints program logo
void logo(void)
{
	printf(
	"###################################################################"
	"#############\n"
	"                   __   __ ___________                 _                        \n"
	"                   \\ \\ / //  ___| ___ \\               | |                       \n"
	"                    \\ V / \\ `--.| |_/ / ___ _ __   ___| |__                     \n"
	"                    /   \\  `--. \\ ___ \\/ _ \\ '_ \\ / __| '_ \\                    \n"
	"                   / /^\\ \\/\\__/ / |_/ /  __/ | | | (__| | | |                   \n"
	"                   \\/   \\/\\____/\\____/ \\___|_| |_|\\___|_| |_|                   \n\n"
	"###################################################################"
	"#############\n");
	center_print("Developed at Argonne National Laboratory", 79);
	center_print("Version: 2", 79);
	printf(
	"###################################################################"
	"#############\n");
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}

// Binary Search function for nuclide grid
// Returns ptr to energy less than the quarry that is closest to the quarry
NuclideGridPoint * binary_search( NuclideGridPoint * A, double quarry, int n )
{
	int min = 0;
	int max = n-1;
	int mid;
	
	// checks to ensure we're not reading off the end of the grid
	if( A[0].energy > quarry )
		return &A[0];
	else if( A[n-1].energy < quarry )
		return &A[n-2];
	
	// Begins binary search	
	while( max >= min )
	{
		mid = min + floor( (max-min) / 2.0);
		if( A[mid].energy < quarry )
			min = mid+1;
		else if( A[mid].energy > quarry )
			max = mid-1;
		else
			return &A[mid];
	}
	return &A[max];
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

/*
int get_bgq_core(void)
{
#ifdef __bgq__
	int core = Kernel_ProcessorCoreID();
	return core;
#else
	return -1;
#endif
}*/
