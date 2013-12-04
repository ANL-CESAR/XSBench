#include "XSbench_header.h"


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


// Binary Search function for nuclide grid
// Returns ptr to energy less than the quarry that is closest to the quarry
int binary_search( NuclideGridPoint * A, double quarry, int n )
{
	int min = 0;
	int max = n-1;
	int mid;
	
	// checks to ensure we're not reading off the end of the grid
	if( A[0].energy > quarry )
		return 0;
	else if( A[n-1].energy < quarry )
		return n-2;
	
	// Begins binary search	
	while( max >= min )
	{
		mid = min + floor( (max-min) / 2.0);
		if( A[mid].energy < quarry )
			min = mid+1;
		else if( A[mid].energy > quarry )
			max = mid-1;
		else
			return mid;
	}
	return max;
}

// (fixed) binary search for energy on unionized energy grid
// returns lower index
int grid_search( int n, double quarry, GridPoint * A)
{
	int lowerLimit = 0;
	int upperLimit = n-1;
	int examinationPoint;
	int length = upperLimit - lowerLimit;

	while( length > 1 )
	{
		examinationPoint = lowerLimit + ( length / 2 );
		
		if( A[examinationPoint].energy > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
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

	while (c = *str++)
		hash = ((hash << 5) + hash) + c;

	return hash % nbins;
}
