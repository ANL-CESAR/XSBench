#include "XSbench_header.h"

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

void gpmatrix_free( NuclideGridPoint ** M )
{
	free( *M );
	free( M );
}


int NGP_compare( const void * a, const void * b )
{
	NuclideGridPoint i, j;
	i = *( (NuclideGridPoint *) a);
	j = *(( NuclideGridPoint *) b);

	if( i.energy > j.energy )
		return 1;
	else if ( i.energy < j.energy)
		return -1;
	else
		return 0;
}

void logo(void)
{
	printf(
	"###################################################################"
	"#############\n"
	"__   __ ___________                 _     \n"
	"\\ \\ / //  ___| ___ \\               | |    \n"
	 " \\ V / \\ `--.| |_/ / ___ _ __   ___| |__  \n"
	  " /   \\  `--. \\ ___ \\/ _ \\ '_ \\ / __| '_ \\ \n"
		"/ /^\\ \\/\\__/ / |_/ /  __/ | | | (__| | | |\n"
		"\\/   \\/\\____/\\____/ \\___|_| |_|\\___|_| |_|\n\n"
	"###################################################################"
	"#############\n"
	"Developed at Argonne National Laboratory\n"
	"Version: 0\n"
	"###################################################################"
	"#############\n");
}

NuclideGridPoint * binary_search( NuclideGridPoint * A, double quarry, int n )
{
	int min = 0;
	int max = n-1;
	int mid;
	
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

double rn(void)
{
	double intpart;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return modf( ((double) tv.tv_usec)/791.9, &intpart);
}
