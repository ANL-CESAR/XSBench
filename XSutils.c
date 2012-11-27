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
	"                   __   __ ___________                 _                        \n"
	"                   \\ \\ / //  ___| ___ \\               | |                       \n"
	"                    \\ V / \\ `--.| |_/ / ___ _ __   ___| |__                     \n"
	"                    /   \\  `--. \\ ___ \\/ _ \\ '_ \\ / __| '_ \\                    \n"
	"                   / /^\\ \\/\\__/ / |_/ /  __/ | | | (__| | | |                   \n"
	"                   \\/   \\/\\____/\\____/ \\___|_| |_|\\___|_| |_|                   \n\n"
	"###################################################################"
	"#############\n");
	center_print("Developed at Argonne National Laboratory", 79);
	center_print("Version: 0", 79);
	printf(
	"###################################################################"
	"#############\n");
}

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
// From "Numerical Recipes", Second Edition
// Press, Teukolsky, Vetterling, Flannery
float ran2( )
{
	int j;
	long k;
	float temp;

	if( *idnum <= 0 ) {
		if( -(*idnum) < 1) *idnum=1;
		else *idnum = -(*idnum);
		idnum2=(*idnum);
		for( j = NTAB + 7; j >= 0; j-- )
		{
			k = (*idnum)/IQ1;
			*idnum = IA1 * (*idnum - k * IQ1 ) - k * IR1;
			if( *idnum < 0 )
				*idnum += IM1;
			if ( j < NTAB )
				iv[j] = *idnum;
		}
		iy=iv[0];
	}
	k=(*idnum)/IQ1;
	*idnum=IA1 * (*idnum-k * IQ1) - k * IR1;
	if( *idnum < 0 )
		*idnum += IM1;
	k = idnum2/IQ2;
	idnum2 = IA2 * (idnum2-k*IQ2) - k * IR2;
	if( idnum2 < 0 )
		idnum2 += IM2;
	j = iy / NDIV;
	iy = iv[j] - idnum2;
	iv[j] = *idnum;
	if( iy < 1 ) iy += IMM1;
	if( (temp=AM*iy) > RNMX)
		return RNMX;
	else
		return temp;
}
*/
