#include "XSbench_header.h"

double * maketree( double * A, double **ptrs, int n )
{
	int i,j;
	int levels = (int) ceil(log2(n));
	printf("Tree is %d levels deep\n", levels);
	double * tree = (double *) malloc( sizeof(double) * (pow(2,levels)-1));
	assert(tree != NULL);
	int idx = 0;
	for( i = 0; i < levels; i++ )
	{
		int elements_in_level = (int) pow(2,i);
		for( j = 1; j <= elements_in_level*2; j+=2 )
		{
			int A_index = (int) ceil(j * n / ( 2.0 * elements_in_level ));
			tree[idx++] = A[A_index-1];
		}
	}
	return tree;
}
