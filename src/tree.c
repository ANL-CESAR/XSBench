#include "XSbench_header.h"

TreeStuff maketree( double * A, int **ptrs, long n )
{
	long levels = (int) ceil(log2(n));

	double * tree = (double *) malloc( sizeof(double) * (pow(2,levels)-1));
	assert(tree  != NULL);
	int **  ptree = (int ** )  malloc( sizeof(int * ) * (pow(2,levels)-1));
	assert(ptree != NULL);

	long idx = 0;
	for( long i = 0; i < levels; i++ )
	{
		long elements_in_level = (int) pow(2,i);
		for( long j = 1; j <= elements_in_level*2; j+=2 )
		{
			long A_index = (long) ceil(j * n / ( 2.0 * elements_in_level ));
			tree[idx] = A[A_index-1];
			ptree[idx++] = ptrs[A_index-1];
		}
	}

	TreeStuff T;
	T.tree = tree;
	T.ptree = ptree;

	return T;
}
