#include "XSbench_header.h"

TreeDataPtrs maketree( double * A, int **ptrs, long n )
{
	long levels = (long) ceil(log2(n));
	long n_tree_elements = (long) pow(2,levels) - 1;

	double * tree = (double *) malloc( sizeof(double) * n_tree_elements);
	assert(tree  != NULL);
	int **  ptree = (int ** )  malloc( sizeof(int * ) * n_tree_elements);
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

	TreeDataPtrs T;
	T.tree = tree;
	T.ptree = ptree;
	T.n_tree_elements = n_tree_elements;

	return T;
}
