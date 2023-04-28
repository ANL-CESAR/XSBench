#include "XSbench_header.h"

node_t* insert_node(node_t* root, FP_PRECISION data) {
    if (root == NULL) {
        node_t* new_node = (node_t*)malloc(sizeof(node_t));
        new_node->data = data;
        new_node->left = NULL;
        new_node->right = NULL;
        return new_node;
    } else if (data < root->data) {
        root->left = insert_node(root->left, data);
    } else if (data > root->data) {
        root->right = insert_node(root->right, data);
    }
    return root;
}

int find_node(node_t* root, FP_PRECISION data) {
    if (root == NULL) {
        return 0;
    } else if (data < root->data) {
        return find_node(root->left, data);
    } else if (data > root->data) {
        return find_node(root->right, data);
    } else {
        return 1;
    }
}

void free_tree(node_t* root) {
    if (root != NULL) {
        free_tree(root->left);
        free_tree(root->right);
        free(root);
    }
}

int FP_PRECISION_compare(const void * a, const void * b)
{
	FP_PRECISION A = *((FP_PRECISION *) a);
	FP_PRECISION B = *((FP_PRECISION *) b);

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

double get_time(void)
{
	#ifdef OPENMP
	return omp_get_wtime();
	#endif

	struct timeval timecheck;

	gettimeofday(&timecheck, NULL);
	long ms = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	double time = (double) ms / 1000.0;

	return time;
}
