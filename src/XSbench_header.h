#ifndef __XSBENCH_HEADER_H__
#define __XSBENCH_HEADER_H__

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<omp.h>
#include<unistd.h>
#include<sys/time.h>

// I/O Specifiers
#define INFO 1
#define DEBUG 1
#define SAVE 1

#ifdef __cplusplus
#define restrict __restrict__
#endif

typedef struct{
	int nthreads;
	long n_isotopes;
	long n_gridpoints;
	int lookups;
	char * HM;
} Inputs;

// io.c function prototypes
void logo(int version);
void center_print(const char *s, int width);
void print_results(Inputs in, int mype, double runtime, int nprocs, unsigned long long vhash);
void print_inputs(Inputs in, int nprocs, int version);
void border_print(void);
void fancy_int(long a);
void print_CLI_error(void);
Inputs read_CLI(int argc, char * argv[]);

// XSutils.c function prototypes
int NGP_compare(const void * a, const void * b);
int d_compare(const void * a, const void * b);
double *** gpmatrix(size_t m, size_t n);
int ** pmatrix(size_t m, size_t n);
void gpmatrix_free(double *** M);
void pmatrix_free(double ** M);
int binary_search(double * A, double quarry, int n);
double rn(unsigned long * seed);
double rn_v(void);
unsigned int hash(unsigned char *str, int nbins);
size_t estimate_mem_usage(Inputs in);
void binary_dump(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid, int * grid_ptrs);
void binary_read(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid, int * grid_ptrs);
double timer();

// GridInit.c function prototypes
void generate_grids(double * nuclide_grids, long n_isotopes, long n_gridpoints);
void generate_grids_v(double * nuclide_grids, long n_isotopes, long n_gridpoints);
void sort_nuclide_grids(double * nuclide_grids, long n_isotopes, long n_gridpoints); 
double * generate_energy_grid(long n_isotopes, long n_gridpoints, double * nuclide_grids);
int * generate_grid_ptrs(long n_isotopes, long n_gridpoints, double * nuclide_grids, double * energy_grid);

// CalculateXS.c function prototypes
void calculate_macro_xs(double p_energy, int mat, long n_isotopes, long n_gridpoints,
			int * restrict num_nucs, double * restrict concs,
			double * restrict energy_grid, double * restrict nuclide_grids,
			int * restrict grid_ptrs, int * restrict mats, int * restrict mats_ptr,
			double * restrict macro_xs_vector);
void calculate_micro_xs(double p_energy, int nuc, long n_isotopes, long n_gridpoints,
			double * restrict energy_grid, double * restrict nuclide_grids,
			int * restrict grid_ptrs, int idx, double * restrict xs_vector);
long grid_search(long n, double quarry, double * A);

// Materials.c function prototypes
int * load_num_nucs(long n_isotopes);
int * load_mats_ptr(int * num_nucs);
int * load_mats(int * num_nucs, int * mats_ptr, int size_mats, long n_isotopes);
double * load_concs(int size_mats);
double * load_concs_v(int size_mats);
int pick_mat(unsigned long * seed);

#endif
