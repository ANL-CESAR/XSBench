#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include<unistd.h>
#include<sys/time.h>

#define INFO 1
#define DEBUG 1

// For the numerical methods rand() algorithm
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


typedef struct{
	double energy;
	double micro_xs;
} NuclideGridPoint;

typedef struct{
	double energy;
	NuclideGridPoint ** xs_ptrs;
} GridPoint;

void logo(void);

NuclideGridPoint ** gpmatrix(size_t m, size_t n);

void gpmatrix_free( NuclideGridPoint ** M );

int NGP_compare( const void * a, const void * b );

void generate_grids( NuclideGridPoint ** nuclide_grids,
                     int n_isotopes, int n_gridpoints );

void sort_nuclide_grids( NuclideGridPoint ** nuclide_grids, int n_isotopes);

GridPoint * generate_energy_grid( int n_isotopes, int n_gridpoints,
                                  NuclideGridPoint ** nuclide_grids);

void set_grid_ptrs( GridPoint * energy_grid, NuclideGridPoint ** nuclide_grids,
                    int n_isotopes, int n_gridpoints );

NuclideGridPoint * binary_search( NuclideGridPoint * A, double quarry, int n );

double calculate_macro_xs( double p_energy, int mat, int n_isotopes,
                           int n_gridpoints, int * num_nucs,
                           double ** concs, GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids, int ** mats );

double calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids, int idx );

int grid_search( int n, double quarry, GridPoint * A);

int * load_num_nucs(void);
int ** load_mats( int * num_nucs );
double ** load_concs( int * num_nucs );
int pick_mat(unsigned long * seed);
double rn(unsigned long * seed);
float ran2( void );
