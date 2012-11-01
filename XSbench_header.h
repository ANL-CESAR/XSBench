#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<omp.h>

#define INFO 1
#define DEBUG 1

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

double calculate_micro_xs( int p_energy, int nuc, int n_isotopes,
                           int n_gridpoints,
                           GridPoint * energy_grid,
                           NuclideGridPoint ** nuclide_grids );

int * load_num_nucs(void);
int ** load_mats( int * num_nucs );
double ** load_concs( int * num_nucs );
int pick_mat(void);
