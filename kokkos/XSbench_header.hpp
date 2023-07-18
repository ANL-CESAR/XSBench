// -*- c-basic-offset: 8; tab-width: 8; indent-tabs-mode: t; -*-
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
#include<assert.h>
#include<stdint.h>

#include <Kokkos_Core.hpp>

// Papi Header
#ifdef PAPI
#include "papi.h"
#endif

// Grid types
#define UNIONIZED 0
#define NUCLIDE 1
#define HASH 2

// Simulation types
#define HISTORY_BASED 1
#define EVENT_BASED 2

// Binary Mode Type
#define NONE 0
#define READ 1
#define WRITE 2

// Starting Seed
#define STARTING_SEED 1070

// Structures
typedef struct{
	double energy;
	double total_xs;
	double elastic_xs;
	double absorbtion_xs;
	double fission_xs;
	double nu_fission_xs;
} NuclideGridPoint;

typedef struct{
	int nthreads;
	long n_isotopes;
	long n_gridpoints;
	int lookups;
	char * HM;
	int grid_type; // 0: Unionized Grid (default)    1: Nuclide Grid
	int hash_bins;
	int particles;
	int simulation_method;
	int binary_mode;
	int kernel_id;
} Inputs;

typedef Kokkos::View<int*> IntView;
typedef Kokkos::View<double*> DoubleView;
typedef Kokkos::View<NuclideGridPoint*> PointView;
typedef Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace,
		     Kokkos::MemoryTraits<Kokkos::Unmanaged>> UIntView;
typedef Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace,
		     Kokkos::MemoryTraits<Kokkos::Unmanaged>> UDoubleView;
typedef Kokkos::View<NuclideGridPoint*, Kokkos::LayoutLeft, Kokkos::HostSpace,
		     Kokkos::MemoryTraits<Kokkos::Unmanaged>> UPointView;

typedef struct{
	IntView* d_num_nucs;				// Length = length_num_nucs;
	DoubleView* d_concs;				// Length = length_concs
	IntView* d_mats;				// Length = length_mats
	DoubleView* d_unionized_energy_array;		// Length = length_unionized_energy_array
	IntView* d_index_grid;				// Length = length_index_grid
	PointView* d_nuclide_grid;			// Length = length_nuclide_grid
	int * num_nucs;
	double * concs;
	int * mats;
	double * unionized_energy_array;
	int * index_grid;
	NuclideGridPoint * nuclide_grid;
	int length_num_nucs;
	int length_concs;
	int length_mats;
	int length_unionized_energy_array;
	long length_index_grid;
	int length_nuclide_grid;
	IntView* d_max_num_nucs;
	int max_num_nucs;
	double * p_energy_samples;
	int length_p_energy_samples;
	int * mat_samples;
	int length_mat_samples;
} SimulationData;

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int(long a);
Inputs read_CLI( int argc, char * argv[] );
void print_CLI_error(void);
void print_inputs(Inputs in, int nprocs, int version);
int print_results( Inputs in, int mype, double runtime, int nprocs, unsigned long long vhash );
void binary_write( Inputs in, SimulationData SD );
SimulationData binary_read( Inputs in );

// Simulation.c
unsigned long long run_event_based_simulation(Inputs in, SimulationData SD, int mype, double* end);
unsigned long long run_history_based_simulation(Inputs in, SimulationData SD, int mype);
KOKKOS_INLINE_FUNCTION
void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
			   long n_gridpoints,
			   DoubleView  egrid, IntView  index_data,
			   PointView *  nuclide_grids,
			   long idx, double *  xs_vector, int grid_type, int hash_bins );
KOKKOS_INLINE_FUNCTION
void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, IntView  num_nucs,
                         DoubleView  concs,
                         DoubleView  egrid, IntView  index_data,
                         PointView  nuclide_grids,
                         IntView  mats,
                         double *  macro_xs_vector, int grid_type, int hash_bins, int max_num_nucs );
KOKKOS_INLINE_FUNCTION
long grid_search( long n, double quarry, DoubleView A, long off);
KOKKOS_INLINE_FUNCTION
long grid_search_nuclide( long n, double quarry, PointView A, long off, long low, long high);
long grid_search_nuclide_old( long n, double quarry, NuclideGridPoint* A, long low, long high);
KOKKOS_INLINE_FUNCTION
int pick_mat( uint64_t * seed );
KOKKOS_FUNCTION
double LCG_random_double(uint64_t * seed);
KOKKOS_INLINE_FUNCTION
uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);
unsigned long long run_event_based_simulation_optimization_1(Inputs in, SimulationData SD, int mype);

// GridInit.c
SimulationData grid_init_do_not_profile( Inputs in, int mype );

// XSutils.c
int NGP_compare( const void * a, const void * b );
int double_compare(const void * a, const void * b);
size_t estimate_mem_usage( Inputs in );

// Materials.c
int * load_num_nucs(long n_isotopes);
int * load_mats( int * num_nucs, long n_isotopes, int * max_num_nucs );
double * load_concs( int * num_nucs, int max_num_nucs );
#endif
