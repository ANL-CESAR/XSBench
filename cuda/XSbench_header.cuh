#ifndef __XSBENCH_HEADER_H__
#define __XSBENCH_HEADER_H__

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<cuda.h>
#include <thrust/reduce.h>
#include<stdint.h>
#include <chrono> 

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

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

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

typedef struct{
	int * num_nucs;                     // Length = length_num_nucs;
	double * concs;                     // Length = length_concs
	int * mats;                         // Length = length_mats
	double * unionized_energy_array;    // Length = length_unionized_energy_array
	int * index_grid;                   // Length = length_index_grid
	NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
	int length_num_nucs;
	int length_concs;
	int length_mats;
	int length_unionized_energy_array;
	long length_index_grid;
	int length_nuclide_grid;
	int max_num_nucs;
	unsigned long * verification;
	int length_verification;
	double * p_energy_samples;
	int length_p_energy_samples;
	int * mat_samples;
	int length_mat_samples;
} SimulationData;

// io.cu
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

// Simulation.cu
unsigned long long run_event_based_simulation_baseline(Inputs in, SimulationData SD, int mype);
__global__ void xs_lookup_kernel_baseline(Inputs in, SimulationData GSD );
__device__ void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           double * __restrict__ egrid, int * __restrict__ index_data,
                           NuclideGridPoint * __restrict__ nuclide_grids,
                           long idx, double * __restrict__ xs_vector, int grid_type, int hash_bins );
__device__ void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, int * __restrict__ num_nucs,
                         double * __restrict__ concs,
                         double * __restrict__ egrid, int * __restrict__ index_data,
                         NuclideGridPoint * __restrict__ nuclide_grids,
                         int * __restrict__ mats,
                         double * __restrict__ macro_xs_vector, int grid_type, int hash_bins, int max_num_nucs );
__device__ long grid_search( long n, double quarry, double * __restrict__ A);
__host__ __device__ long grid_search_nuclide( long n, double quarry, NuclideGridPoint * A, long low, long high);
__device__ int pick_mat( uint64_t * seed );
__host__ __device__ double LCG_random_double(uint64_t * seed);
__device__ uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);

unsigned long long run_event_based_simulation_optimization_1(Inputs in, SimulationData GSD, int mype);
__global__ void sampling_kernel(Inputs in, SimulationData GSD );
__global__ void xs_lookup_kernel_optimization_1(Inputs in, SimulationData GSD );

unsigned long long run_event_based_simulation_optimization_2(Inputs in, SimulationData GSD, int mype);
__global__ void xs_lookup_kernel_optimization_2(Inputs in, SimulationData GSD, int m );

unsigned long long run_event_based_simulation_optimization_3(Inputs in, SimulationData GSD, int mype);
__global__ void xs_lookup_kernel_optimization_3(Inputs in, SimulationData GSD, int m );

unsigned long long run_event_based_simulation_optimization_4(Inputs in, SimulationData GSD, int mype);
__global__ void xs_lookup_kernel_optimization_4(Inputs in, SimulationData GSD, int m, int n_lookups, int offset );

unsigned long long run_event_based_simulation_optimization_5(Inputs in, SimulationData GSD, int mype);
__global__ void xs_lookup_kernel_optimization_5(Inputs in, SimulationData GSD, int n_lookups, int offset );

unsigned long long run_event_based_simulation_optimization_6(Inputs in, SimulationData GSD, int mype);

// GridInit.cu
SimulationData grid_init_do_not_profile( Inputs in, int mype );
SimulationData move_simulation_data_to_device( Inputs in, int mype, SimulationData SD );

// XSutils.cu
int NGP_compare( const void * a, const void * b );
int double_compare(const void * a, const void * b);
double rn_v(void);
size_t estimate_mem_usage( Inputs in );
double get_time(void);

// Materials.cu
int * load_num_nucs(long n_isotopes);
int * load_mats( int * num_nucs, long n_isotopes, int * max_num_nucs );
double * load_concs( int * num_nucs, int max_num_nucs );
#endif
