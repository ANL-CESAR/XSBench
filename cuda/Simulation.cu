#include "XSbench_header.cuh"

////////////////////////////////////////////////////////////////////////////////////
// BASELINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// All "baseline" code is at the top of this file. The baseline code is a simple
// port of the original CPU OpenMP code to CUDA with few significant changes or
// optimizations made. Following these functions are a number of optimized variants,
// which each deploy a different combination of optimizations strategies. By
// default, XSBench will only run the baseline implementation. Optimized variants
// must be specifically selected using the "-k <optimized variant ID>" command
// line argument.
////////////////////////////////////////////////////////////////////////////////////

unsigned long long run_event_based_simulation_baseline(Inputs in, SimulationData GSD, int mype)
{
	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Running baseline event-based simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);

	xs_lookup_kernel_baseline<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

// In this kernel, we perform a single lookup with each thread. Threads within a warp
// do not really have any relation to each other, and divergence due to high nuclide count fuel
// material lookups are costly. This kernel constitutes baseline performance.
__global__ void xs_lookup_kernel_baseline(Inputs in, SimulationData GSD )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	const int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= in.lookups )
		return;

	// Set the initial seed value
	uint64_t seed = STARTING_SEED;	

	// Forward seed to lookup index (we need 2 samples per lookup)
	seed = fast_forward_LCG(seed, 2*i);

	// Randomly pick an energy and material for the particle
	double p_energy = LCG_random_double(&seed);
	int mat         = pick_mat(&seed); 
		
	double macro_xs_vector[5] = {0};
		
	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs(
			p_energy,        // Sampled neutron energy (in lethargy)
			mat,             // Sampled material type index neutron is in
			in.n_isotopes,   // Total number of isotopes in simulation
			in.n_gridpoints, // Number of gridpoints per isotope in simulation
			GSD.num_nucs,     // 1-D array with number of nuclides per material
			GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
			GSD.unionized_energy_array, // 1-D Unionized energy array
			GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
			GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
			GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
			macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
			in.grid_type,    // Lookup type (nuclide, hash, or unionized)
			in.hash_bins,    // Number of hash bins used (if using hash lookup type)
			GSD.max_num_nucs  // Maximum number of nuclides present in any material
			);

	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we have each thread
	// write to its thread_id index in an array, which we will reduce
	// with a thrust reduction kernel after the main simulation kernel.
	double max = -1.0;
	int max_idx = 0;
	for(int j = 0; j < 5; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	GSD.verification[i] = max_idx+1;
}

// Calculates the microscopic cross section for a given nuclide & energy
__device__ void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           double * __restrict__ egrid, int * __restrict__ index_data,
                           NuclideGridPoint * __restrict__ nuclide_grids,
                           long idx, double * __restrict__ xs_vector, int grid_type, int hash_bins ){
	// Variables
	double f;
	NuclideGridPoint * low, * high;

	// If using only the nuclide grid, we must perform a binary search
	// to find the energy location in this particular nuclide's grid.
	if( grid_type == NUCLIDE )
	{
		// Perform binary search on the Nuclide Grid to find the index
		idx = grid_search_nuclide( n_gridpoints, p_energy, &nuclide_grids[nuc*n_gridpoints], 0, n_gridpoints-1);

		// pull ptr from nuclide grid and check to ensure that
		// we're not reading off the end of the nuclide's grid
		if( idx == n_gridpoints - 1 )
			low = &nuclide_grids[nuc*n_gridpoints + idx - 1];
		else
			low = &nuclide_grids[nuc*n_gridpoints + idx];
	}
	else if( grid_type == UNIONIZED) // Unionized Energy Grid - we already know the index, no binary search needed.
	{
		// pull ptr from energy grid and check to ensure that
		// we're not reading off the end of the nuclide's grid
		if( index_data[idx * n_isotopes + nuc] == n_gridpoints - 1 )
			low = &nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc] - 1];
		else
			low = &nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc]];
	}
	else // Hash grid
	{
		// load lower bounding index
		int u_low = index_data[idx * n_isotopes + nuc];

		// Determine higher bounding index
		int u_high;
		if( idx == hash_bins - 1 )
			u_high = n_gridpoints - 1;
		else
			u_high = index_data[(idx+1)*n_isotopes + nuc] + 1;

		// Check edge cases to make sure energy is actually between these
		// Then, if things look good, search for gridpoint in the nuclide grid
		// within the lower and higher limits we've calculated.
		double e_low  = nuclide_grids[nuc*n_gridpoints + u_low].energy;
		double e_high = nuclide_grids[nuc*n_gridpoints + u_high].energy;
		int lower;
		if( p_energy <= e_low )
			lower = 0;
		else if( p_energy >= e_high )
			lower = n_gridpoints - 1;
		else
			lower = grid_search_nuclide( n_gridpoints, p_energy, &nuclide_grids[nuc*n_gridpoints], u_low, u_high);

		if( lower == n_gridpoints - 1 )
			low = &nuclide_grids[nuc*n_gridpoints + lower - 1];
		else
			low = &nuclide_grids[nuc*n_gridpoints + lower];
	}
	
	high = low + 1;
	
	// calculate the re-useable interpolation factor
	f = (high->energy - p_energy) / (high->energy - low->energy);

	// Total XS
	xs_vector[0] = high->total_xs - f * (high->total_xs - low->total_xs);
	
	// Elastic XS
	xs_vector[1] = high->elastic_xs - f * (high->elastic_xs - low->elastic_xs);
	
	// Absorbtion XS
	xs_vector[2] = high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs);
	
	// Fission XS
	xs_vector[3] = high->fission_xs - f * (high->fission_xs - low->fission_xs);
	
	// Nu Fission XS
	xs_vector[4] = high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs);
}

// Calculates macroscopic cross section based on a given material & energy 
__device__ void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, int * __restrict__ num_nucs,
                         double * __restrict__ concs,
                         double * __restrict__ egrid, int * __restrict__ index_data,
                         NuclideGridPoint * __restrict__ nuclide_grids,
                         int * __restrict__ mats,
                         double * __restrict__ macro_xs_vector, int grid_type, int hash_bins, int max_num_nucs ){
	int p_nuc; // the nuclide we are looking up
	long idx = -1;	
	double conc; // the concentration of the nuclide in the material

	// cleans out macro_xs_vector
	for( int k = 0; k < 5; k++ )
		macro_xs_vector[k] = 0;

	// If we are using the unionized energy grid (UEG), we only
	// need to perform 1 binary search per macroscopic lookup.
	// If we are using the nuclide grid search, it will have to be
	// done inside of the "calculate_micro_xs" function for each different
	// nuclide in the material.
	if( grid_type == UNIONIZED )
		idx = grid_search( n_isotopes * n_gridpoints, p_energy, egrid);	
	else if( grid_type == HASH )
	{
		double du = 1.0 / hash_bins;
		idx = p_energy / du;
	}
	
	// Once we find the pointer array on the UEG, we can pull the data
	// from the respective nuclide grids, as well as the nuclide
	// concentration data for the material
	// Each nuclide from the material needs to have its micro-XS array
	// looked up & interpolatied (via calculate_micro_xs). Then, the
	// micro XS is multiplied by the concentration of that nuclide
	// in the material, and added to the total macro XS array.
	// (Independent -- though if parallelizing, must use atomic operations
	//  or otherwise control access to the xs_vector and macro_xs_vector to
	//  avoid simulataneous writing to the same data structure)
	for( int j = 0; j < num_nucs[mat]; j++ )
	{
		double xs_vector[5];
		p_nuc = mats[mat*max_num_nucs + j];
		conc = concs[mat*max_num_nucs + j];
		calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		                    n_gridpoints, egrid, index_data,
		                    nuclide_grids, idx, xs_vector, grid_type, hash_bins );
		for( int k = 0; k < 5; k++ )
			macro_xs_vector[k] += xs_vector[k] * conc;
	}
}


// binary search for energy on unionized energy grid
// returns lower index
__device__ long grid_search( long n, double quarry, double * __restrict__ A)
{
	long lowerLimit = 0;
	long upperLimit = n-1;
	long examinationPoint;
	long length = upperLimit - lowerLimit;

	while( length > 1 )
	{
		examinationPoint = lowerLimit + ( length / 2 );
		
		if( A[examinationPoint] > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
}

// binary search for energy on nuclide energy grid
__host__ __device__ long grid_search_nuclide( long n, double quarry, NuclideGridPoint * A, long low, long high)
{
	long lowerLimit = low;
	long upperLimit = high;
	long examinationPoint;
	long length = upperLimit - lowerLimit;

	while( length > 1 )
	{
		examinationPoint = lowerLimit + ( length / 2 );
		
		if( A[examinationPoint].energy > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
}

// picks a material based on a probabilistic distribution
__device__ int pick_mat( uint64_t * seed )
{
	// I have a nice spreadsheet supporting these numbers. They are
	// the fractions (by volume) of material in the core. Not a 
	// *perfect* approximation of where XS lookups are going to occur,
	// but this will do a good job of biasing the system nonetheless.

	// Also could be argued that doing fractions by weight would be 
	// a better approximation, but volume does a good enough job for now.

	double dist[12];
	dist[0]  = 0.140;	// fuel
	dist[1]  = 0.052;	// cladding
	dist[2]  = 0.275;	// cold, borated water
	dist[3]  = 0.134;	// hot, borated water
	dist[4]  = 0.154;	// RPV
	dist[5]  = 0.064;	// Lower, radial reflector
	dist[6]  = 0.066;	// Upper reflector / top plate
	dist[7]  = 0.055;	// bottom plate
	dist[8]  = 0.008;	// bottom nozzle
	dist[9]  = 0.015;	// top nozzle
	dist[10] = 0.025;	// top of fuel assemblies
	dist[11] = 0.013;	// bottom of fuel assemblies
	
	double roll = LCG_random_double(seed);

	// makes a pick based on the distro
	for( int i = 0; i < 12; i++ )
	{
		double running = 0;
		for( int j = i; j > 0; j-- )
			running += dist[j];
		if( roll < running )
			return i;
	}

	return 0;
}

__host__ __device__ double LCG_random_double(uint64_t * seed)
{
	// LCG parameters
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
}	

__device__ uint64_t fast_forward_LCG(uint64_t seed, uint64_t n)
{
	// LCG parameters
	const uint64_t m = 9223372036854775808ULL; // 2^63
	uint64_t a = 2806196910506780709ULL;
	uint64_t c = 1ULL;

	n = n % m;

	uint64_t a_new = 1;
	uint64_t c_new = 0;

	while(n > 0) 
	{
		if(n & 1)
		{
			a_new *= a;
			c_new = c_new * a + c;
		}
		c *= (a + 1);
		a *= a;

		n >>= 1;
	}

	return (a_new * seed + c_new) % m;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZED VARIANT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// This section contains a number of optimized variants of some of the above
// functions, which each deploy a different combination of optimizations strategies
// specific to GPU. By default, XSBench will not run any of these variants. They
// must be specifically selected using the "-k <optimized variant ID>" command
// line argument.
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// Optimization 1 -- Basic kernel splitting of sampling & lookup routines
////////////////////////////////////////////////////////////////////////////////////
// This optimization requires a little extra data to store all material IDs and 
// energies for the sampled particles between kernel calls. By itself, this
// optimization is likely actually a bit of a slowdown compared to the baseline
// kernel. However, it will be used by better optimization kernels down the line.
////////////////////////////////////////////////////////////////////////////////////
unsigned long long run_event_based_simulation_optimization_1(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 1 - basic sample/lookup kernel splitting";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	xs_lookup_kernel_optimization_1<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

__global__ void sampling_kernel(Inputs in, SimulationData GSD )
{
	// The lookup ID.
	const int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= in.lookups )
		return;

	// Set the initial seed value
	uint64_t seed = STARTING_SEED;	

	// Forward seed to lookup index (we need 2 samples per lookup)
	seed = fast_forward_LCG(seed, 2*i);

	// Randomly pick an energy and material for the particle
	double p_energy = LCG_random_double(&seed);
	int mat         = pick_mat(&seed); 

	// Store sample data in state array
	GSD.p_energy_samples[i] = p_energy;
	GSD.mat_samples[i] = mat;
}

__global__ void xs_lookup_kernel_optimization_1(Inputs in, SimulationData GSD )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	const int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= in.lookups )
		return;
		
	double macro_xs_vector[5] = {0};
		
	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs(
			GSD.p_energy_samples[i],        // Sampled neutron energy (in lethargy)
			GSD.mat_samples[i],             // Sampled material type index neutron is in
			in.n_isotopes,   // Total number of isotopes in simulation
			in.n_gridpoints, // Number of gridpoints per isotope in simulation
			GSD.num_nucs,     // 1-D array with number of nuclides per material
			GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
			GSD.unionized_energy_array, // 1-D Unionized energy array
			GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
			GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
			GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
			macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
			in.grid_type,    // Lookup type (nuclide, hash, or unionized)
			in.hash_bins,    // Number of hash bins used (if using hash lookup type)
			GSD.max_num_nucs  // Maximum number of nuclides present in any material
			);

	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we have each thread
	// write to its thread_id index in an array, which we will reduce
	// with a thrust reduction kernel after the main simulation kernel.
	double max = -1.0;
	int max_idx = 0;
	for(int j = 0; j < 5; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	GSD.verification[i] = max_idx+1;
}

////////////////////////////////////////////////////////////////////////////////////
// Optimization 2 -- Kernel Splitting + Material-Specific Lookup Kernels
////////////////////////////////////////////////////////////////////////////////////
// This one builds on the first optimization. It uses multiple kernels, one
// for each material type, to better balance the workload across threads within
// a warp. This works because each material will have a different number of 
// isotopes, with some having a ton, meaning that SIMD efficiency can be rather
// low by default. Better efficiency may be gained in further optimizations by
// sorting the lookups first.
////////////////////////////////////////////////////////////////////////////////////
unsigned long long run_event_based_simulation_optimization_2(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 2 - Material Lookup Kernels";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// Launch all material kernels individually
	for( int m = 0; m < 12; m++ )
		xs_lookup_kernel_optimization_2<<<nblocks, nthreads>>>( in, GSD, m );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

__global__ void xs_lookup_kernel_optimization_2(Inputs in, SimulationData GSD, int m )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	const int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= in.lookups )
		return;
	
	// Check that our material type matches the kernel material
	int mat = GSD.mat_samples[i];
	if( mat != m )
		return;

	double macro_xs_vector[5] = {0};
		
	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs(
			GSD.p_energy_samples[i],        // Sampled neutron energy (in lethargy)
			mat,             // Sampled material type index neutron is in
			in.n_isotopes,   // Total number of isotopes in simulation
			in.n_gridpoints, // Number of gridpoints per isotope in simulation
			GSD.num_nucs,     // 1-D array with number of nuclides per material
			GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
			GSD.unionized_energy_array, // 1-D Unionized energy array
			GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
			GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
			GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
			macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
			in.grid_type,    // Lookup type (nuclide, hash, or unionized)
			in.hash_bins,    // Number of hash bins used (if using hash lookup type)
			GSD.max_num_nucs  // Maximum number of nuclides present in any material
			);

	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we have each thread
	// write to its thread_id index in an array, which we will reduce
	// with a thrust reduction kernel after the main simulation kernel.
	double max = -1.0;
	int max_idx = 0;
	for(int j = 0; j < 5; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	GSD.verification[i] = max_idx+1;
}


////////////////////////////////////////////////////////////////////////////////////
// Optimization 3 -- Kernel Splitting + Fuel or Not-Fuel Lookups
////////////////////////////////////////////////////////////////////////////////////
// This optimization alters Optimization 2. Instead of executing a kernel call for
// ALL different material types, only two different calls are made. One for fuel,
// and one for all the other materials. As the fuel material has by far the most
// isotopes, it takes much longer than the rest.
////////////////////////////////////////////////////////////////////////////////////
unsigned long long run_event_based_simulation_optimization_3(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 3 - Fuel or Other Lookup Kernels";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// Launch all material kernels individually
	xs_lookup_kernel_optimization_3<<<nblocks, nthreads>>>( in, GSD, 0 );
	xs_lookup_kernel_optimization_3<<<nblocks, nthreads>>>( in, GSD, 1 );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

__global__ void xs_lookup_kernel_optimization_3(Inputs in, SimulationData GSD, int is_fuel )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	const int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= in.lookups )
		return;
	
	int mat = GSD.mat_samples[i];

	// If this is the fuel kernel, AND this is a fuel lookup, then perform a lookup
	// OR if this is not the fuel kernel, AND this is not a fuel lookup, then perform the lookup
	if( ((is_fuel == 1) && (mat == 0)) || ((is_fuel == 0) && (mat != 0 ) ))
	{
		double macro_xs_vector[5] = {0};
			
		// Perform macroscopic Cross Section Lookup
		calculate_macro_xs(
				GSD.p_energy_samples[i],        // Sampled neutron energy (in lethargy)
				mat,             // Sampled material type index neutron is in
				in.n_isotopes,   // Total number of isotopes in simulation
				in.n_gridpoints, // Number of gridpoints per isotope in simulation
				GSD.num_nucs,     // 1-D array with number of nuclides per material
				GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
				GSD.unionized_energy_array, // 1-D Unionized energy array
				GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
				GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
				GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
				macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
				in.grid_type,    // Lookup type (nuclide, hash, or unionized)
				in.hash_bins,    // Number of hash bins used (if using hash lookup type)
				GSD.max_num_nucs  // Maximum number of nuclides present in any material
				);

		// For verification, and to prevent the compiler from optimizing
		// all work out, we interrogate the returned macro_xs_vector array
		// to find its maximum value index, then increment the verification
		// value by that index. In this implementation, we have each thread
		// write to its thread_id index in an array, which we will reduce
		// with a thrust reduction kernel after the main simulation kernel.
		double max = -1.0;
		int max_idx = 0;
		for(int j = 0; j < 5; j++ )
		{
			if( macro_xs_vector[j] > max )
			{
				max = macro_xs_vector[j];
				max_idx = j;
			}
		}
		GSD.verification[i] = max_idx+1;
	}
}


////////////////////////////////////////////////////////////////////////////////////
// Optimization 4 -- Kernel Splitting + All Material Lookups + Full Sort
////////////////////////////////////////////////////////////////////////////////////
// This optimization builds on optimization 2, adding in a full sort before
// hand so that the warps should be densely packed together. This should maximize
// SIMD efficiency of the kernel, but may incur an added cost for the sort.
////////////////////////////////////////////////////////////////////////////////////
unsigned long long run_event_based_simulation_optimization_4(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 4 - All Material Lookup Kernels + Material Sort";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// Count the number of fuel material lookups that need to be performed (fuel id = 0)
	int n_lookups_per_material[12];
	for( int m = 0; m < 12; m++ )
		n_lookups_per_material[m] = thrust::count(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, m);

	// Sort materials
	thrust::sort_by_key(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, GSD.p_energy_samples);
	
	// Launch all material kernels individually
	int offset = 0;
	for( int m = 0; m < 12; m++ )
	{
		nthreads = 32;
		nblocks = ceil((double) n_lookups_per_material[m] / (double) nthreads);
		xs_lookup_kernel_optimization_4<<<nblocks, nthreads>>>( in, GSD, m, n_lookups_per_material[m], offset );
		offset += n_lookups_per_material[m];
	}
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

__global__ void xs_lookup_kernel_optimization_4(Inputs in, SimulationData GSD, int m, int n_lookups, int offset )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= n_lookups )
		return;

	i += offset;

	// Check that our material type matches the kernel material
	int mat = GSD.mat_samples[i];
	if( mat != m )
		return;

	double macro_xs_vector[5] = {0};
		
	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs(
			GSD.p_energy_samples[i],        // Sampled neutron energy (in lethargy)
			mat,             // Sampled material type index neutron is in
			in.n_isotopes,   // Total number of isotopes in simulation
			in.n_gridpoints, // Number of gridpoints per isotope in simulation
			GSD.num_nucs,     // 1-D array with number of nuclides per material
			GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
			GSD.unionized_energy_array, // 1-D Unionized energy array
			GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
			GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
			GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
			macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
			in.grid_type,    // Lookup type (nuclide, hash, or unionized)
			in.hash_bins,    // Number of hash bins used (if using hash lookup type)
			GSD.max_num_nucs  // Maximum number of nuclides present in any material
			);

	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we have each thread
	// write to its thread_id index in an array, which we will reduce
	// with a thrust reduction kernel after the main simulation kernel.
	double max = -1.0;
	int max_idx = 0;
	for(int j = 0; j < 5; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	GSD.verification[i] = max_idx+1;
}

////////////////////////////////////////////////////////////////////////////////////
// Optimization 5 -- Kernel Splitting + Fuel/Other Lookups + Fuel/Other Partition
////////////////////////////////////////////////////////////////////////////////////
// This optimization is similar to optimization 4, but instead of sorting
// fully by material, we just sort by fuel or not fuel. Similarly, instead of
// launching kernels for all materials, similar to optimization 3 we only launch
// kernels for the fuel and other mateirals.
////////////////////////////////////////////////////////////////////////////////////

// Comparator for partitioning stage
struct is_mat_fuel{
	__host__ __device__
		bool operator()(const int & a)
		{
			return a == 0;
		}
};

unsigned long long run_event_based_simulation_optimization_5(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 5 - Fuel/No Fuel Lookup Kernels + Fuel/No Fuel Sort";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// Count the number of fuel material lookups that need to be performed (fuel id = 0)
	int n_fuel_lookups = thrust::count(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, 0);

	// Partition fuel into the first part of the array
	thrust::partition(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, GSD.p_energy_samples, is_mat_fuel());

	// Launch all material kernels individually (asynchronous is allowed)
	nblocks = ceil( (double) n_fuel_lookups / (double) nthreads);
	xs_lookup_kernel_optimization_5<<<nblocks, nthreads>>>( in, GSD, n_fuel_lookups, 0 );

	nblocks = ceil( (double) (in.lookups - n_fuel_lookups) / (double) nthreads);
	xs_lookup_kernel_optimization_5<<<nblocks, nthreads>>>( in, GSD, in.lookups-n_fuel_lookups, n_fuel_lookups );

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}

__global__ void xs_lookup_kernel_optimization_5(Inputs in, SimulationData GSD, int n_lookups, int offset )
{
	// The lookup ID. Used to set the seed, and to store the verification value
	int i = blockIdx.x *blockDim.x + threadIdx.x;

	if( i >= n_lookups )
		return;
	
	i += offset;

	double macro_xs_vector[5] = {0};
		
	// Perform macroscopic Cross Section Lookup
	calculate_macro_xs(
			GSD.p_energy_samples[i],        // Sampled neutron energy (in lethargy)
			GSD.mat_samples[i],             // Sampled material type index neutron is in
			in.n_isotopes,   // Total number of isotopes in simulation
			in.n_gridpoints, // Number of gridpoints per isotope in simulation
			GSD.num_nucs,     // 1-D array with number of nuclides per material
			GSD.concs,        // Flattened 2-D array with concentration of each nuclide in each material
			GSD.unionized_energy_array, // 1-D Unionized energy array
			GSD.index_grid,   // Flattened 2-D grid holding indices into nuclide grid for each unionized energy level
			GSD.nuclide_grid, // Flattened 2-D grid holding energy levels and XS_data for all nuclides in simulation
			GSD.mats,         // Flattened 2-D array with nuclide indices defining composition of each type of material
			macro_xs_vector, // 1-D array with result of the macroscopic cross section (5 different reaction channels)
			in.grid_type,    // Lookup type (nuclide, hash, or unionized)
			in.hash_bins,    // Number of hash bins used (if using hash lookup type)
			GSD.max_num_nucs  // Maximum number of nuclides present in any material
			);

	// For verification, and to prevent the compiler from optimizing
	// all work out, we interrogate the returned macro_xs_vector array
	// to find its maximum value index, then increment the verification
	// value by that index. In this implementation, we have each thread
	// write to its thread_id index in an array, which we will reduce
	// with a thrust reduction kernel after the main simulation kernel.
	double max = -1.0;
	int max_idx = 0;
	for(int j = 0; j < 5; j++ )
	{
		if( macro_xs_vector[j] > max )
		{
			max = macro_xs_vector[j];
			max_idx = j;
		}
	}
	GSD.verification[i] = max_idx+1;
}

////////////////////////////////////////////////////////////////////////////////////
// Optimization 6 -- Kernel Splitting + All Material Lookups + Full Sort
//                   + Energy Sort
////////////////////////////////////////////////////////////////////////////////////
// This optimization builds on optimization 4, adding in a second sort by energy.
// It is extremely fast, as now most of the threads within a warp will be hitting
// the same indices in the lookup grids. This greatly reduces thread divergence and
// greatly improves cache efficiency and re-use.
//
// However, it is unlikely that this exact optimization would be possible in a real
// application like OpenMC. One major difference is that particle objects are quite
// large, often having 50+ variable fields, such that sorting them in memory becomes
// rather expensive. Instead, the best possible option would probably be to create
// intermediate indexing (per Hamilton et. al 2019), and run the kernels indirectly.
////////////////////////////////////////////////////////////////////////////////////
unsigned long long run_event_based_simulation_optimization_6(Inputs in, SimulationData GSD, int mype)
{
	const char * optimization_name = "Optimization 6 - Material & Energy Sorts + Material-specific Kernels";
	
	if( mype == 0)	printf("Simulation Kernel:\"%s\"\n", optimization_name);
	
	////////////////////////////////////////////////////////////////////////////////
	// Allocate Additional Data Structures Needed by Optimized Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Allocating additional device data required by kernel...\n");
	size_t sz;
	size_t total_sz = 0;

	sz = in.lookups * sizeof(double);
	gpuErrchk( cudaMalloc((void **) &GSD.p_energy_samples, sz) );
	total_sz += sz;
	GSD.length_p_energy_samples = in.lookups;

	sz = in.lookups * sizeof(int);
	gpuErrchk( cudaMalloc((void **) &GSD.mat_samples, sz) );
	total_sz += sz;
	GSD.length_mat_samples = in.lookups;
	
	if( mype == 0)	printf("Allocated an additional %.0lf MB of data on GPU.\n", total_sz/1024.0/1024.0);

	////////////////////////////////////////////////////////////////////////////////
	// Configure & Launch Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Beginning optimized simulation...\n");

	int nthreads = 32;
	int nblocks = ceil( (double) in.lookups / 32.0);
	
	sampling_kernel<<<nblocks, nthreads>>>( in, GSD );
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// Count the number of fuel material lookups that need to be performed (fuel id = 0)
	int n_lookups_per_material[12];
	for( int m = 0; m < 12; m++ )
		n_lookups_per_material[m] = thrust::count(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, m);

	// Sort by material first
	thrust::sort_by_key(thrust::device, GSD.mat_samples, GSD.mat_samples + in.lookups, GSD.p_energy_samples);

	// Now, sort each material by energy
	int offset = 0;
	for( int m = 0; m < 12; m++ )
	{
		thrust::sort_by_key(thrust::device, GSD.p_energy_samples + offset, GSD.p_energy_samples + offset + n_lookups_per_material[m], GSD.mat_samples + offset);
		offset += n_lookups_per_material[m];
	}
	
	// Launch all material kernels individually
	offset = 0;
	for( int m = 0; m < 12; m++ )
	{
		nthreads = 32;
		nblocks = ceil((double) n_lookups_per_material[m] / (double) nthreads);
		xs_lookup_kernel_optimization_4<<<nblocks, nthreads>>>( in, GSD, m, n_lookups_per_material[m], offset );
		offset += n_lookups_per_material[m];
	}
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	////////////////////////////////////////////////////////////////////////////////
	// Reduce Verification Results
	////////////////////////////////////////////////////////////////////////////////
	if( mype == 0)	printf("Reducing verification results...\n");

	unsigned long verification_scalar = thrust::reduce(thrust::device, GSD.verification, GSD.verification + in.lookups, 0);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	return verification_scalar;
}
