#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#include <thrust/system/omp/execution_policy.h>
#include <thrust/sort.h>
#include<omp.h>

extern "C"
void parallel_sort_by_material( int * mat_samples, double * p_energy_samples, int size )
{
	//thrust::sort(thrust::omp::par, H.begin(), H.end());
	thrust::sort_by_key(thrust::omp::par, mat_samples, mat_samples + size, p_energy_samples);
}

extern "C"
void parallel_sort_by_energy( int * mat_samples, double * p_energy_samples, int size )
{
	//thrust::sort(thrust::omp::par, H.begin(), H.end());
	thrust::sort_by_key(thrust::omp::par, p_energy_samples, p_energy_samples + size, mat_samples);
}
