#include "XSbench_header.h"

////////////////////////////////////////////////////////////////////////////////////
// BASELINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// All "baseline" code is at the top of this file. The baseline code is a simple
// implementation of the algorithm, with only minor CPU optimizations in place.
// Following these functions are a number of optimized variants,
// which each deploy a different combination of optimizations strategies. By
// default, XSBench will only run the baseline implementation. Optimized variants
// have not yet neen implemented in this OpenCL port.
////////////////////////////////////////////////////////////////////////////////////

unsigned long long run_event_based_simulation(Inputs in, SimulationData SD, int mype, double * sim_runtime)
{
	if( mype == 0)	
		printf("Initializing OpenCL data structures and JIT compiling kernel...\n");

	double start = get_time();

	int * verification_array_host = (int *) malloc( in.lookups * sizeof(int));
	
	////////////////////////////////////////////////////////////////////////////////
	// OpenCL Boilerplate Setup
	////////////////////////////////////////////////////////////////////////////////

	// Let's start setting up our openCL boilerplate
	// Load the kernel source code into the array source_str
	FILE *fp;
	char *source_str;
	size_t source_size;

	fp = fopen("kernel.cl", "r");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel.\n");
		exit(1);
	}
	source_str = (char*) malloc(MAX_SOURCE_SIZE);
	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose( fp );

	// Get platform and device information
	cl_platform_id platform_id = NULL;
	cl_device_id device_id = NULL;   
	cl_uint ret_num_devices;
	cl_uint ret_num_platforms;
	cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
	check(ret);
	ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
	//ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
	//ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
	check(ret);

	// Print info about where we are running
	print_single_info(platform_id, device_id);

	// Create an OpenCL context
	cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
	check(ret);

	// Create a command queue
	cl_command_queue command_queue = clCreateCommandQueueWithProperties(context, device_id, 0, &ret);
	check(ret);

	////////////////////////////////////////////////////////////////////////////////
	// OpenCL Move Memory To Device Buffers
	////////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////////
	// SUMMARY: Simulation Data Structure Manifest for "SD" Object
	// Here we list all heap arrays (and lengths) in SD that would need to be
	// offloaded manually if using an accelerator with a seperate memory space
	////////////////////////////////////////////////////////////////////////////////
	// int * num_nucs;                     // Length = length_num_nucs;
	// double * concs;                     // Length = length_concs
	// int * mats;                         // Length = length_mats
	// double * unionized_energy_array;    // Length = length_unionized_energy_array
	// int * index_grid;                   // Length = length_index_grid
	// NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
	// 
	// Note: "unionized_energy_array" and "index_grid" can be of zero length
	//        depending on lookup method.
	//
	// Note: "Lengths" are given as the number of objects in the array, not the
	//       number of bytes.
	////////////////////////////////////////////////////////////////////////////////
	
	// Create memory buffers on the device for each vector and move data over
	size_t sz = SD.length_num_nucs * sizeof(int);
	cl_mem num_nucs_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	check(ret);
	ret = clEnqueueWriteBuffer(command_queue, num_nucs_d, CL_TRUE, 0, sz, SD.num_nucs, 0, NULL, NULL);
	check(ret);

	sz = SD.length_concs * sizeof(double);
	cl_mem concs_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	check(ret);
	ret = clEnqueueWriteBuffer(command_queue, concs_d, CL_TRUE, 0, sz, SD.concs, 0, NULL, NULL);
	check(ret);
	
	sz = SD.length_mats * sizeof(int);
	cl_mem mats_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	check(ret);
	ret = clEnqueueWriteBuffer(command_queue, mats_d, CL_TRUE, 0, sz, SD.mats, 0, NULL, NULL);
	check(ret);
	
	// This buffer is not used if we are using the nuclide grid or hash grid methods,
	// so only allocate it if we need it.
	sz = SD.length_unionized_energy_array * sizeof(double);
	cl_mem unionized_energy_array_d;
	if( sz > 0 )
	{
		unionized_energy_array_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
		check(ret);
		ret = clEnqueueWriteBuffer(command_queue, unionized_energy_array_d, CL_TRUE, 0, sz, SD.unionized_energy_array, 0, NULL, NULL);
		check(ret);
	}
	
	// This buffer is not used if we are using the nuclide grid only method,
	// so only allocate it if we need it.
	sz = SD.length_index_grid * sizeof(int);
	cl_mem index_grid_d;
	if( sz > 0 )
	{
		// If using the unionized grid, this will be our largest allocation. As OpenCL
		// devices can have (for no good reason) very small maximum allocation sizes,
		// we need to check if our allocation will go over the limit ourselves. Currently,
		// this is not policed for us, OpenCL gives CL_SUCCESS and returns a non-zero ptr
		// even though it hasn't succeeded, so we check manually.
		cl_ulong max_opencl_allocation_size;
		clGetDeviceInfo(device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_opencl_allocation_size), &max_opencl_allocation_size, NULL);

		cl_ulong index_grid_requested_allocation_size = sz;
		assert( index_grid_requested_allocation_size <= max_opencl_allocation_size);

		index_grid_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
		check(ret);
		assert( ret == CL_SUCCESS);
		assert( index_grid_d != NULL);
		ret = clEnqueueWriteBuffer(command_queue, index_grid_d, CL_TRUE, 0, sz, SD.index_grid, 0, NULL, NULL);
		check(ret);
	}
	
	sz = SD.length_nuclide_grid * sizeof(NuclideGridPoint);
	cl_mem nuclide_grid_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	check(ret);
	ret = clEnqueueWriteBuffer(command_queue, nuclide_grid_d, CL_TRUE, 0, sz, SD.nuclide_grid, 0, NULL, NULL);
	check(ret);
	
	sz = in.lookups * sizeof(int);
	cl_mem verification_array = clCreateBuffer(context, CL_MEM_READ_WRITE,  sz, NULL, &ret);
	check(ret);
	
	////////////////////////////////////////////////////////////////////////////////
	// OpenCL Build and Intiailize Kernel Program
	////////////////////////////////////////////////////////////////////////////////
	
	// Create a program from the kernel source
	cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
	check(ret);

	// Build the program
	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	check(ret);
	
	printCompilerError( program, device_id );

	// Create the OpenCL kernel
	cl_kernel kernel = clCreateKernel(program, "macro_xs_lookup_kernel", &ret);
	check(ret);

	// Set the arguments of the kernel
	ret = clSetKernelArg(kernel, 0, sizeof(Inputs), (void *)&in);
	check(ret);
	ret = clSetKernelArg(kernel, 1, sizeof(int), (void *)&SD.max_num_nucs);
	check(ret);
	ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&num_nucs_d);
	check(ret);
	ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&concs_d);
	check(ret);
	if( SD.length_unionized_energy_array > 0 )
	{
		ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&unionized_energy_array_d);
		check(ret);
	}
	else
	{
		ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), NULL);
		check(ret);
	}
	if( SD.length_index_grid > 0 )
	{
		ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&index_grid_d);
		check(ret);
	}
	else
	{
		ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), NULL);
		check(ret);
	}
	ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&nuclide_grid_d);
	check(ret);
	ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mats_d);
	check(ret);
	ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&verification_array);
	check(ret);
	
	double stop = get_time();
	if( mype == 0) printf("OpenCL initialization time: %.3lf seconds\n", stop-start);
	start = stop;
	
	////////////////////////////////////////////////////////////////////////////////
	// Run Simulation Kernel
	////////////////////////////////////////////////////////////////////////////////
	
	if( mype == 0) printf("Running event based simulation...\n");

	// Execute the OpenCL kernel on the list
	size_t global_item_size = in.lookups; // Process the entire lists
	size_t local_item_size = 64; // Divide work items into groups of 64
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
	check(ret);
	
	////////////////////////////////////////////////////////////////////////////////
	// Retrieve verification data from device and reduce it
	////////////////////////////////////////////////////////////////////////////////
	
	// Read the memory buffer C on the device to the local variable C
	ret = clEnqueueReadBuffer(command_queue, verification_array, CL_TRUE, 0, in.lookups * sizeof(int), verification_array_host, 0, NULL, NULL);
	check(ret);
	
	if( mype == 0) printf("Reducing verification value...\n");
	
	unsigned long long verification = 0;

	for( int l = 0; l < in.lookups; l++ )
		verification += verification_array_host[l];

	stop = get_time();
	*sim_runtime = stop-start;
	if( mype == 0) printf("Simulation + Verification Reduction Runtime: %.3lf seconds\n", *sim_runtime);
	
	////////////////////////////////////////////////////////////////////////////////
	// OpenCL cleanup
	////////////////////////////////////////////////////////////////////////////////

	ret = clFlush(command_queue);
	check(ret);
	ret = clFinish(command_queue);
	check(ret);
	ret = clReleaseKernel(kernel);
	check(ret);
	ret = clReleaseProgram(program);
	check(ret);
	ret = clReleaseMemObject(num_nucs_d);
	check(ret);
	if( SD.length_unionized_energy_array > 0)
	{
		ret = clReleaseMemObject(unionized_energy_array_d);
		check(ret);
	}
	if( SD.length_index_grid > 0)
	{
		ret = clReleaseMemObject(index_grid_d);
		check(ret);
	}
	ret = clReleaseMemObject(nuclide_grid_d);
	check(ret);
	ret = clReleaseMemObject(mats_d);
	check(ret);
	ret = clReleaseMemObject(verification_array);
	check(ret);
	ret = clReleaseCommandQueue(command_queue);
	check(ret);
	ret = clReleaseContext(context);
	check(ret);

	return verification;
}

// binary search for energy on nuclide energy grid
long grid_search_nuclide( long n, double quarry, NuclideGridPoint * A, long low, long high)
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

double LCG_random_double(uint64_t * seed)
{
	// LCG parameters
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
}	
