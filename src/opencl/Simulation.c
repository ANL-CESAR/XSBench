#include "XSbench_header.h"

#define CL_TARGET_OPENCL_VERSION 200
#include <CL/cl.h>
#define MAX_SOURCE_SIZE (0x100000)

const char *getErrorString(cl_int error)
{
switch(error){
    // run-time and JIT compiler errors
    case 0: return "CL_SUCCESS";
    case -1: return "CL_DEVICE_NOT_FOUND";
    case -2: return "CL_DEVICE_NOT_AVAILABLE";
    case -3: return "CL_COMPILER_NOT_AVAILABLE";
    case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5: return "CL_OUT_OF_RESOURCES";
    case -6: return "CL_OUT_OF_HOST_MEMORY";
    case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8: return "CL_MEM_COPY_OVERLAP";
    case -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
    }
}

void check(cl_int error)
{
	printf("%s\n", getErrorString(error));
}

void printCompilerError( cl_program program, cl_device_id device )
{
	cl_int status;

	size_t logSize;
	char * log;

	status = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);

	//check( status );

	log = (char *) malloc(logSize);
	if( !log) {
		exit(-1);
	}

	status = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
	//check(status);

	printf("%s\n", log);
}

void print_single_info( cl_platform_id platform, cl_device_id device)
{
	char* value;
	size_t valueSize;
	cl_uint maxComputeUnits;
	// print device name
	clGetDeviceInfo(device, CL_DEVICE_NAME, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_NAME, valueSize, value, NULL);
	printf("Device: %s\n", value);
	free(value);
	// print hardware device version
	clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_VERSION, valueSize, value, NULL);
	printf(" %d Hardware version: %s\n", 1, value);
	free(value);
	// print software driver version
	clGetDeviceInfo(device, CL_DRIVER_VERSION, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DRIVER_VERSION, valueSize, value, NULL);
	printf(" %d Software version: %s\n", 2, value);
	free(value);
	// print c version supported by compiler for device
	clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
	printf(" %d OpenCL C version: %s\n", 3, value);
	free(value);
	// print parallel compute units
	clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,
			sizeof(maxComputeUnits), &maxComputeUnits, NULL);
	printf(" %d Parallel compute units: %d\n", 4, maxComputeUnits);
}

void print_opencl_info(void)
{
	int i, j;
	char* value;
	size_t valueSize;
	cl_uint platformCount;
	cl_platform_id* platforms;
	cl_uint deviceCount;
	cl_device_id* devices;
	cl_uint maxComputeUnits;
	// get all platforms
	clGetPlatformIDs(0, NULL, &platformCount);
	platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
	clGetPlatformIDs(platformCount, platforms, NULL);
	for (i = 0; i < platformCount; i++) {
		// get all devices
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
		devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
		// for each device print critical attributes
		for (j = 0; j < deviceCount; j++) {
			// print device name
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
			value = (char*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
			printf("%d. Device: %s\n", j+1, value);
			free(value);
			// print hardware device version
			clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
			value = (char*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
			printf(" %d.%d Hardware version: %s\n", j+1, 1, value);
			free(value);
			// print software driver version
			clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
			value = (char*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
			printf(" %d.%d Software version: %s\n", j+1, 2, value);
			free(value);
			// print c version supported by compiler for device
			clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
			value = (char*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
			printf(" %d.%d OpenCL C version: %s\n", j+1, 3, value);
			free(value);
			// print parallel compute units
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
					sizeof(maxComputeUnits), &maxComputeUnits, NULL);
			printf(" %d.%d Parallel compute units: %d\n", j+1, 4, maxComputeUnits);
		}
		free(devices);
	}
	free(platforms);	
}

////////////////////////////////////////////////////////////////////////////////////
// BASELINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
// All "baseline" code is at the top of this file. The baseline code is a simple
// implementation of the algorithm, with only minor CPU optimizations in place.
// Following these functions are a number of optimized variants,
// which each deploy a different combination of optimizations strategies. By
// default, XSBench will only run the baseline implementation. Optimized variants
// must be specifically selected using the "-k <optimized variant ID>" command
// line argument.
////////////////////////////////////////////////////////////////////////////////////

unsigned long long run_event_based_simulation(Inputs in, SimulationData SD, int mype)
{
	if( mype == 0)	
		printf("Beginning event based simulation...\n");

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
	ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
	//ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
	//ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);

	// Print info about where we are running
	print_single_info(platform_id, device_id);

	// Create an OpenCL context
	cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);

	// Create a command queue
	cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

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
	ret = clEnqueueWriteBuffer(command_queue, num_nucs_d, CL_TRUE, 0, sz, SD.num_nucs, 0, NULL, NULL);

	sz = SD.length_concs * sizeof(double);
	cl_mem concs_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	ret = clEnqueueWriteBuffer(command_queue, concs_d, CL_TRUE, 0, sz, SD.concs, 0, NULL, NULL);
	
	sz = SD.length_mats * sizeof(int);
	cl_mem mats_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	ret = clEnqueueWriteBuffer(command_queue, mats_d, CL_TRUE, 0, sz, SD.mats, 0, NULL, NULL);
	
	sz = SD.length_unionized_energy_array * sizeof(double);
	cl_mem unionized_energy_array_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	ret = clEnqueueWriteBuffer(command_queue, unionized_energy_array_d, CL_TRUE, 0, sz, SD.unionized_energy_array, 0, NULL, NULL);
	
	sz = SD.length_index_grid * sizeof(int);
	cl_mem index_grid_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	ret = clEnqueueWriteBuffer(command_queue, index_grid_d, CL_TRUE, 0, sz, SD.index_grid, 0, NULL, NULL);
	
	sz = SD.length_nuclide_grid * sizeof(double);
	cl_mem nuclide_grid_d = clCreateBuffer(context, CL_MEM_READ_ONLY,  sz, NULL, &ret);
	ret = clEnqueueWriteBuffer(command_queue, nuclide_grid_d, CL_TRUE, 0, sz, SD.nuclide_grid, 0, NULL, NULL);
	
	sz = in.lookups * sizeof(int);
	cl_mem verification_array = clCreateBuffer(context, CL_MEM_READ_WRITE,  sz, NULL, &ret);
	
	////////////////////////////////////////////////////////////////////////////////
	// OpenCL Prepare and Launch Kernel
	////////////////////////////////////////////////////////////////////////////////
	
	// Create a program from the kernel source
	cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

	// Build the program
	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	
	printCompilerError( program, device_id );

	// Create the OpenCL kernel
	cl_kernel kernel = clCreateKernel(program, "macro_xs_lookup_kernel", &ret);

	// Set the arguments of the kernel
	ret = clSetKernelArg(kernel, 0, sizeof(Inputs), (void *)&in);
	ret = clSetKernelArg(kernel, 1, sizeof(int), (void *)&SD.max_num_nucs);
	ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&num_nucs_d);
	ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&unionized_energy_array_d);
	ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&index_grid_d);
	ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&nuclide_grid_d);
	ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mats_d);
	ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&verification_array);

	// Execute the OpenCL kernel on the list
	size_t global_item_size = in.lookups; // Process the entire lists
	size_t local_item_size = 64; // Divide work items into groups of 64
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
	
	////////////////////////////////////////////////////////////////////////////////
	// Retrieve verification data from device & cleanup OpenCL objects
	////////////////////////////////////////////////////////////////////////////////
	
	// Read the memory buffer C on the device to the local variable C
	ret = clEnqueueReadBuffer(command_queue, verification_array, CL_TRUE, 0, in.lookups * sizeof(int), verification_array_host, 0, NULL, NULL);

	// Clean up
	ret = clFlush(command_queue);
	ret = clFinish(command_queue);
	ret = clReleaseKernel(kernel);
	ret = clReleaseProgram(program);
	ret = clReleaseMemObject(num_nucs_d);
	ret = clReleaseMemObject(unionized_energy_array_d);
	ret = clReleaseMemObject(index_grid_d);
	ret = clReleaseMemObject(nuclide_grid_d);
	ret = clReleaseMemObject(mats_d);
	ret = clReleaseMemObject(verification_array);
	ret = clReleaseCommandQueue(command_queue);
	ret = clReleaseContext(context);

	////////////////////////////////////////////////////////////////////////////////
	// Redice Verification Value
	////////////////////////////////////////////////////////////////////////////////
	unsigned long long verification = 0;

	for( int l = 0; l < in.lookups; l++ )
		verification += verification_array_host[l];

	return verification;
}

// Calculates the microscopic cross section for a given nuclide & energy
void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints,
                           double * restrict egrid, int * restrict index_data,
                           NuclideGridPoint * restrict nuclide_grids,
                           long idx, double * restrict xs_vector, int grid_type, int hash_bins ){
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
	
	//test
	/*	
	if( omp_get_thread_num() == 0 )
	{
		printf("Lookup: Energy = %lf, nuc = %d\n", p_energy, nuc);
		printf("e_h = %lf e_l = %lf\n", high->energy , low->energy);
		printf("xs_h = %lf xs_l = %lf\n", high->elastic_xs, low->elastic_xs);
		printf("total_xs = %lf\n\n", xs_vector[1]);
	}
	*/
	
}

// Calculates macroscopic cross section based on a given material & energy 
void calculate_macro_xs( double p_energy, int mat, long n_isotopes,
                         long n_gridpoints, int * restrict num_nucs,
                         double * restrict concs,
                         double * restrict egrid, int * restrict index_data,
                         NuclideGridPoint * restrict nuclide_grids,
                         int * restrict mats,
                         double * restrict macro_xs_vector, int grid_type, int hash_bins, int max_num_nucs ){
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
	
	//test
	/*
	for( int k = 0; k < 5; k++ )
		printf("Energy: %lf, Material: %d, XSVector[%d]: %lf\n",
		       p_energy, mat, k, macro_xs_vector[k]);
	*/
}


// (fixed) binary search for energy on unionized energy grid
// returns lower index
long grid_search( long n, double quarry, double * restrict A)
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

// picks a material based on a probabilistic distribution
int pick_mat( unsigned long * seed )
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

double LCG_random_double(uint64_t * seed)
{
	// LCG parameters
	const uint64_t m = 9223372036854775808ULL; // 2^63
	const uint64_t a = 2806196910506780709ULL;
	const uint64_t c = 1ULL;
	*seed = (a * (*seed) + c) % m;
	return (double) (*seed) / (double) m;
	//return ldexp(*seed, -63);

}	

uint64_t fast_forward_LCG(uint64_t seed, uint64_t n)
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

