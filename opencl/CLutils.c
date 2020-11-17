#include "XSbench_header.h"

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
  if( error != 0 )
    printf("%s\n", getErrorString(error));
}

void printCompilerError( cl_program program, cl_device_id device )
{
  cl_int status;

  size_t logSize;
  char * log;

  status = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);

  check( status );

  log = (char *) malloc(logSize);
  if( !log) {
    exit(-1);
  }

  status = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
  check(status);

  if( strlen(log) > 0 )
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
  printf("Number of platforms detected = %d\n", platformCount);
  platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
  clGetPlatformIDs(platformCount, platforms, NULL);
  for (i = 0; i < platformCount; i++) {
    printf("Platorm %d:\n", i);
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
      printf("\tDevice %d: %s\n", j+1, value);
      free(value);
      // print hardware device version
      clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
      printf("\t\t%d.%d Hardware version: %s\n", j+1, 1, value);
      free(value);
      // print software driver version
      clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
      printf("\t\t%d.%d Software version: %s\n", j+1, 2, value);
      free(value);
      // print c version supported by compiler for device
      clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
      value = (char*) malloc(valueSize);
      clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
      printf("\t\t%d.%d OpenCL C version: %s\n", j+1, 3, value);
      free(value);
      // print parallel compute units
      clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
          sizeof(maxComputeUnits), &maxComputeUnits, NULL);
      printf("\t\t%d.%d Parallel compute units: %d\n", j+1, 4, maxComputeUnits);
    }
    free(devices);
  }
  free(platforms);  
}

OpenCLInfo initialize_device(int user_platform_id, int user_device_id)
{
  border_print();
  center_print("AVAILABLE OPENCL PLATFORMS & DEVICES", 79);
  border_print();
  print_opencl_info();
  border_print();

  center_print("DEVICE INITIALIZATION", 79);
  border_print();

  // Get platform and device information
  /*
  cl_platform_id platform_id = NULL;
  cl_device_id device_id = NULL;   
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  check(ret);
  ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
  check(ret);
  */

  int platform_idx = 0;
  if( user_platform_id != -1 )
    platform_idx = user_platform_id;

  // Get # of Platforms
  cl_uint ret_num_platforms;
  cl_int ret = clGetPlatformIDs(0, NULL, &ret_num_platforms);
  check(ret);

  // Allocate space for platform information
  cl_platform_id * platform_ids = (cl_platform_id *) malloc(ret_num_platforms * sizeof(cl_platform_id));

  // Fill in data on Platforms
  ret = clGetPlatformIDs(ret_num_platforms, platform_ids, NULL);

  cl_platform_id target_platform_id = platform_ids[platform_idx];   // '1' is the target platform index. needs to be changed accordingly.
  
  int device_idx = CL_DEVICE_TYPE_DEFAULT;
  if( user_device_id != -1 )
    device_idx = user_device_id;
  
  cl_uint ret_num_devices;
  cl_device_id device_id = NULL;   
  ret = clGetDeviceIDs( target_platform_id, device_idx, 1, &device_id, &ret_num_devices);

  // Print info about where we are running
  printf("Selected Device (platform id = %d, device id = %d)\n", platform_idx, device_idx);
  print_single_info(target_platform_id, device_id);
  printf("Note: platform ID can be specified with the \"-P\" flag and device id with the \"-D\" flag\n");

  // Create an OpenCL context
  cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
  check(ret);

  // Create a command queue
  //cl_command_queue command_queue = clCreateCommandQueueWithProperties(context, device_id, 0, &ret);
  cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
  check(ret);

  OpenCLInfo CL;
  CL.platform_id = target_platform_id;
  CL.device_id = device_id;
  CL.context = context;
  CL.command_queue = command_queue;
  return CL;
}

cl_mem copy_array_to_device(OpenCLInfo * CL, cl_mem_flags mem_flags, void * array, size_t sz)
{
  cl_int ret;
  cl_mem d_array = clCreateBuffer(CL->context, mem_flags,  sz, NULL, &ret);
  check(ret);
  ret = clEnqueueWriteBuffer(CL->command_queue, d_array, CL_TRUE, 0, sz, array, 0, NULL, NULL);
  check(ret);
  return d_array;
}

void copy_array_from_device(OpenCLInfo * CL, cl_mem * d_array, void * h_array, size_t sz)
{
  cl_int ret;
  ret = clEnqueueReadBuffer(CL->command_queue, *d_array, CL_TRUE, 0, sz, h_array, 0, NULL, NULL);
  check(ret);
}

void set_kernel_arguments(cl_kernel * kernel, int argc, size_t * arg_sz, void ** args)
{
  cl_int ret;
  for( int i = 0; i < argc; i++ )
  {
    ret = clSetKernelArg(*kernel, i, arg_sz[i], args[i]);
    check(ret);
  }
}

cl_kernel compile_kernel(OpenCLInfo * CL, char * kernel_name)
{
  printf("Compiling %s...\n", kernel_name);

  // Load the kernel source code into the array source_str
  FILE *fp;
  char *source_str;
  size_t source_size;
  char kernel_fname[256];
  strcpy(kernel_fname, kernel_name);
  strcat(kernel_fname, ".cl");

  fp = fopen(kernel_fname, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  source_str = (char*) malloc(MAX_SOURCE_SIZE);
  source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
  assert(source_size > 0 );
  fclose( fp );

  // Create a program from the kernel source
  cl_int ret;
  cl_program program = clCreateProgramWithSource(CL->context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
  check(ret);

  // Build the program
  ret = clBuildProgram(program, 1, &CL->device_id, NULL, NULL, NULL);
  check(ret);

  printCompilerError( program, CL->device_id );

  // Create the OpenCL kernel
  cl_kernel kernel = clCreateKernel(program, kernel_name, &ret);
  check(ret);

  return kernel;
} 

void clear_array(OpenCLInfo * CL, cl_mem * buffer, size_t sz)
{
  float fill = 0.0;
  cl_int ret = clEnqueueFillBuffer(
      CL->command_queue,
      *buffer,
      (void *) &fill,
      sizeof(float),
      0,
      sz,
      0,
      NULL,
      NULL);
  check(ret);
}

