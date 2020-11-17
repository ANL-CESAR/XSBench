
typedef struct{
  cl_platform_id platform_id;
  cl_device_id device_id;
  cl_context context;
  cl_command_queue command_queue;
} OpenCLInfo;

const char *getErrorString(cl_int error);
void check(cl_int error);
void printCompilerError( cl_program program, cl_device_id device );
void print_single_info( cl_platform_id platform, cl_device_id device);
void print_opencl_info(void);
OpenCLInfo initialize_device(int user_platform_id, int user_device_id);
cl_mem copy_array_to_device(OpenCLInfo * CL, cl_mem_flags mem_flags, void * array, size_t sz);
void copy_array_from_device(OpenCLInfo * CL, cl_mem * d_array, void * h_array, size_t sz);
void set_kernel_arguments(cl_kernel * kernel, int argc, size_t * arg_sz, void ** args);
cl_kernel compile_kernel(OpenCLInfo * CL, char * kernel_name);
void clear_array(OpenCLInfo * CL, cl_mem * buffer, size_t sz);
