#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

#define USING_OPENMP 0
#define USING_CUDA   !(USING_OPENMP)

#if USING_OPENMP
const long outer_dim = 128;
const long inner_dim = 128;
#elif USING_CUDA
 /* const long outer_dim = 937500; const long inner_dim = 16; */
 const long outer_dim = 468750; const long inner_dim = 32;
/* const long outer_dim = 234376; const long inner_dim = 64; */
/* const long outer_dim = 156251; const long inner_dim = 96; */
/* const long outer_dim = 117188; const long inner_dim = 128; */
#endif

int main( int argc, char* argv[] )
{
  // =====================================================================
  // Initialization & Command Line Read-In
  // =====================================================================
  int version = 13;
  int i;

  int mype = 0;
  int nprocs = -1;

  // Timing
  struct timeval start, end;
  double wall_time;

  //---V_sum-------------------------------------------------------------------
  // These are the sum of the function values, evaluated in kernel.
  // They are the cumulative result of the random lookups.

  // Vectors for sums of F(x_i).  Dimensions will be V_sums[0:outer_dim].
  // In kernel, Each outer unit j will reduce is results to V_sum[i].
  // In main, we will need to reduce V_sums to get a single V_sum
  double *V_sums;

  // Sum of all F(x_i) from kernel.  Outside of kernel, V_sums will be reduced
  // to get V_sum
  double V_sum[5] = {0, 0, 0, 0, 0};
  //---------------------------------------------------------------------------

  //---OCCA declarations--------------------------------------------------------
#if USING_OPENMP
  const char *device_infos = "mode = OpenMP";
#elif USING_CUDA
  const char *device_infos = "mode = CUDA, deviceID = 0";
#endif

  occaKernel lookup_touch, lookup_kernel;
  occaDevice device;

  // For XSBench
  occaMemory dev_num_nucs, dev_energy_grid, dev_grid_ptrs, dev_nuclide_vector,
             dev_mats, dev_mats_idx, dev_concs;

  // For verification
  occaMemory dev_V_sums;

  occaKernelInfo lookupInfo = occaGenKernelInfo();
  occaKernelInfoAddDefine(lookupInfo, "inner_dim", occaLong(inner_dim));
  occaKernelInfoAddDefine(lookupInfo, "outer_dim", occaLong(outer_dim));
#ifdef VERIFICATION
  // occaKernelInfoAddDefine(lookupInfo, "VERIFICATION", occaInt(1));
#endif
  //---------------------------------------------------------------------------

  device = occaGetDevice(device_infos);

#if USING_OPENMP
  lookup_touch = occaBuildKernelFromSource(device,
                                           "lookup_kernel.okl","lookup_touch",
                                           lookupInfo);
  lookup_kernel = occaBuildKernelFromSource(device,
                                            "lookup_kernel.okl", "lookup_kernel",
                                            lookupInfo);
#elif USING_CUDA
  lookup_touch = occaBuildKernelFromSource(device,
                                           "cuda_lookup_kernel.okl", "lookup_touch",
                                           lookupInfo);
  lookup_kernel = occaBuildKernelFromSource(device,
                                            "cuda_lookup_kernel.okl", "lookup_kernel",
                                            lookupInfo);
#endif

#ifdef MPI
  MPI_Status stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#endif

  // rand() is only used in the serial initialization stages.
  // A custom RNG is used in parallel portions.
#ifdef VERIFICATION
  srand(26);
#else
  srand(time(NULL));
#endif

  // Process CLI Fields -- store in "Inputs" structure
  Inputs in = read_CLI( argc, argv );

  // Print-out of Input Summary
  if( mype == 0 )
    print_inputs( in, nprocs, version );

  // =====================================================================
  // Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
  // =====================================================================

  // Allocate & fill energy grids
#ifndef BINARY_READ
  if( mype == 0) printf("Generating Nuclide Energy Grids...\n");
#endif

  NuclideGridPoint ** nuclide_grids = gpmatrix(in.n_isotopes,in.n_gridpoints);

#ifdef VERIFICATION
  generate_grids_v( nuclide_grids, in.n_isotopes, in.n_gridpoints );
#else
  generate_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );
#endif

  // Sort grids by energy
#ifndef BINARY_READ
  if( mype == 0) printf("Sorting Nuclide Energy Grids...\n");
  sort_nuclide_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );
#endif

  // Prepare Unionized Energy Grid Framework
  int * grid_ptrs = generate_ptr_grid(in.n_isotopes, in.n_gridpoints);
#ifndef BINARY_READ
  GridPoint * energy_grid = generate_energy_grid( in.n_isotopes,
      in.n_gridpoints, nuclide_grids, grid_ptrs );
#else
  GridPoint * energy_grid = (GridPoint *)malloc( in.n_isotopes *
      in.n_gridpoints * sizeof( GridPoint ) );
  for( i = 0; i < in.n_isotopes*in.n_gridpoints; i++ )
    energy_grid[i].xs_ptrs = i*in.n_isotopes;
#endif

  // Double Indexing. Filling in energy_grid with pointers to the
  // nuclide_energy_grids.
#ifndef BINARY_READ
  set_grid_ptrs( energy_grid, nuclide_grids, grid_ptrs, in.n_isotopes, in.n_gridpoints );
#endif

#ifdef BINARY_READ
  if( mype == 0 ) printf("Reading data from \"XS_data.dat\" file...\n");
  binary_read(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
#endif

  // Get material data
  if( mype == 0 )
    printf("Loading Mats...\n");

  int size_mats;
  if (in.n_isotopes == 68)
    size_mats = 197;
  else
    size_mats = 484;

  int *num_nucs  = load_num_nucs(in.n_isotopes);
  int *mats_idx  = load_mats_idx(num_nucs);
  int *mats      = load_mats( num_nucs, mats_idx, size_mats, in.n_isotopes );

#ifdef VERIFICATION
  double *concs = load_concs_v(size_mats);
#else
  double *concs = load_concs(size_mats);
#endif

#ifdef BINARY_DUMP
  if( mype == 0 ) printf("Dumping data to binary file...\n");
  binary_dump(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
  if( mype == 0 ) printf("Binary file \"XS_data.dat\" written! Exiting...\n");
  return 0;
#endif

  // =====================================================================
  // Prepare verification arrays
  // =====================================================================

  V_sums = (double *) calloc( 5 * in.lookups, sizeof(double) );

  // =====================================================================
  // OCCA mallocs and memcopies
  // =====================================================================

  printf("Allocating and copying to device memory...\n");
  // REMEMBER: memcopy is part of malloc (last arg gets copied to device)
  dev_num_nucs       = occaDeviceMalloc(device, 12*sizeof(int), num_nucs);
  dev_nuclide_vector = occaDeviceMalloc(device,
                       in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint),
                       NULL);
  dev_energy_grid    = occaDeviceMalloc(device,
                       in.n_isotopes*in.n_gridpoints*sizeof(GridPoint),
                       NULL);
  dev_grid_ptrs      = occaDeviceMalloc(device,
                       in.n_isotopes*in.n_isotopes*in.n_gridpoints*sizeof(int),
                       NULL);
  dev_mats           = occaDeviceMalloc(device, size_mats*sizeof(int), mats);
  dev_mats_idx       = occaDeviceMalloc(device, 12*sizeof(int), mats_idx);
  dev_concs          = occaDeviceMalloc(device, size_mats*sizeof(double), concs);
  dev_V_sums         = occaDeviceMalloc(device, 5*in.lookups*sizeof(double), V_sums);

  // Call kernel to apply "proper" first-touch on large arrays
  occaKernelRun(lookup_touch,
                dev_energy_grid,
                dev_grid_ptrs,
                dev_nuclide_vector,
                occaLong(in.n_isotopes),
                occaLong(in.n_gridpoints));

  // Properly initialize arrays
  occaCopyPtrToMem(dev_nuclide_vector, nuclide_grids[0], occaAutoSize, occaNoOffset);
  occaCopyPtrToMem(dev_energy_grid   , energy_grid     , occaAutoSize, occaNoOffset);
  occaCopyPtrToMem(dev_grid_ptrs     , grid_ptrs       , occaAutoSize, occaNoOffset);

  // =====================================================================
  // Cross Section (XS) Parallel Lookup Simulation Begins
  // =====================================================================

  // Begin timer
  occaDeviceFinish(device);
  gettimeofday(&start, NULL);

  printf("Beginning kernel...\n");
  occaKernelRun( lookup_kernel,
      dev_num_nucs,
      dev_energy_grid,
      dev_grid_ptrs,
      dev_nuclide_vector,
      dev_mats,
      dev_mats_idx,
      dev_concs,
      occaLong(in.lookups),
      occaLong(in.n_isotopes),
      occaLong(in.n_gridpoints),
      dev_V_sums
      );
  printf("Kernel complete.\n" );

  occaDeviceFinish(device);
  gettimeofday(&end, NULL);
  wall_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.;

  printf("Copying from device memory...\n");
  // Device-to-host memcopy
  occaCopyMemToPtr(V_sums, dev_V_sums, 5*in.lookups*sizeof(double), 0);

  // Reduce sums
  for(i = 0; i < in.lookups; ++i){
    V_sum[0] += V_sums[5*i + 0];
    V_sum[1] += V_sums[5*i + 1];
    V_sum[2] += V_sums[5*i + 2];
    V_sum[3] += V_sums[5*i + 3];
    V_sum[4] += V_sums[5*i + 4];
  }

  const double V_total_sum = ( V_sum[0] +
                               V_sum[1] +
                               V_sum[2] +
                               V_sum[3] +
                               V_sum[4] );

  // Print / Save Results and Exit
  print_results( in, mype, wall_time, nprocs, V_total_sum ); //last arge should be vhahs

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
