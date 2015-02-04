#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main( int argc, char* argv[] )
{
  // =====================================================================
  // Declarations
  // =====================================================================
  int version = 1;  

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

  occaKernelInfo lookupInfo;
  occaKernel lookup_touch, lookup_kernel;
  occaDevice device;

  occaMemory dev_num_nucs, dev_energy_grid, dev_grid_ptrs, dev_nuclide_vector,
             dev_mats, dev_mats_idx, dev_concs;

  occaMemory dev_V_sums;

  // =====================================================================
  // Initialize MPI
  // =====================================================================

#ifdef MPI
  MPI_Status stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#endif

  // =====================================================================
  // Read command-line input
  // =====================================================================

  // Process CLI Fields and print summary
  Inputs in = read_CLI( argc, argv );
  if( mype == 0 )
    print_inputs( in, nprocs, version );

  // =====================================================================
  // Initialize OCCA
  // =====================================================================

  lookupInfo = occaGenKernelInfo();
  occaKernelInfoAddDefine(lookupInfo, "inner_dim", occaLong(in.inner_dim));
  occaKernelInfoAddDefine(lookupInfo, "outer_dim", occaLong(in.outer_dim));
#ifdef VERIFICATION
  // occaKernelInfoAddDefine(lookupInfo, "VERIFICATION", occaInt(1));
#endif

  device = occaGetDevice(in.device_info);

  lookup_touch = occaBuildKernelFromSource(device, "hybridLookupKernel.okl","lookup_touch", lookupInfo);
  lookup_kernel = occaBuildKernelFromSource(device, in.kernel, "lookup_kernel", lookupInfo);

  // =====================================================================
  // Initialize RNG
  // =====================================================================

  // rand() is only used in the serial initialization stages.
  // A custom RNG is used in parallel portions.
#ifdef VERIFICATION
  srand(26);
#else
  srand(time(NULL));
#endif

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
  // Initialize verification arrays
  // =====================================================================

  V_sums = (double *) calloc( 5 * in.lookups, sizeof(double) );

  // =====================================================================
  // Allocate OCCA memory
  // =====================================================================

  printf("Allocating and copying to device memory...\n");
  // REMEMBER: memcopy is part of malloc (last arg gets copied to device)

  if (strcasecmp(in.mode, "OpenMP") == 0) {
    dev_num_nucs       = occaDeviceWrapMemory(device, num_nucs, 12*sizeof(int));
    dev_nuclide_vector = occaDeviceWrapMemory(device, nuclide_grids[0],
        in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint));
    dev_energy_grid    = occaDeviceWrapMemory(device, energy_grid,
        in.n_isotopes*in.n_gridpoints*sizeof(GridPoint));
    dev_grid_ptrs      = occaDeviceWrapMemory(device, grid_ptrs,
        in.n_isotopes*in.n_isotopes*in.n_gridpoints*sizeof(int));
    dev_mats           = occaDeviceWrapMemory(device, mats    , size_mats*sizeof(int));
    dev_mats_idx       = occaDeviceWrapMemory(device, mats_idx, 12*sizeof(int));
    dev_concs          = occaDeviceWrapMemory(device, concs   , size_mats*sizeof(double));
    dev_V_sums         = occaDeviceWrapMemory(device, V_sums  , 5*in.lookups*sizeof(double));
  }
  else {
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
  }

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

  // =====================================================================
  // Finalize and verify
  // =====================================================================

  printf("Copying from device memory...\n");
  // Device-to-host memcopy
  occaCopyMemToPtr(V_sums, dev_V_sums, 5*in.lookups*sizeof(double), 0);

  // Reduce sums
  for(int i = 0; i < in.lookups; ++i){
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
