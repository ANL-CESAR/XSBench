#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

const long outer_dim = 128;
const long inner_dim = 128;  

int main( int argc, char* argv[] )
{
  // =====================================================================
  // Initialization & Command Line Read-In
  // =====================================================================
  int version = 13;
  int mype = 0;
  int max_procs = omp_get_num_procs();
  int i, thread, mat;
  unsigned long seed;
  double omp_start, omp_end, p_energy;
  unsigned long long vhash = 0;
  int nprocs;

  // Timing
  struct timeval start, end;
  double wall_time;

  const char *mode = "CUDA";
  int platformID = 0;
  int deviceID   = 0;

  occaKernel lookup_kernel;
  occaDevice device;

  occaMemory dev_energy_grid, dev_grid_ptrs, dev_nuclide_vector, dev_mats,
             dev_mats_idx, dev_concs;

  occaKernelInfo lookupInfo = occaGenKernelInfo();
  occaKernelInfoAddDefine(lookupInfo, "inner_dim", occaLong(inner_dim));
  occaKernelInfoAddDefine(lookupInfo, "outer_dim", occaLong(outer_dim));

  device = occaGetDevice(mode, platformID, deviceID);
  lookup_kernel = occaBuildKernelFromSource(device, "lookup_kernel.okl",
      "lookup_kernel", lookupInfo);


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

  // Set number of OpenMP Threads
  omp_set_num_threads(in.nthreads); 

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

  dev_nuclide_vector = occaDeviceMalloc(device,
      in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint), NULL);
  dev_energy_grid = occaDeviceMalloc(device,
      in.n_isotopes*in.n_gridpoints*sizeof(GridPoint), NULL);
  dev_grid_ptrs = occaDeviceMalloc(device,
      in.n_isotopes*in.n_isotopes*in.n_gridpoints*sizeof(int), NULL);
  dev_mats = occaDeviceMalloc(device, size_mats*sizeof(int), NULL);
  dev_mats_idx = occaDeviceMalloc(device, 12*sizeof(int), NULL);
  dev_concs = occaDeviceMalloc(device, size_mats*sizeof(double), NULL);

  occaCopyPtrToMem(dev_nuclide_vector, nuclide_grids[0],
      in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint), 0);
  occaCopyPtrToMem(dev_energy_grid, energy_grid,
      in.n_isotopes*in.n_gridpoints*sizeof(GridPoint), 0);
  occaCopyPtrToMem(dev_grid_ptrs, grid_ptrs,
      in.n_isotopes*in.n_isotopes*in.n_gridpoints*sizeof(int), 0);
  occaCopyPtrToMem(dev_mats, mats, size_mats*sizeof(int), 0);
  occaCopyPtrToMem(dev_mats_idx, mats_idx, 12*sizeof(int), 0);
  occaCopyPtrToMem(dev_concs, concs, size_mats*sizeof(double), 0);

  // =====================================================================
  // Cross Section (XS) Parallel Lookup Simulation Begins
  // =====================================================================

  // Begin timer
  occaDeviceFinish(device);
  gettimeofday(&start, NULL);

  occaKernelRun( lookup_kernel,
    dev_energy_grid,
    dev_grid_ptrs,
    dev_nuclide_vector,
    dev_mats,
    dev_mats_idx,
    dev_concs,
    occaLong(in.lookups),
    occaLong(in.n_isotopes),
    occaLong(in.n_gridpoints)
    );

  occaDeviceFinish(device);
  gettimeofday(&end, NULL);
  wall_time = (end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec);
  
  printf("\n" );
  printf("Simulation complete.\n" );

  // Print / Save Results and Exit
  print_results( in, mype, wall_time, nprocs, vhash );

#ifdef BENCHMARK
}
#endif

#ifdef MPI
MPI_Finalize();
#endif

return 0;
}
