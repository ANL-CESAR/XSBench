#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

#ifdef _OPENACC
#include "openacc.h"
#endif

int main( int argc, char* argv[] )
{
  // =====================================================================
  // Initialization & Command Line Read-In
  // =====================================================================
  int version = 13;             // Version number
  int thread, mat;              // OMP thread index, material index
  unsigned long seed;           // RNG seed for OMP version
  double tick, tock, p_energy;  // Start time, end time, particle energy
  double dval = 0;              // A dummy value, used to reduce xs values
  unsigned long long vhash = 0; // The verfication hash
  int nprocs;                   // Number of MPI procs
  int mype = 0;                 // MPI rank

  char HM[6];                   // Size of HM benchmark problem

  // Fractions (by volume) of materials in the reactor core.
  // These are used as probabilities to approximate where xs lookups will occur.
  double dist[12] = {           
    0.140,	// fuel
    0.052,	// cladding
    0.275,	// cold, borated water
    0.134,	// hot, borated water
    0.154,	// RPV
    0.064,	// Lower, radial reflector
    0.066,	// Upper reflector / top plate
    0.055,	// bottom plate
    0.008,	// bottom nozzle
    0.015,	// top nozzle
    0.025,	// top of fuel assemblies
    0.013 	// bottom of fuel assemblies
  };

#ifdef MPI
  MPI_Status stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#endif

#ifdef _OPENACC
  acc_init(acc_get_device_type());
#endif

  // rand() is only used in the serial initialization stages.
  // A custom RNG is used in parallel portions.
#ifdef VERIFICATION
  srand(26);
#else
  srand(time(NULL));
#endif

  // Process CLI Fields -- store in "Inputs" structure
  // Duplicate as constant values to resolve data dependencies in _OPENACC loops.
  Inputs in = read_CLI(argc, argv);
  const int nthreads = in.nthreads;
  const long n_isotopes = in.n_isotopes;
  const long n_gridpoints = in.n_gridpoints; 
  const int lookups = in.lookups;

  // Set number of OpenMP Threads
  omp_set_num_threads(in.nthreads);

  // Print-out of Input Summary
  if(mype == 0)
    print_inputs(in, nprocs, version);

  // =====================================================================
  // Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
  // =====================================================================

  // === Allocate & fill energy grids

#ifndef BINARY_READ
  if(mype == 0) printf("Generating Nuclide Energy Grids...\n");
#endif

  // allocates nuclide_grids[0:n_isotopes][0:n_gridpoints]
  NuclideGridPoint (* restrict nuclide_grids)[(long) n_gridpoints] = 
    (NuclideGridPoint (*)[(long) n_gridpoints]) 
    malloc(n_isotopes * n_gridpoints * sizeof(NuclideGridPoint));

  // fill grids deterministically or randomly
#ifdef VERIFICATION
  generate_grids_v(n_isotopes, n_gridpoints, nuclide_grids);
#else
  generate_grids(n_isotopes, n_gridpoints, nuclide_grids);	
#endif

  // Sort grids by energy
#ifndef BINARY_READ
  if(mype == 0) printf("Sorting Nuclide Energy Grids...\n");
  sort_nuclide_grids(n_isotopes, n_gridpoints, nuclide_grids);
#endif

  // Prepare Unionized Energy Grid Framework
  int * restrict grid_ptrs = generate_ptr_grid(n_isotopes, n_gridpoints);
#ifndef BINARY_READ
  GridPoint * restrict energy_grid = generate_energy_grid(n_isotopes,
      n_gridpoints, nuclide_grids, grid_ptrs);
#else
  GridPoint * restrict energy_grid = (GridPoint *)malloc(n_isotopes *
      n_gridpoints * sizeof(GridPoint));
  for(i = 0; i < n_isotopes*n_gridpoints; i++)
    energy_grid[i].xs_ptrs = i*n_isotopes;
#endif

  // Double Indexing. Filling in energy_grid with pointers to the
  // nuclide_energy_grids.
#ifndef BINARY_READ
  set_grid_ptrs(energy_grid, grid_ptrs, n_isotopes, n_gridpoints, nuclide_grids);
#endif

#ifdef BINARY_READ
  if(mype == 0) printf("Reading data from \"XS_data.dat\" file...\n");
  binary_read(n_isotopes, n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
#endif

  // Get material data
  if(mype == 0) printf("Loading Mats...\n");

  int size_mats;
  if (n_isotopes == 68) 
    size_mats = 197;
  else
    size_mats = 484;

  // The number of nuclides in each material
  int * restrict num_nucs = load_num_nucs(n_isotopes);
  // The indices of each material 
  int * restrict mats_idx = load_mats_idx(num_nucs);
  // The nuclide identities of each material
  int * restrict mats     = load_mats(num_nucs, mats_idx, size_mats, n_isotopes);

  // The concentrations of nuclides in each material
#ifdef VERIFICATION
  double * restrict concs = load_concs_v(size_mats);
#else
  double * restrict concs = load_concs(size_mats);
#endif

#ifdef _OPENACC

  // === In OMP, the random numbers are generated on the fly; and the results
  // are hashed on the fly.  In _OPENACC, we pre-generate random numbers and
  // post-process the results on the host.

  // Generate a stream of random numbers to copyin to device
  double * restrict rands = malloc(2*lookups*sizeof(double));
  for(int i=0; i<lookups; i++){
  #ifdef VERIFICATION
    rands[2*i] = rn_v();
    rands[2*i+1] = rn_v();
  #else
    rands[2*i] = (double) rand() / (double) RAND_MAX;
    rands[2*i+1] = (double) rand() / (double) RAND_MAX;
  #endif
  }

  // Allocate arrays for results to copyout from device
  #ifdef VERIFICATION
  int n_v_ints = lookups;
  int n_v_doubles = 6*lookups;
  int * restrict v_ints = malloc(n_v_ints*sizeof(int));
  double * restrict v_doubles = malloc(n_v_doubles*sizeof(double));
  #else
  int n_v_ints = 1;
  int n_v_doubles = 1;
  int * restrict v_ints = malloc(n_v_ints*sizeof(int));
  double * restrict v_doubles = malloc(n_v_doubles*sizeof(double));
  #endif

#endif

#ifdef BINARY_DUMP
  if(mype == 0) printf("Dumping data to binary file...\n");
  binary_dump(n_isotopes, n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
  if(mype == 0) printf("Binary file \"XS_data.dat\" written! Exiting...\n");
  return 0;
#endif

#ifdef _OPENACC
#pragma acc enter data \
  copyin( \
      vhash, \
      dval, \
      n_isotopes, \
      n_gridpoints, \
      lookups, \
      energy_grid[0:n_isotopes*n_gridpoints], \
      nuclide_grids[0:n_isotopes][0:n_gridpoints], \
      grid_ptrs[0:n_isotopes*n_isotopes*n_gridpoints], \
      mats[0:size_mats], \
      mats_idx[0:12], \
      concs[0:size_mats], \
      num_nucs[0:12], \
      dist[0:12], \
      rands[0:2*lookups] \
      ) \
  create( \
      v_ints[0:n_v_ints], \
      v_doubles[0:n_v_doubles] \
      ) 
#endif

  // =====================================================================
  // Cross Section (XS) Parallel Lookup Simulation Begins
  // =====================================================================

  if(mype == 0){
    printf("\n");
    border_print();
    center_print("SIMULATION", 79);
    border_print();
  }

#ifndef _OPENACC
  tick = omp_get_wtime();
#pragma omp parallel default(none) \
  private(thread, p_energy, mat, seed) \
  shared( in, energy_grid, nuclide_grids, \
      mats, concs, num_nucs, mype, vhash, \
      mats_idx, dist, grid_ptrs, dval)
#else
  tick = timer();
#pragma acc kernels \
  present( \
      vhash, \
      dval, \
      v_ints[0:n_v_ints], \
      v_doubles[0:n_v_doubles], \
      n_isotopes, \
      n_gridpoints, \
      lookups, \
      energy_grid[0:n_isotopes*n_gridpoints], \
      nuclide_grids[0:n_isotopes][0:n_gridpoints], \
      grid_ptrs[0:n_isotopes*n_isotopes*n_gridpoints], \
      mats[0:size_mats], \
      mats_idx[0:12], \
      concs[0:size_mats], \
      num_nucs[0:12], \
      dist[0:12], \
      rands[0:2*lookups] \
      )
#endif
  {

#ifndef _OPENACC
    // In OMP, initialize a private RNG seed for each thread
    thread = omp_get_thread_num();
    seed   = (thread+1)*19+17;
#endif

    // Dummy variable for bounds of XS lookup loop
    const int _lookups = lookups;

    // === The XS lookup loop
#ifndef _OPENACC
#pragma omp for schedule(dynamic)
#else
#pragma acc loop independent gang, vector(32) private(seed, mat) reduction(+:dval, vhash)
#endif
    for(int i = 0; i < _lookups; i++)
    {
#ifndef _OPENACC
      // Status text
      if( INFO && mype == 0 && thread == 0 && i % 1000 == 0 )
        printf("\rCalculating XS's... (%.0lf%% completed)",
            (i / ( (double)in.lookups / (double) in.nthreads ))
            / (double) in.nthreads * 100.0);
#endif

      // Randomly pick an energy and material for the particle
      double roll;
#ifndef _OPENACC
      // In OMP, generate random numbers on the fly
  #ifdef VERIFICATION
#pragma omp critical
      {
        p_energy = rn_v();
        roll = rn_v();
      }
  #else
      p_energy = rn(&seed);
      roll = rn(&seed);
  #endif
#else
      // In _OPENACC, use pre-generated random numbers
      p_energy = rands[2*i];
      roll = rands[2*i+1];
#endif

      // Use distribution to pick a material
      // ( inlined from pick_mat(mat_roll))
      for(mat = 0; mat < 12; mat++)
      {
        double running = 0;
        for(int j = mat; j > 0; j-- )
          running += dist[j];
        if( roll < running )
          break;
      }
      mat = mat % 12;

      // This returns the macro_xs_vector, but we're not going
      // to do anything with it in this program, so return value
      // is written over.
      // INLINE: calculate_macro_xs( p_energy, mat, n_isotopes,
      //     n_gridpoints, num_nucs, concs,
      //     energy_grid, grid_ptrs, nuclide_grids, mats, mats_idx,
      //     macro_xs_vector );
      double macro_xs_0 = 0;
      double macro_xs_1 = 0;
      double macro_xs_2 = 0;
      double macro_xs_3 = 0;
      double macro_xs_4 = 0;

      // binary search for energy on unionized energy grid (UEG)
      long idx = grid_search(n_isotopes * n_gridpoints, p_energy, energy_grid);	

      // Once we find the pointer array on the UEG, we can pull the data
      // from the respective nuclide grids, as well as the nuclide
      // concentration data for the material
      // Each nuclide from the material needs to have its micro-XS array
      // looked up & interpolatied (via calculate_micro_xs). Then, the
      // micro XS is multiplied by the concentration of that nuclide
      // in the material, and added to the total macro XS array.
#ifdef _OPENACC
// This inner loop is parallelizable. 
// However, we have found that performance is better if it is executed sequentialy.
#pragma acc loop seq reduction(+:macro_xs_0, macro_xs_1, macro_xs_2, macro_xs_3, macro_xs_4)
#endif
      for(int j = 0; j < num_nucs[mat]; j++)
      {
        
        // the nuclide we are looking up
        int p_nuc = mats[mats_idx[mat] + j];      
        // the concentration of the nuclide in the material
        double conc = concs[mats_idx[mat] + j];  
        // Interpolation factor
        double f;
        // Bounding energy gridpoints
        NuclideGridPoint * low, * high;

        // pull ptr from energy grid and check to ensure that
        // we're not reading off the end of the nuclide's grid
        if( grid_ptrs[energy_grid[idx].xs_ptrs + p_nuc] == n_gridpoints - 1 )
          low = &nuclide_grids[p_nuc][grid_ptrs[energy_grid[idx].xs_ptrs + p_nuc] - 1];
        else
          low = &nuclide_grids[p_nuc][grid_ptrs[energy_grid[idx].xs_ptrs + p_nuc]];

        high = low + 1;

        // calculate the re-useable interpolation factor
        f = (high->energy - p_energy) / (high->energy - low->energy);

        // Total XS
        macro_xs_0 += conc * (high->total_xs - f * (high->total_xs - low->total_xs));

        // Elastic XS
        macro_xs_1 += conc * (high->elastic_xs - f * (high->elastic_xs - low->elastic_xs));

        // Absorbtion XS
        macro_xs_2 += conc * (high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs));

        // Fission XS
        macro_xs_3 += conc * (high->fission_xs - f * (high->fission_xs - low->fission_xs));

        // Nu Fission XS
        macro_xs_4 += conc * (high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs));

      } // END: for( int j = 0; j < num_nucs[mat]; j++ )

      // Accumulate results into a dummy variable for reduction
      dval += (mat + p_energy + macro_xs_0 + macro_xs_1 + macro_xs_2 + macro_xs_3 + macro_xs_4);

      // Verification hash calculation
      // This method provides a consistent hash accross
      // architectures and compilers.
#ifdef VERIFICATION
  #ifndef _OPENACC
      // In OMP, hash results on-the-fly
      char line[256];
      sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
          p_energy, mat,
          macro_xs_0,
          macro_xs_1,
          macro_xs_2,
          macro_xs_3,
          macro_xs_4);
      unsigned long long vhash_local = hash(line, 10000);
#pragma omp atomic
      vhash += vhash_local;
  #else
      // In _OPENACC, results are stored and hashed on the host
      v_ints[i] = mat;
      v_doubles[6*i] = p_energy;
      v_doubles[6*i+1] = macro_xs_0;
      v_doubles[6*i+2] = macro_xs_1;
      v_doubles[6*i+3] = macro_xs_2;
      v_doubles[6*i+4] = macro_xs_3;
      v_doubles[6*i+5] = macro_xs_4;
  #endif
#endif
    } // END: for(int i = 0; i < _lookups; i++)
  } // END:  #pragma omp parallel OR #pragma acc kernels
#ifndef _OPENACC
  tock = omp_get_wtime();
#else
  tock = timer();
#pragma acc exit data \
  copyout(\
      vhash, \
      dval, \
      v_ints[0:n_v_ints], \
      v_doubles[0:n_v_doubles] \
      )
#endif

// For _OPENACC, hash the results
#ifdef _OPENACC
  #ifdef VERIFICATION
  for(int i = 0; i < lookups; i++){
    char line[256];
    sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
        v_doubles[6*i], v_ints[i], v_doubles[6*i+1], v_doubles[6*i+2], v_doubles[6*i+3], v_doubles[6*i+4], v_doubles[6*i+5]);
    unsigned long long vhash_local = hash(line, 10000);
    vhash += vhash_local;
  }
  #endif
#endif

  // Print / Save Results and Exit
  print_results(in, mype, tock-tick, nprocs, dval, vhash);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
