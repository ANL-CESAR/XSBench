#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main(int argc, char* argv[])
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
  double acc_start, acc_end;

  //Inputs
  int nthreads;
  long n_isotopes;
  long n_gridpoints;
  int lookups;
  char HM[6];

  double *nuclide_grids;
  double *energy_grid;
  int *grid_ptrs;
  int *index_data;
  int size_mats, *num_nucs, *mats_ptr, *mats;
  double *concs;
  int bench_n; // benchmark loop index
  double macro_xs_vector[5];
  char line[256]; // verification hash
  unsigned long long vhash_local; // verification hash

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
  read_CLI(argc, argv, &nthreads, &n_isotopes, &n_gridpoints, &lookups, HM);

  // Set number of OpenMP Threads
  omp_set_num_threads(nthreads); 

  // Print-out of Input Summary
  if(mype == 0) print_inputs(nthreads, n_isotopes, n_gridpoints, lookups, HM, nprocs, version);

  // =====================================================================
  // Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
  // =====================================================================

  // Allocate & fill energy grids
#ifndef BINARY_READ
  if(mype == 0) printf("Generating Nuclide Energy Grids...\n");
#endif

  nuclide_grids = (double *) malloc(n_isotopes *n_gridpoints * 6 * sizeof(double));

#ifdef VERIFICATION
  generate_grids_v(nuclide_grids,n_isotopes,n_gridpoints);	
#else
  generate_grids(nuclide_grids,n_isotopes,n_gridpoints);	
#endif

  // Sort grids by energy
#ifndef BINARY_READ
  if(mype == 0) printf("Sorting Nuclide Energy Grids...\n");
  sort_nuclide_grids(nuclide_grids,n_isotopes,n_gridpoints);
#endif

  // Prepare Unionized Energy Grid Framework
  // Double Indexing. Filling in energy_grid with pointers to the
  // nuclide_energy_grids.
#ifndef BINARY_READ
  energy_grid = generate_energy_grid(n_isotopes,n_gridpoints, nuclide_grids);
  grid_ptrs = generate_grid_ptrs(n_isotopes,n_gridpoints, nuclide_grids, energy_grid);	
#else
  energy_grid = malloc(n_isotopes*n_gridpoints*sizeof(double));
  grid_ptrs = (int *) malloc(n_isotopes*n_gridpoints*n_isotopes*sizeof(int));
#endif

#ifdef BINARY_READ
  if(mype == 0) printf("Reading data from \"XS_data.dat\" file...\n");
  binary_read(n_isotopes,n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
#endif

  // Get material data
  if(mype == 0) printf("Loading Mats...\n");
  if(n_isotopes == 68) size_mats = 197;
  else size_mats = 484;
  num_nucs  = load_num_nucs(n_isotopes);
  mats_ptr  = load_mats_ptr(num_nucs);
  mats      = load_mats(num_nucs, mats_ptr, size_mats,n_isotopes);

#ifdef VERIFICATION
  concs = load_concs_v(size_mats);
#else
  concs = load_concs(size_mats);
#endif

#ifdef BINARY_DUMP
  if(mype == 0) printf("Dumping data to binary file...\n");
  binary_dump(n_isotopes,n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
  if(mype == 0) printf("Binary file \"XS_data.dat\" written! Exiting...\n");
  return 0;
#endif

  // =====================================================================
  // Cross Section (XS) Parallel Lookup Simulation Begins
  // =====================================================================

  // Outer benchmark loop can loop through all possible # of threads
#ifdef BENCHMARK
  for(bench_n = 1; bench_n <=omp_get_num_procs(); bench_n++)
  {
    nthreads = bench_n;
    omp_set_num_threads(nthreads);
#endif

    if(mype == 0)
    {
      printf("\n");
      border_print();
      center_print("SIMULATION", 79);
      border_print();
    }

#ifndef OPENACC
    omp_start = omp_get_wtime();
#else
    acc_start = timer();
#endif

    //initialize papi with one thread (master) here
#ifdef PAPI
    if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
      fprintf(stderr, "PAPI library init error!\n");
      exit(1);
    }
#endif	

#ifndef OPENACC
#pragma omp parallel default(none) \
    private(i, thread, p_energy, mat, seed, vhash_local, line, macro_xs_vector) \
    shared( max_procs, nthreads, n_isotopes, n_gridpoints, lookups, HM, energy_grid, \
        nuclide_grids, grid_ptrs, mats_ptr, mats, concs, num_nucs, mype, vhash) 
#else
#pragma acc data \
    copy(vhash) \
    copyin(lookups, n_isotopes, n_gridpoints, \
        num_nucs[0:n_isotopes], concs[0:size_mats], mats[0:size_mats], mats_ptr[0:12], \
        energy_grid[0:n_isotopes*n_gridpoints], \
        grid_ptrs[0:n_isotopes*n_isotopes*n_gridpoints], \
        nuclide_grids[0:n_isotopes*n_gridpoints*6])
#endif
    {
      // Initialize parallel PAPI counters
#ifdef PAPI
      int eventset = PAPI_NULL; 
      int num_papi_events;
#pragma omp critical
      {
        counter_init(&eventset, &num_papi_events);
      }
#endif

#ifndef OPENACC
      thread = omp_get_thread_num();
      seed   = (thread+1)*19+17;
#else
      seed = 13; //what to do for openacc?
#endif

      // XS Lookup Loop
#ifndef OPENACC
#pragma omp for schedule(dynamic)
#else
//#pragma acc parallel for independent firstprivate(seed)  private(macro_xs_vector, p_energy, mat, vhash_local, line)
#pragma acc parallel for
#endif
      for(i=0; i<lookups; i++)
      {
#ifndef OPENACC
        // Status text
        if( INFO && mype == 0 && thread == 0 && i % 1000 == 0 )
          printf("\rCalculating XS's... (%.0lf%% completed)",
              (i / ( (double)lookups / (double)nthreads ))
              / (double)nthreads * 100.0);
#endif

        // Randomly pick an energy and material for the particle
#ifdef VERIFICATION
#ifndef OPENACC
#pragma omp critical
#endif
        {
          mat = pick_mat(&seed); 
          p_energy = rn_v();
        }
#else
        // INLINE: mat = pick_mat(&seed); 
        {

          double roll, running;
          int i, j;
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

          //double roll = (double) rand() / (double) RAND_MAX;
#ifdef VERIFICATION
          roll = rn_v();
#else
          roll = rn(&seed);
#endif
          // makes a pick based on the distro
          for( mat = 0; mat < 12; mat++ )
          {
            running = 0;
            for( j = mat; j > 0; j-- )
              running += dist[j];
            if( roll < running )
              break;
          }
          mat = mat % 12;
        } // END INLINE: mat = pick_mat(&seed); 

        p_energy = rn(&seed);
#endif

        // This returns the macro_xs_vector, but we're not going
        // to do anything with it in this program, so return value
        // is written over.
        // INLINE: calculate_macro_xs(p_energy, mat, n_isotopes,
        //   n_gridpoints, num_nucs, concs, energy_grid, nuclide_grids,
        //   grid_ptrs, mats, mats_ptr, macro_xs_vector);
        {
          double xs_vector[5];
          int p_nuc; // the nuclide we are looking up
          long idx = 0;	
          double conc; // the concentration of the nuclide in the material
          int j, k;

          // cleans out macro_xs_vector
          for(k=0; k<5; k++)
            macro_xs_vector[k] = 0;

          // binary search for energy on unionized energy grid (UEG)
          // INLINE: idx = grid_search(n_isotopes*n_gridpoints, p_energy, energy_grid);	
          {
            long lowerLimit = 0;
            long upperLimit = n_isotopes*n_gridpoints-1;
            long examinationPoint;
            long length = upperLimit - lowerLimit;

            while(length > 1){
              examinationPoint = lowerLimit + (length/2);	
              if(energy_grid[examinationPoint] > p_energy) upperLimit = examinationPoint;
              else lowerLimit = examinationPoint;
              length = upperLimit - lowerLimit;
            }
            idx = lowerLimit;
          } // END INLINE: idx = grid_search(...)

          // Once we find the pointer array on the UEG, we can pull the data
          // from the respective nuclide grids, as well as the nuclide
          // concentration data for the material
          // Each nuclide from the material needs to have its micro-XS array
          // looked up & interpolatied (via calculate_micro_xs). Then, the
          // micro XS is multiplied by the concentration of that nuclide
          // in the material, and added to the total macro XS array.
          for(j=0; j<num_nucs[mat]; j++){
            p_nuc = mats[mats_ptr[mat] + j];
            conc = concs[mats_ptr[mat] + j];

            // INLINE: calculate_micro_xs(p_energy, p_nuc, n_isotopes,
            //   n_gridpoints, energy_grid, nuclide_grids, grid_ptrs, idx,
            //   xs_vector);
            {	
              // Variables
              double f;
              double *low, *high;

              // pull ptr from energy grid and check to ensure that
              // we're not reading off the end of the nuclide's grid
              if(grid_ptrs[n_isotopes*idx + p_nuc] == n_gridpoints - 1){
                low = &nuclide_grids[p_nuc*n_gridpoints*6 + (grid_ptrs[n_isotopes*idx + p_nuc] - 1)*6];
                high = &nuclide_grids[p_nuc*n_gridpoints*6 + (grid_ptrs[n_isotopes*idx + p_nuc])*6];
              }
              else{
                low = &nuclide_grids[p_nuc*n_gridpoints*6 + (grid_ptrs[n_isotopes*idx + p_nuc])*6];
                high = &nuclide_grids[p_nuc*n_gridpoints*6 + (grid_ptrs[n_isotopes*idx + p_nuc] + 1)*6];
              }	

              // calculate the re-useable interpolation factor
              f = (high[0] - p_energy) / (high[0] - low[0]);

              // Total XS
              xs_vector[0] = high[1] - f * (high[1] - low[1]);

              // Elastic XS
              xs_vector[1] = high[2] - f * (high[2] - low[2]);

              // Absorbtion XS
              xs_vector[2] = high[3] - f * (high[3] - low[3]);

              // Fission XS
              xs_vector[3] = high[4] - f * (high[4] - low[4]);

              // Nu Fission XS
              xs_vector[4] = high[5] - f * (high[5] - low[5]);

            } // END INLINE: calculate_micro_xs(...)

            for(k=0; k<5; k++)
              macro_xs_vector[k] += xs_vector[k] * conc;

          } // END: for(j=0; j<num_nucs[mat]; j++)

        } // INLINE END: calculate_macro_xs(...)

        // Verification hash calculation
        // This method provides a consistent hash accross
        // architectures and compilers.
#ifdef VERIFICATION
        sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
            p_energy, mat,
            macro_xs_vector[0],
            macro_xs_vector[1],
            macro_xs_vector[2],
            macro_xs_vector[3],
            macro_xs_vector[4]);
        vhash_local = hash((unsigned char *)line, 10000);
#ifndef OPENACC
#pragma omp atomic
#endif
        vhash += vhash_local;
#endif
      } // END: for(i=0; i<lookups; i++)

      // Prints out thread local PAPI counters
#ifdef PAPI
      if( mype == 0 && thread == 0 )
      {
        printf("\n");
        border_print();
        center_print("PAPI COUNTER RESULTS", 79);
        border_print();
        printf("Count          \tSmybol      \tDescription\n");
      }
      {
#pragma omp barrier
      }
      counter_stop(&eventset, num_papi_events);
#endif
    }

#ifndef PAPI
    if( mype == 0) printf("\nSimulation complete.\n" );
#endif

#ifndef OPENACC
    omp_end = omp_get_wtime();
    print_results(nthreads, n_isotopes, n_gridpoints, lookups, HM, mype, omp_end-omp_start, nprocs, vhash);
#else
    acc_end = timer();
    print_results(nthreads, n_isotopes, n_gridpoints, lookups, HM, mype, acc_end-acc_start, nprocs, vhash);
#endif

#ifdef BENCHMARK
  }
#endif

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
