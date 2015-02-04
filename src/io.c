#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

// Prints program logo
void logo(int version)
{
  border_print();
  printf(
      "        ____  _______________        _  _______ ____                  __   \n"  
      "       / __ \\/ ____/ ____/   |      | |/ / ___// __ )___  ____  _____/ /_  \n"
      "      / / / / /   / /   / /| |______|   /\\__ \\/ __  / _ \\/ __ \\/ ___/ __ \\ \n"
      "     / /_/ / /___/ /___/ ___ /_____/   |___/ / /_/ /  __/ / / / /__/ / / / \n"
      "     \\____/\\____/\\____/_/  |_|    /_/|_/____/_____/\\___/_/ /_/\\___/_/ /_/  \n\n"
      );
  border_print();
  center_print("Developed at Argonne National Laboratory and Rice University", 79);
  // center_print("Developed at:", 79);
  // center_print("Argonne National Laboratory", 79);
  // center_print("Rice University", 79);
  char v[100];
  sprintf(v, "Version: %d", version);
  center_print(v, 79);
  border_print();
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
  int length = strlen(s);
  int i;
  for (i=0; i<=(width-length)/2; i++) {
    fputs(" ", stdout);
  }
  fputs(s, stdout);
  fputs("\n", stdout);
}

void print_results( Inputs in, int mype, double runtime, int nprocs, double V_sum )
{
  // Calculate Lookups per sec
  int lookups_per_sec = (int) ((double) in.lookups / runtime);

  // If running in MPI, reduce timing statistics and calculate average
#ifdef MPI
  int total_lookups = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&lookups_per_sec, &total_lookups, 1, MPI_INT,
      MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  // Print output
  if( mype == 0 )
  {
    border_print();
    center_print("RESULTS", 79);
    border_print();

    // Print the results
    printf("Threads:     %d\n", in.n_threads);
#ifdef MPI
    printf("MPI ranks:   %d\n", nprocs);
#endif
#ifdef MPI
    printf("Total Lookups/s:            ");
    fancy_int(total_lookups);
    printf("Avg Lookups/s per MPI rank: ");
    fancy_int(total_lookups / nprocs);
#else
    printf("Runtime:     %.3lf seconds\n", runtime);
    printf("Lookups:     "); fancy_int(in.lookups);
    printf("Lookups/s:   ");
    fancy_int(lookups_per_sec);
#endif
#ifdef VERIFICATION
    printf("Verification sum: %0.5f\n", V_sum / in.lookups);
#endif
    border_print();

    // For bechmarking, output lookup/s data to file
    if( SAVE )
    {
      FILE * out = fopen( "results.txt", "a" );
      fprintf(out, "%d\t%d\n", in.n_threads, lookups_per_sec);
      fclose(out);
    }
  }
}

void print_inputs(Inputs in, int nprocs, int version )
{
  // Calculate Estimate of Memory Usage
  int mem_tot = estimate_mem_usage( in );
  logo(version);
  center_print("INPUT SUMMARY", 79);
  border_print();
#ifdef VERIFICATION
  printf("Verification Mode:            on\n");
#endif
  printf("Materials:                    %d\n", 12);
  printf("H-M Benchmark Size:           %s\n", in.HM);
  printf("Total Nuclides:               %ld\n", in.n_isotopes);
  printf("Gridpoints (per Nuclide):     ");
  fancy_int(in.n_gridpoints);
  printf("Unionized Energy Gridpoints:  ");
  fancy_int(in.n_isotopes*in.n_gridpoints);
  printf("XS Lookups:                   "); fancy_int(in.lookups);
#ifdef MPI
  printf("MPI Ranks:                    %d\n", nprocs);
  printf("OMP Threads per MPI Rank:     %d\n", in.n_threads);
  printf("Mem Usage per MPI Rank (MB):  "); fancy_int(mem_tot);
#else
  printf("Threads:                      %d\n", in.n_threads);
  printf("Est. Memory Usage (MB):       "); fancy_int(mem_tot);
#endif
  border_print();
  center_print("INITIALIZATION", 79);
  border_print();
}

void border_print(void)
{
  printf(
      "==================================================================="
      "=============\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int( long a )
{
  if( a < 1000 )
    printf("%ld\n",a);

  else if( a >= 1000 && a < 1000000 )
    printf("%ld,%03ld\n", a / 1000, a % 1000);

  else if( a >= 1000000 && a < 1000000000 )
    printf("%ld,%03ld,%03ld\n",a / 1000000,(a % 1000000) / 1000,a % 1000 );

  else if( a >= 1000000000 )
    printf("%ld,%03ld,%03ld,%03ld\n",
        a / 1000000000,
        (a % 1000000000) / 1000000,
        (a % 1000000) / 1000,
        a % 1000 );
  else
    printf("%ld\n",a);
}

void print_CLI_error(void)
{
  printf("Usage: ./XSBench <options>\n");
  printf("Options include:\n");
  printf("  -t <threads>     Number of OpenMP threads to run\n");
  printf("  -s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)\n");
  printf("  -g <gridpoints>  Number of gridpoints per nuclide (overrides -s defaults)\n");
  printf("  -l <lookups>     Number of Cross-section (XS) lookups\n");
  printf("Default is equivalent to: -s large -l 15000000\n");
  printf("See readme for full description of default run values\n");
  exit(4);
}

Inputs read_CLI( int argc, char *const argv[] )
{
  Inputs input;

  // defaults to max threads on the system
  if (getenv("OMP_NUM_THREADS") != NULL)
    input.n_threads = atoi(getenv("OMP_NUM_THREADS"));
  else
    input.n_threads = omp_get_num_procs();

  // defaults to 355 (corresponding to H-M Large benchmark)
  input.n_isotopes = 355;

  // defaults to 11303 (corresponding to H-M Large benchmark)
  input.n_gridpoints = 11303;

  // defaults to 15,000,000
  input.lookups = 15000000;

  // defaults to H-M Large benchmark
  input.HM = (char *) malloc(128 * sizeof(char));
  strcpy(input.HM, "large");

  // defaults to OpenMP mode
  input.mode = (char *) malloc(128 * sizeof(char));
  strcpy(input.mode, "mode = OpenMP");

  // defaults to hybrid kernel
  input.kernel = (char *) malloc(128 * sizeof(char));
  strcpy(input.kernel, "hybridLookupKernel.okl");

  // Check if user sets gridpoints
  bool user_g = false;

  // Get input
  int opt;
  while ((opt = getopt(argc, argv, ":t:l:g:o:i:s:m:k:")) != -1) {
    switch (opt) {
      case 't': input.n_threads = atoi(optarg); break;
      case 'l': input.lookups = atoi(optarg); break;
      case 'g': input.n_gridpoints = atol(optarg); user_g = 1; break;
      case 'o': input.outer_dim = atol(optarg); break;
      case 'i': input.inner_dim = atol(optarg); break;
      case 's': free(input.HM); input.HM = optarg; break;
      case 'm': free(input.mode); input.mode = optarg; break;
      case 'k': free(input.kernel); input.kernel = optarg; break;
      default:  print_CLI_error();
    }
  }

  // Validate numerical input
  if( (input.n_threads < 1) | (input.lookups < 1) | (input.n_gridpoints < 1) |
      (input.outer_dim < 1) | (input.inner_dim < 1))
    print_CLI_error();

  // Validate HM size
  if( strcasecmp(input.HM, "small") != 0 &&
      strcasecmp(input.HM, "large") != 0 &&
      strcasecmp(input.HM, "XL") != 0 &&
      strcasecmp(input.HM, "XXL") != 0 )
    print_CLI_error();

  // Set HM size specific parameters
  // (defaults to large)
  if( strcasecmp(input.HM, "small") == 0 )
    input.n_isotopes = 68;
  else if( strcasecmp(input.HM, "XL") == 0 && user_g == false )
    input.n_gridpoints = 238847; // sized to make 120 GB XS data
  else if( strcasecmp(input.HM, "XXL") == 0 && user_g == false )
    input.n_gridpoints = 238847 * 2.1; // 252 GB XS data

  // Return input struct
  return input;
}
