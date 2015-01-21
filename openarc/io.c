#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

// Prints program logo
void logo(int version)
{
	char v[100];
	border_print();
	printf(
	"                   __   __ ___________                 _                        \n"
	"                   \\ \\ / //  ___| ___ \\               | |                       \n"
	"                    \\ V / \\ `--.| |_/ / ___ _ __   ___| |__                     \n"
	"                    /   \\  `--. \\ ___ \\/ _ \\ '_ \\ / __| '_ \\                    \n"
	"                   / /^\\ \\/\\__/ / |_/ /  __/ | | | (__| | | |                   \n"
	"                   \\/   \\/\\____/\\____/ \\___|_| |_|\\___|_| |_|                   \n\n"
    );
	border_print();
	center_print("Developed at Argonne National Laboratory", 79);
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

void print_results( int nthreads, long n_isotopes, long n_gridpoints, int
    lookups, char *HM,  int mype, double runtime, int nprocs, unsigned long
    long vhash )
{
	// Calculate Lookups per sec
	int lookups_per_sec = (int) ((double) lookups / runtime);
	FILE * out;	

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
		printf("Threads:     %d\n", nthreads);
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
		printf("Lookups:     "); fancy_int(lookups);
		printf("Lookups/s:   ");
		fancy_int(lookups_per_sec);
		#endif
		#ifdef VERIFICATION
		printf("Verification checksum: %llu\n", vhash);
		#endif
		border_print();

		// For bechmarking, output lookup/s data to file
		if( SAVE )
		{
			out = fopen( "results.txt", "a" );
			fprintf(out, "%d\t%d\n", nthreads, lookups_per_sec);
			fclose(out);
		}
	}
}

void print_inputs(int nthreads, long n_isotopes, 
    long n_gridpoints, int lookups, char *HM, int nprocs, int version )
{
	// Calculate Estimate of Memory Usage
	int mem_tot = estimate_mem_usage( n_isotopes, n_gridpoints );
	logo(version);
	center_print("INPUT SUMMARY", 79);
	border_print();
	#ifdef VERIFICATION
	printf("Verification Mode:            on\n");
	#endif
	printf("Materials:                    %d\n", 12);
	printf("H-M Benchmark Size:           %s\n", HM);
	printf("Total Nuclides:               %ld\n", n_isotopes);
	printf("Gridpoints (per Nuclide):     ");
	fancy_int(n_gridpoints);
	printf("Unionized Energy Gridpoints:  ");
	fancy_int(n_isotopes*n_gridpoints);
	printf("XS Lookups:                   "); fancy_int(lookups);
	#ifdef MPI
	printf("MPI Ranks:                    %d\n", nprocs);
	printf("OMP Threads per MPI Rank:     %d\n", nthreads);
	printf("Mem Usage per MPI Rank (MB):  "); fancy_int(mem_tot);
	#else
	printf("Threads:                      %d\n", nthreads);
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

void read_CLI( int argc, char * argv[], int *nthreads, long *n_isotopes, long
    *n_gridpoints, int *lookups, char *HM )
{
	int user_g, i;
	char * arg;	

	// defaults to max threads on the system	
	*nthreads = omp_get_num_procs();
	
	// defaults to 355 (corresponding to H-M Large benchmark)
	*n_isotopes = 355;
	
	// defaults to 11303 (corresponding to H-M Large benchmark)
	*n_gridpoints = 11303;
	
	// defaults to 15,000,000
	*lookups = 15000000;
	
	// defaults to H-M Large benchmark
  strcpy(HM, "large");
	
	// Check if user sets these
	user_g = 0;
	
	// Collect Raw Input
	for( i = 1; i < argc; i++ )
	{
		arg = argv[i];

		// nthreads (-t)
		if( strcmp(arg, "-t") == 0 )
		{
			if( ++i < argc )
				*nthreads = atoi(argv[i]);
			else
				print_CLI_error();
		}
		// n_gridpoints (-g)
		else if( strcmp(arg, "-g") == 0 )
		{	
			if( ++i < argc )
			{
				user_g = 1;
				*n_gridpoints = atol(argv[i]);
			}
			else
				print_CLI_error();
		}
		// lookups (-l)
		else if( strcmp(arg, "-l") == 0 )
		{
			if( ++i < argc )
				*lookups = atoi(argv[i]);
			else
				print_CLI_error();
		}
		// HM (-s)
		else if( strcmp(arg, "-s") == 0 )
		{	
			if( ++i < argc )
				strcpy(HM, argv[i]);
			else
				print_CLI_error();
		}
		else
			print_CLI_error();
	}

	// Validate Input

	// Validate nthreads
	if( *nthreads < 1 )
		print_CLI_error();
	
	// Validate n_isotopes
	if( *n_isotopes < 1 )
		print_CLI_error();
	
	// Validate n_gridpoints
	if( *n_gridpoints < 1 )
		print_CLI_error();

	// Validate lookups
	if( *lookups < 1 )
		print_CLI_error();
	
	// Validate HM size
	if( strcasecmp(HM, "small") != 0 &&
		strcasecmp(HM, "large") != 0 &&
		strcasecmp(HM, "XL") != 0 &&
		strcasecmp(HM, "XXL") != 0 )
		print_CLI_error();
	
	// Set HM size specific parameters
	// (defaults to large)
	if( strcasecmp(HM, "small") == 0 )
		*n_isotopes = 68;
	else if( strcasecmp(HM, "XL") == 0 && user_g == 0 )
		*n_gridpoints = 238847; // sized to make 120 GB XS data
	else if( strcasecmp(HM, "XXL") == 0 && user_g == 0 )
		*n_gridpoints = 238847 * 2.1; // 252 GB XS data

	return;
}
