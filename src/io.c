#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

void print_results( Inputs in, int mype, double runtime, int nprocs,
	unsigned long long vhash )
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
		printf("Threads:     %d\n", in.nthreads);
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
		printf("Verification checksum: %llu\n", vhash);
		#endif
		border_print();

		// For bechmarking, output lookup/s data to file
		if( SAVE )
		{
			FILE * out = fopen( "results.txt", "a" );
			fprintf(out, "%d\t%d\n", in.nthreads, lookups_per_sec);
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
	printf("OMP Threads per MPI Rank:     %d\n", nthreads);
	printf("Mem Usage per MPI Rank (MB):  "); fancy_int(mem_tot);
	#else
	printf("Threads:                      %d\n", in.nthreads);
	printf("Est. Memory Usage (MB):       "); fancy_int(mem_tot);
	#endif
	border_print();
	center_print("INITIALIZATION", 79);
	border_print();
}
