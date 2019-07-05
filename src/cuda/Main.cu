#include "XSbench_header.cuh"


int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	int version = 19;
	int mype = 0;
	double omp_start, omp_end;
	int nprocs = 1;
	unsigned long long verification;

	// rand() is only used in the serial initialization stages.
	// A custom RNG is used in parallel portions.
	srand(26);

	// Process CLI Fields -- store in "Inputs" structure
	Inputs in = read_CLI( argc, argv );

	// Print-out of Input Summary
	if( mype == 0 )
		print_inputs( in, nprocs, version );

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// This is not reflective of a real Monte Carlo simulation workload,
	// therefore, do not profile this region!
	// =====================================================================
	
	SimulationData GSD = grid_init_do_not_profile( in, mype );

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation
	// This is the section that should be profiled, as it reflects a 
	// realistic continuous energy Monte Carlo macroscopic cross section
	// lookup kernel.
	// =====================================================================

	if( mype == 0 )
	{
		printf("\n");
		border_print();
		center_print("SIMULATION", 79);
		border_print();
	}

	// Start Simulation Timer
	omp_start = get_time();

	// Run simulation
	if( in.simulation_method == EVENT_BASED )
		verification = run_event_based_simulation(in, GSD, mype);
	else
	{
		printf("History-based simulation not implemented in CUDA code. Use event-based method with \"-m event\" argument\n");
		exit(1);
	}

	if( mype == 0)	
	{	
		printf("\n" );
		printf("Simulation complete.\n" );
	}

	// End Simulation Timer
	omp_end = get_time();

	// =====================================================================
	// Output Results & Finalize
	// =====================================================================

	// Final Hash Step
	verification = verification % 1000000;

	// Print / Save Results and Exit
	print_results( in, mype, omp_end-omp_start, nprocs, verification );

	return 0;
}
