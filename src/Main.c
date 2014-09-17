#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main(int argc, char* argv[])
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	Inputs in;
	int version = 13;
	int mype = 0;
	int max_procs = omp_get_num_procs();
	int i, thread, mat;
	unsigned long seed;
	double omp_start, omp_end, p_energy;
	unsigned long long vhash = 0;
	int nprocs;
	//double elapsed_time;

	//Inputs in;
	double *** nuclide_grids;
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
	in = read_CLI(argc, argv);

	//!!!!!
	in.n_isotopes = 68;

	// Set number of OpenMP Threads
	omp_set_num_threads(in.nthreads); 

	// Print-out of Input Summary
	if(mype == 0) print_inputs(in, nprocs, version);

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	#ifndef BINARY_READ
	if(mype == 0) printf("Generating Nuclide Energy Grids...\n");
	#endif
	
	nuclide_grids = gpmatrix(in.n_isotopes,in.n_gridpoints);

	#ifdef VERIFICATION
	generate_grids_v(nuclide_grids, in.n_isotopes, in.n_gridpoints);	
	#else
	generate_grids(nuclide_grids, in.n_isotopes, in.n_gridpoints);	
	#endif

	// Sort grids by energy
	#ifndef BINARY_READ
	if(mype == 0) printf("Sorting Nuclide Energy Grids...\n");
	sort_nuclide_grids(nuclide_grids, in.n_isotopes, in.n_gridpoints);
	#endif

	// Prepare Unionized Energy Grid Framework
	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	#ifndef BINARY_READ
	energy_grid = generate_energy_grid(in.n_isotopes, in.n_gridpoints, nuclide_grids);
	grid_ptrs = generate_grid_ptrs(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid);	
	#else
	energy_grid = malloc(in.n_isotopes*in.n_gridpoints*sizeof(double));
	grid_ptrs = (int *) malloc(in.n_isotopes*in.n_gridpoints*in.n_isotopes*sizeof(int));
	#endif

	#ifdef BINARY_READ
	if(mype == 0) printf("Reading data from \"XS_data.dat\" file...\n");
	binary_read(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
	#endif
	
	// Get material data
	if(mype == 0) printf("Loading Mats...\n");
	if(in.n_isotopes == 68) size_mats = 197;
	else size_mats = 484;
	num_nucs  = load_num_nucs(in.n_isotopes);
	mats_ptr  = load_mats_ptr(num_nucs);
	mats      = load_mats(num_nucs, mats_ptr, size_mats, in.n_isotopes);

	#ifdef VERIFICATION
	concs = load_concs_v(size_mats);
	#else
	concs = load_concs(size_mats);
	#endif

	#ifdef BINARY_DUMP
	if(mype == 0) printf("Dumping data to binary file...\n");
	binary_dump(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
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
		in.nthreads = bench_n;
		omp_set_num_threads(in.nthreads);
 	#endif

	if(mype == 0)
	{
		printf("\n");
		border_print();
		center_print("SIMULATION", 79);
		border_print();
	}

	omp_start = omp_get_wtime();

	#pragma omp parallel default(none) \
	private(i, thread, p_energy, mat, seed, vhash_local, line, macro_xs_vector) \
	shared( max_procs, in, energy_grid, nuclide_grids, grid_ptrs, \
	        mats_ptr, mats, concs, num_nucs, mype, vhash) 

	//elapsed_time = timer();
//	#pragma acc data \
	copy(in, vhash) \
	copyin(num_nucs[0:in.n_isotopes], concs[0:size_mats], mats[0:size_mats], mats_ptr[0:12], \
	       energy_grid[0:in.n_isotopes*in.n_gridpoints]) \
	       /*grid_ptrs[0:in.n_isotopes*in.n_isotopes*in.n_gridpoints], \
	       nuclide_grids[0:in.n_isotopes*in.n_gridpoints*6])*/ \
	create(p_energy, mat, macro_xs_vector, vhash_local, line, seed)
	{	

		thread = omp_get_thread_num();
		seed   = (thread+1)*19+17;
		//seed = 13; //what to do for openacc?

		// XS Lookup Loop
		#pragma omp for schedule(dynamic)
		//#pragma acc parallel \
		private(macro_xs_vector, p_energy, mat, seed, vhash_local, line)
		for(i=0; i<in.lookups; i++)
		{
			// Status text
			if( INFO && mype == 0 && thread == 0 && i % 1000 == 0 )
				printf("\rCalculating XS's... (%.0lf%% completed)",
						(i / ( (double)in.lookups / (double) in.nthreads ))
						/ (double) in.nthreads * 100.0);

			// Randomly pick an energy and material for the particle
			#ifdef VERIFICATION
			#pragma omp critical
			{
				mat = pick_mat(&seed); 
				p_energy = rn_v();
			}
			#else
			mat = pick_mat(&seed); 
			p_energy = rn(&seed);
			#endif
		
			// This returns the macro_xs_vector, but we're not going
			// to do anything with it in this program, so return value
			// is written over.
			calculate_macro_xs(p_energy, mat, in.n_isotopes, in.n_gridpoints,
					   num_nucs, concs, energy_grid, nuclide_grids,
					   grid_ptrs, mats, mats_ptr, macro_xs_vector);

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
			vhash_local = hash(line, 10000);
			#pragma omp atomic
			vhash += vhash_local;
			#endif
		}
	}

	omp_end = omp_get_wtime();

	if(mype == 0) printf("\nSimulation complete.\n" );

	//elapsed_time = timer() - elapsed_time;
	//printf("Acclerator elapsed time: %lf seconds\n", elapsed_time);

	// Print / Save Results and Exit
	print_results(in, mype, omp_end-omp_start, nprocs, vhash);

	#ifdef BENCHMARK
	}
	#endif

	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
