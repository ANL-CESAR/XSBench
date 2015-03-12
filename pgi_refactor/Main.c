#include "XSbench_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	int version = 13;
	int i, thread, mat;
	RNG_INT seed;
	double tick, tock, p_energy;
	VHASH_TYPE vval = 0;
	unsigned long long vhash = 0;
	int nprocs;
	double roll;

	char HM[6];

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

	double macro_xs_vector[5];

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
	Inputs in = read_CLI(argc, argv);
	const int nthreads = in.nthreads;
	const long n_isotopes = in.n_isotopes;
	const long n_gridpoints = in.n_gridpoints; //why not a const?
	const int lookups = in.lookups;

	// Print-out of Input Summary
	print_inputs(in, nprocs, version);

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	#ifndef BINARY_READ
	printf("Generating Nuclide Energy Grids...\n");
	#endif

	// allocate nuclide_grids[0:n_isotopes][0:n_gridpoints]
	NuclideGridPoint (* restrict nuclide_grids)[(long) n_gridpoints] = 
		(NuclideGridPoint (*)[(long) n_gridpoints]) 
		malloc(n_isotopes * n_gridpoints * sizeof(NuclideGridPoint));

	#ifdef VERIFICATION
	generate_grids_v(n_isotopes, n_gridpoints, nuclide_grids);
	#else
	generate_grids(n_isotopes, n_gridpoints, nuclide_grids);	
	#endif

	// Sort grids by energy
	#ifndef BINARY_READ
	printf("Sorting Nuclide Energy Grids...\n");
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
	printf("Reading data from \"XS_data.dat\" file...\n");
	binary_read(n_isotopes, n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
	#endif

	// Get material data
	printf("Loading Mats...\n");

	int size_mats;
	if (n_isotopes == 68) 
		size_mats = 197;
	else
		size_mats = 484;

	int * restrict num_nucs = load_num_nucs(n_isotopes);
	int * restrict mats_idx = load_mats_idx(num_nucs);
	int * restrict mats     = load_mats(num_nucs, mats_idx, size_mats, n_isotopes);

	#ifdef VERIFICATION
	double * restrict concs = load_concs_v(size_mats);
	#else
	double * restrict concs = load_concs(size_mats);
	#endif

	// Generate a stream of random numbers to copy in to device
	double * restrict rands = malloc(2*lookups*sizeof(double));
	for(i=0; i<lookups; i++){
		rands[2*i] = rn_v();
		rands[2*i+1] = rn_v();
	}

	// Create arrays to store values for verification
	int * restrict v_ints = malloc(lookups*sizeof(int));
	double * restrict v_doubles = malloc(6*lookups*sizeof(double));

	#ifdef BINARY_DUMP
	printf("Dumping data to binary file...\n");
	binary_dump(n_isotopes, n_gridpoints, nuclide_grids, energy_grid, grid_ptrs);
	printf("Binary file \"XS_data.dat\" written! Exiting...\n");
	return 0;
	#endif

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================

	printf("\n");
	border_print();
	center_print("SIMULATION", 79);
	border_print();

	tick = timer();

	#ifdef ACC
	#pragma acc data \
	copy(vhash, vval, v_ints[0:lookups], v_doubles[0:6*lookups]) \
	copyin( \
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
	     rands[0:2*lookups] ) \
	create(macro_xs_vector[0:5])
	#endif
	{
		#ifdef ACC
		#pragma acc kernels 
		#endif
		{	
			// XS Lookup Loop
			const int _lookups = lookups;
			#ifdef ACC
			#pragma acc loop independent gang private(seed, macro_xs_vector, mat) reduction(+:vval, vhash)
			#endif
			for(i = 0; i < _lookups; i++)
			{
				//seed = ((i % 10) +1)*19+17;

				// Randomly pick an energy and material for the particle
				//p_energy = rn(&seed);
				//roll = rn(&seed);
				p_energy = rands[2*i];
				roll = rands[2*i+1];

				// INLINE:  pick_mat(mat_roll)
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
				double xs_vector[5];
				int p_nuc; // the nuclide we are looking up
				long idx = 0;	
				double conc; // the concentration of the nuclide in the material

				// cleans out macro_xs_vector
				for( int k = 0; k < 5; k++ )
					macro_xs_vector[k] = 0;

				// binary search for energy on unionized energy grid (UEG)
				idx = grid_search(n_isotopes * n_gridpoints, p_energy, energy_grid);	

				// Once we find the pointer array on the UEG, we can pull the data
				// from the respective nuclide grids, as well as the nuclide
				// concentration data for the material
				// Each nuclide from the material needs to have its micro-XS array
				// looked up & interpolatied (via calculate_micro_xs). Then, the
				// micro XS is multiplied by the concentration of that nuclide
				// in the material, and added to the total macro XS array.
				#ifdef ACC
				#pragma acc loop worker private(xs_vector, p_nuc, conc)
				#endif
				for(int j = 0; j < num_nucs[mat]; j++)
				{

					p_nuc = mats[mats_idx[mat] + j];
					conc = concs[mats_idx[mat] + j];

					// INLINE: calculate_micro_xs( p_energy, p_nuc, n_isotopes,
					//     n_gridpoints, energy_grid, grid_ptrs,
					//     nuclide_grids, idx, xs_vector );
					double f;
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
					xs_vector[0] = high->total_xs - f * (high->total_xs - low->total_xs);

					// Elastic XS
					xs_vector[1] = high->elastic_xs - f * (high->elastic_xs - low->elastic_xs);

					// Absorbtion XS
					xs_vector[2] = high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs);

					// Fission XS
					xs_vector[3] = high->fission_xs - f * (high->fission_xs - low->fission_xs);

					// Nu Fission XS
					xs_vector[4] = high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs);

					for(int k = 0; k < 5; k++)
						macro_xs_vector[k] += xs_vector[k] * conc;

				} // END: for( int j = 0; j < num_nucs[mat]; j++ )

				// Verification hash calculation
				// This method provides a consistent hash accross
				// architectures and compilers.
				vval += (mat + p_energy + macro_xs_vector[0] + macro_xs_vector[1] + macro_xs_vector[2] + macro_xs_vector[3] + macro_xs_vector[4]);
				#ifdef VERIFICATION
				v_ints[i] = mat;
				v_doubles[6*i] = p_energy;
				for(int k = 0; k < 5; k++)
					v_doubles[6*i+k+1] = macro_xs_vector[k];
				#endif
			} // END: for( i = 0; i < _lookups; i++ )
		} // END: #pragma acc parallel OR #pragma omp parallel
	} // END: #pragma acc data

	tock = timer();

	#ifdef VERIFICATION
	FILE *fp = fopen("out", "w");
	for(int i = 0; i < lookups; i++){
		fprintf(fp, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf\n",
		     v_doubles[6*i], v_ints[i], v_doubles[6*i+1], v_doubles[6*i+2], v_doubles[6*i+3], v_doubles[6*i+4], v_doubles[6*i+5]);
		char line[256];
		sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
		     v_doubles[6*i], v_ints[i], v_doubles[6*i+1], v_doubles[6*i+2], v_doubles[6*i+3], v_doubles[6*i+4], v_doubles[6*i+5]);
		unsigned long long vhash_local = hash(line, 10000);
		vhash += vhash_local;
	}
	#endif

	// Print / Save Results and Exit
	print_results(in, 0, tock-tick, nprocs, vval, vhash);

	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
