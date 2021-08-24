// Separate allocation and copy element of NuclideGridPoint
static aml_final_mapper_decl(NuclideGridPoint_mapper,
														 AML_MAPPER_FLAG_COPY | AML_MAPPER_FLAG_SPLIT,
                             NuclideGridPoint);

// Separate allocation but no copy element of int arrays.
static aml_final_mapper_decl(int_alloc_mapper, AML_MAPPER_FLAG_SPLIT, int);

// Bulk allocation and copy element of double arrays.
static aml_final_mapper_decl(double_copy_mapper, AML_MAPPER_FLAG_COPY, double);
// Bulk allocation and but no copy of double arrays.
static aml_final_mapper_decl(double_alloc_mapper, 0, double);
// Separate allocation but no copy element of unsigned long arrays.
static aml_final_mapper_decl(ulong_alloc_mapper, AML_MAPPER_FLAG_SPLIT, unsigned long);

// Top level structure.
// Copy but do not allocate.
// The pointer with space will be provided. Only child fields will require
// allocation.
aml_mapper_decl(SimulationData_mapper,
                AML_MAPPER_FLAG_COPY | AML_MAPPER_FLAG_SHALLOW, SimulationData,
                num_nucs, length_num_nucs, &int_alloc_mapper,
								concs, length_concs, &double_copy_mapper,
								mats, length_mats, &int_alloc_mapper,
								unionized_energy_array, length_unionized_energy_array, &double_copy_mapper,
								index_grid, length_index_grid, &int_alloc_mapper,
								nuclide_grid, length_nuclide_grid, &NuclideGridPoint_mapper,
								p_energy_samples, length_p_energy_samples, &double_alloc_mapper,
								mat_samples, length_mat_samples, &int_alloc_mapper);
