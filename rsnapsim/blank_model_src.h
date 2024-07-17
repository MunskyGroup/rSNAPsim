extern void generic_ssa_cpp(int* particle_array, int* state_array, int* resource_array, // arrays to fill
                    int* particle_x0, int* state_x0,  int* resource_x0, int* rxns, int* rxnids, double* elongs, int* probes,
                    double* parameters, int npars,       // initial simulation state
                    int max_rib, int n_constant_reactions, int n_ribosome_reactions,  // constants of the simulation
                    int n_colors, int n_states, int n_resources, int Nt, double* time_vector, double tf, double burnin, int seed, int L, int error_code);