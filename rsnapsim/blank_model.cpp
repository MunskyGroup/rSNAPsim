#include <iostream>
#include <ctime>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <cstdlib>
#include <random>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;


void decompress_X_locations(const Eigen::MatrixXi& X_compressed,  Eigen::MatrixXi& X_lat, int max_L, int N_particles ) {
    // Decompress the X vector of notation [x1, x2, x3,...0] to binary vector X_full of notation 1xL [0, 1,0,0,1...L]

    X_lat.setZero(); // reset x_full binary vector 
    for (int k = 0; k < N_particles; k++){
        X_lat(0, X_compressed(0, k)-1 ) = 1;// update based on the x_compressed vector
    }
}


void compress_X_locations(Eigen::MatrixXi& X_compressed, const Eigen::MatrixXi& X_full, int max_L, int max_particles  ) {
    // Compress the x_full vector into the format [x1, x2, x3..... 0] for N ribs and front leading ribosomes
    X_compressed.setZero();
    int N = 0; 
    for (int i = max_L; i>= 0; i--){
        if (X_full(0,i) == 1){
        X_compressed(0,N ) = i + 1;
        N++;
        }
    }
}

/*
void custom_record( Eigen::MatrixXd X_compressed,  Eigen::MatrixXd X_states, double* parameters,   double t   ){

    int a = 1;

}
*/


void update_spatial_X( const Eigen::MatrixXi& X_compressed, Eigen::MatrixXi&  X_spatial, int NR, int frame_L ){
		X_spatial.setZero();
       
		for (int i = 0; i < NR; i++){
			X_spatial(0,i) = (X_compressed(0,i) + (X_compressed(0,i)>frame_L) + (X_compressed(0,i) > (2*frame_L-1)) )%frame_L; 
        }
}



Eigen::MatrixXd stepping_forward_propensities( const Eigen::MatrixXi& X_compressed, const Eigen::MatrixXi&  X_spatial, int n_ribosomes, const Eigen::MatrixXd& forward_rates, int R) {
    int loc_not_free; 
    int ribs_ahead;
    int less_than_equal_zero;
    int less_than_R; 
    Eigen::VectorXd wn_stepping;
    wn_stepping.setZero(n_ribosomes);

    for (int i = 0; i< n_ribosomes; i++){
        // how many ribosomes ahead of current ribosome?
        // if this is zero automatically fill the stepping reaction


        ribs_ahead = (X_spatial.block(0,0,1,n_ribosomes).array() - X_spatial(0,i) >  0  ).cast<int>().sum() ; 
        if (ribs_ahead == 0){
            loc_not_free = 0;
            wn_stepping(i) = forward_rates(0,X_compressed(0,i)-1);

        } 

        else{ // otherwise see if its free to step forward
            
            // how many are less than zero 
            less_than_equal_zero = (X_spatial.block(0,0,1,n_ribosomes).array()-X_spatial(0,i) <= 0).cast<int>().sum();
            // how many are less than ribosomal exclusion 
            less_than_R = (X_spatial.block(0,0,1,n_ribosomes).array()-X_spatial(0,i) < R).cast<int>().sum();

         
            if (less_than_R > less_than_equal_zero){ // if theres a ribosome thats positive and < R (ie one is excluding stepping forward)
                loc_not_free = 1;
            }
            else{
                loc_not_free = 0;
                wn_stepping(i) = forward_rates(0,X_compressed(0,i)-1 );
            }

        }
    }
    
    return wn_stepping;
}

Eigen::VectorXd propensity_function( const Eigen::MatrixXi& X_compressed,  const Eigen::MatrixXi& X_states, const Eigen::MatrixXi& X_spatial, const Eigen::MatrixXi& X_full,
                                      double* parameters, double* forward_rates,  double t, int n_total_reactions, int n_ribosomes, int R, int max_length, double& probe_binding_rate   ){

                                          
    Eigen::VectorXd wn;
    wn.setZero(n_total_reactions);

    Eigen::MatrixXd forward_rate_matrix(1,max_length+1);	// make forward_rates
    forward_rate_matrix.setZero();

	for(int i=0; i < max_length; i++)	{
		forward_rate_matrix(0,i) = forward_rates[i];
	}   
  
    //INSERT_GENERATED_PROPENSITY_HERE

    Eigen::VectorXd wn_full;
    wn_full.setZero(n_total_reactions + n_ribosomes);

    if (n_ribosomes > 0) {  
        Eigen::VectorXd wn_stepping = stepping_forward_propensities(X_compressed, X_spatial, n_ribosomes, forward_rate_matrix, R );
        wn_full << wn, wn_stepping;
    }
    else{
        wn_full << wn;
    }

    return wn_full;

}


void update_probe_matrix_leave(Eigen::MatrixXi& X_probes, int row_to_delete, int N_particles, int n_total_probes){
    // shift all rows up
    X_probes.block(row_to_delete, 0, N_particles - row_to_delete , n_total_probes ) << X_probes.block(row_to_delete+1, 0, N_particles - row_to_delete , n_total_probes );
    //set last row to zero
    X_probes.block(N_particles, 0, 1 , n_total_probes ).setZero();
}

// update the probe_matrix
void update_probe_matrix_step(Eigen::MatrixXi& X_probes, const Eigen::MatrixXi& X_compressed, int* probe_locations, double probe_binding_rate, int n_total_probes, int rib_stepped_too, int who_stepped, double randnum){
    
    for (int i = 0; i < n_total_probes; i ++){
        if (X_compressed(0,who_stepped)-1 == probe_locations[i] ){

            if (randnum < probe_binding_rate){
                X_probes(who_stepped, i) += 1;
            }

        }
    }
}

// get the current intensity from the probes 
Eigen::VectorXi get_current_intensity( const Eigen::MatrixXi& X_probes, int Ncolors, int* Nprobes, int Nribs ){
    Eigen::VectorXi current_intensity;

    current_intensity.setZero(Ncolors);

    int mm = 0;
    for (int i = 0; i < Nribs; i ++){
        mm = 0;
        for (int j = 0; j < Ncolors; j++){
            for (int k = 0; k < Nprobes[j]; k++){
                current_intensity(j) += X_probes(i, mm);
                mm++; 
            }
        }
    }
    return current_intensity;
}

void generic_ssa_cpp(int* result, int* intensity, int* states, int* Stoich_states, int* Stoich_lattice, double* forward_rates, double* parameters, int* xi_lattice, int* xi_state, double* time_vector, 
                    double tf, int seed, int Nt, int n_rxns, int n_total_rxns, int n_states, int length, int max_particles, int used_frames, int* probe_location_matrix,   int Ncolors){


    int* Nprobes = new int[Ncolors];
    // reconstruct the probe information
    int k = 0;
    int probe_val = 0;
    int n_total_probes = 0;
    for (int i = 0; i < Ncolors; i++){
        probe_val = 0;
        for (int j = 0; j < length; j++){

            if (probe_location_matrix[k] == 1){
                probe_val++;
                n_total_probes++;
            }
            k++;
        }
        Nprobes[i] = probe_val;
    }

    int* probe_locations =  new int[n_total_probes];
    k = 0;
    int m = 0;
    for (int i = 0; i < Ncolors; i++){
        probe_val = 0;
        for (int  j = 0; j < length; j++){

            if (probe_location_matrix[k] == 1){
                probe_locations[m] = j;  
                m++;
            }
            k++;
        }

    }


	std::mt19937_64 rng; // initalize twister RNG
    rng.seed(seed);
	std::uniform_real_distribution<double> unif(0, 1); // define the uniform dist

    double t = 0; // time 
    t=time_vector[0]; // set the time to the first time point
    //tf = time_vector[Nt-1];  // get the final time point from Nt

    int reaction_ind = 0; // reaction index
    int time_ind = 0; // time index for recording
    int N_particles = 0; // number of ribosomes in the system at a time
    int N_particles_old = 0;
    int R = 10;  // ribosomal exclusion (its minus one so 10 == 9)
    int rib_stepped_too = 0;
    double probe_binding_rate = 1;

    //X_lattice - binary vector of 1 or 0 for ribosome occupying this location for all viable locations in the model 
    MatrixXi X_lat(1,length);	// make initial X_lattice vector
    X_lat.setZero();
    
    // initialize the X_lattice to xi_lattice
	for(int i=0; i < length; i++)	{
		X_lat(0,i) = xi_lattice[i];
	}
    N_particles = X_lat.sum(); // get the number of particles starting off
	//Eigen::MatrixXd Xdouble = X.cast<double>(); // cast to double for dot product 

    // X compressed - vector of 1 to max ribosomes recording the indexes of where x_lattice == 1
    Eigen::MatrixXi X_comp(1,max_particles);
    X_comp.setZero();

    if (N_particles > 0){
    compress_X_locations(X_comp, X_lat, length, N_particles); // get X_compressed
    }

    int frame_Ls = 0;
    if (used_frames == 1){
        frame_Ls = length;
    }
    if (used_frames == 2){
        frame_Ls = (length+1)/2;
    }
    if (used_frames == 3){
        frame_Ls = (length+2)/3 ;
    }

    Eigen::MatrixXi X_probes(max_particles, n_total_probes);
    X_probes.setZero();


    Eigen::MatrixXi X_spatial(1,max_particles);
    X_spatial.setZero();
    if (N_particles > 0){
        update_spatial_X(X_comp, X_spatial, N_particles, frame_Ls);  // update spatial X vector 
    }
    

    MatrixXi X_states(1,n_states);	// make initial x state vector
    X_states.setZero();

	for(int i=0; i < n_states; i++)	{
		X_states(0,i) = xi_state[i];
	}

	// convert things to eigen format
	
	double a0, r1, r2;
	
	Eigen::Map<Eigen::MatrixXi> X_array(result,Nt,max_particles);
    Eigen::Map<Eigen::MatrixXi> X_intensity(intensity,Nt,Ncolors);
    Eigen::Map<Eigen::MatrixXi> states_record(states,Nt,n_states);


    Eigen::Map<Eigen::MatrixXi> S_state( Stoich_states, n_states, n_rxns ); // have to read in row wise this is  dumb
    Eigen::Map<Eigen::MatrixXi> S_lat( Stoich_lattice, length, n_total_rxns - n_rxns ); // have to read in row wise this is  dumb


    Eigen::MatrixXi St_state(n_rxns,n_states);
    St_state << S_state.transpose();

    Eigen::MatrixXi St_lat(n_total_rxns - n_rxns, length) ;
    St_lat << S_lat.transpose();

    MatrixXd forward_rate_matrix(1,length);	// make forward_rates
    forward_rate_matrix.setZero(1,length);

	for(int i=0; i < length; i++)	{
		forward_rate_matrix(0,i) = forward_rates[i];
	}


	while (t < tf){
		
		// record reactions:	
		
		// generate new propensities:

        Eigen::VectorXd wn(1,n_total_rxns + N_particles);
   
        wn = propensity_function(X_comp,X_states, X_spatial, X_lat, parameters, forward_rates, t, n_total_rxns, N_particles, R, length, probe_binding_rate );
		a0 = wn.sum();

        // Generate some random numbers.
		
        r1 =  unif(rng);
        r2 =  unif(rng);

		while((r1==0)){
			r1 =  unif(rng);			
		}
		
        // Update the time vector
        t -= log(r1)/a0;

        reaction_ind = 1;
        while (wn.head(reaction_ind).sum() < r2*a0)
        {	
            reaction_ind +=1;
        }

        while( (time_ind<Nt) && (t>time_vector[time_ind])) {

            X_array.row(time_ind) << X_comp.row(0);
            X_intensity.row(time_ind) << get_current_intensity( X_probes, Ncolors, Nprobes, N_particles  ).transpose();
            states_record.row(time_ind) << X_states.row(0);

            time_ind+=1;
			
         }	

        if (reaction_ind-1 >= n_total_rxns){  // ribosome took a step
            X_comp(0,reaction_ind-n_total_rxns-1 ) +=1;
            rib_stepped_too = X_comp(0,reaction_ind-n_total_rxns-1 ) -1;
            update_spatial_X(X_comp, X_spatial, N_particles, frame_Ls);  // update spatial X vector 
            decompress_X_locations(X_comp, X_lat, length, N_particles  );  // update X_lat
            // update the probe_matrix if a ribosome stepped
            update_probe_matrix_step(X_probes, X_comp, probe_locations, probe_binding_rate, n_total_probes , rib_stepped_too, reaction_ind-n_total_rxns-1, unif(rng));

            
        }
        else{
            if (reaction_ind-1 >= n_rxns){  // not a state changing reaction
                N_particles_old = N_particles;
                X_lat.topLeftCorner(1,length) = X_lat.topLeftCorner(1,length) + St_lat.row(reaction_ind-1 - n_rxns); // x_lat_stoich
                N_particles = X_lat.sum(); // update number of particles

                // update the probe vector before we update xcomp, we need to know which ribosome is leaving (if ones leaving)
                if (N_particles_old > N_particles){
                    int leaving_reaction = reaction_ind-1 - n_rxns;
                    int index_leaving = 0;
                    int row_to_delete = 0;
                    for (m = 0; m < length; m++){
                        if ( St_lat( leaving_reaction, m) == -1){
                            index_leaving = m;
                        }
                    }
                    for (m = 0; m < N_particles_old; m++){
                        if (index_leaving == X_comp(0, m) - 1){
                            row_to_delete = m;
                        }
                    }

                    update_probe_matrix_leave(X_probes, row_to_delete, N_particles, n_total_probes);
                }
                compress_X_locations(X_comp, X_lat, length, max_particles  );  // update X_compressed
                update_spatial_X(X_comp, X_spatial, N_particles, frame_Ls);  // update spatial X vector 



            }
            else{   // state changed
                X_states.topLeftCorner(1,n_states) = X_states.topLeftCorner(1,n_states) + St_state.row(reaction_ind-1); // update the state vector accordingly
            }

        }

	}


}
    
