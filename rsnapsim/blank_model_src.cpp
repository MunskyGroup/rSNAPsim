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
using Eigen::VectorXi;
using namespace std;


int min_int(int a, int b)
{
    if (a == b) {
        return a;
    }
    if (a < b) {
        return a;
    } 
    if (a > b){
        return b;
    }
    return a;
}


double min_double(double a, double b)
{
    if (a == b) {
        return a;
    }
    if (a < b) {
        return a;
    } 
    if (a > b){
        return b;
    }
    return a;
}



/*
propensity here
*/
void BLANKPROP( Eigen::VectorXd& wn, double* parameters, double t, const Eigen::MatrixXi& rib_arr, 
                                     const Eigen::MatrixXd& elong_mat, const Eigen::VectorXi& occupied,
                                     const Eigen::VectorXi& lattice_arr, const Eigen::MatrixXi& probe_mat,
                                     const Eigen::VectorXi& state_arr, const Eigen::VectorXi& resource_arr, int NR,
                                    int n_constant_rxns, int n_ribosome_rxns, int L, int footprint, int max_rib,
                                    int rib_arr_col_size,
                                    int n_states,
                                    int n_resources,
                                    int n_parameters,
                                    double tc){

    wn.setZero();
    int k;
     
    //INSERT_GENERATED_CONSTANT_PROPENSITY_HERE

    k = n_constant_rxns;
    for (int i = 0; i < NR; i++){
        //for (int j = 0; j < n_ribosome_rxns; j++){

            //INSERT_GENERATED_RIBOSOME_PROPENSITY_HERE


        //}
        k += n_ribosome_rxns;
    }

}


/*
REACTION MATRIX STRUCTURE:
Nrxns rows by Ncolors + Nstates + Nresources + 7

rxn row: [TYPE, EXCLUDE?, FRAME, LOC, Dexist, Dframe, Dloc, Dprobe1.... DprobeN, Dstate1... DstateN, Dresource1.... DresourceN]


RIBOSOME MATRIX STRUCTURE (agent array of all particles on the tasep)

Max_ribosomes by Nrxns + Ncolors + 4

ribosome row: [ID, Exists 0 or 1, frame, location, Ncolor1 ... Ncolor2, N times rxn1 .... N times rxn N]

*/



void generic_ssa_cpp(int* particle_array, int* state_array, int* resource_array, // arrays to fill
                    int* particle_x0, int* state_x0,  int* resource_x0, int* rxns, int* rxnids, double* elongs, int* probes,
                    double* parameters, int npars,       // initial simulation state
                    int max_rib, int n_constant_reactions, int n_ribosome_reactions,  // constants of the simulation
                    int n_colors, int n_states, int n_resources, int Nt, double* time_vector, double tf, double burnin, int seed, int L, int error_code){


 //int* Stoich_states, int* Stoich_lattice, double* forward_rates,
 // double* parameters, int* xi_lattice, int* xi_state, double* time_vector, 
 //                   double tf, int seed, int Nt, int n_rxns, int n_total_rxns,
  //                   int n_states, int length, int max_particles,
   //                   int used_frames, int* probe_location_matrix,
    //                  int Ncolors, double burnin){

    /*
    Initialize constants
    footprint = size of a particle
    max_rib = maximum number of particles
    n_colors = number of colors
    n_states = number of states
    n_resources = number of resources
    n_constant_reactions = number of constant reactions
    n_ribosome_reactions = number of per ribosome reactions
    L = length of the mRNA
    */

   int footprint = 9;           // PARTICLE SIZE
   int rib_ind = 0;             // ribosome that reactions are happening to

   //int reaction_taken = 0; //
   //int dexist = 0;
   int ribosome_moved = 0;      // boolean did a ribosome move this simulation step
   int ribosome_left = 0;       // boolean did a ribosome leave this simulation step
   int tindex = 0;              // index for storing simulation data back to python arrays

   double tc = 0 - burnin;      // current time
   int NR = 0;                  // number of ribosomes
   int k = 0;                   //
   int event = 0;               // which reaction is happening
   int n_rxns = n_constant_reactions + n_ribosome_reactions; // total number of reactions
   int rib_arr_size = n_colors+4+n_rxns;                     // Size of the particle/agent array
   int rxn_arr_size = n_colors+7+n_states+n_resources;       // size of the reaction matrix (agent instruction array)

   int pf = 1; // probe function, set to 1 for now, this can be a function late

    /*
    Initialize all vectors / matrices:
    rib_arr = particle array of ribosomes moving (agents)
    lattice_arr = occupation lattice (1 x L)
    occupied = list of ribosome locations (1 x max_ribosomes)
    state_arr = state vector
    resource_arr = resource vector

    */


// make initial ribosome_array / dependent arrays: lattice and occupied
    MatrixXi rib_arr(max_rib, 4 + n_colors + n_rxns);	
    rib_arr.setZero();

    VectorXi lattice_arr(L);	
    lattice_arr.setZero();

    VectorXi occupied(max_rib);	
    occupied.setZero();

    // read in the particle initial state
    k = 0;
    for (int i = 0; i < max_rib; i++){
        for (int j = 0; j < (4 + n_colors + n_rxns); j++){
            rib_arr(i, j) = particle_x0[k];
            k++;
        }
    }    
    NR = rib_arr.col(1).sum(); // get the number of particles starting off

    occupied << rib_arr.col(3); // occupation vector is the third column of the agent array
    for (int i = 0; i < max_rib; i++){ // for each non zero occupation value, set the lattice arr location to 1
        if (occupied[i] > 0){
            lattice_arr[occupied[i]] = 1;
        }
    }


// make initial state_array
    VectorXi state_arr(n_states);	
    state_arr.setZero();
    // read in the particle initial state
    k = 0;
    for (int i = 0; i < n_states; i++){
        state_arr[i] = state_x0[k];
        k++;
    }    

// make initial resource_array
    VectorXi resource_arr(n_resources);	
    resource_arr.setZero();
    // read in the particle initial state
    k = 0;
    for (int i = 0; i < n_resources; i++){
        resource_arr[i] = resource_x0[k];
        k++;
    }    


    /*
    initialize Reaction matrix, elongation matrix, probe_matrix
    */

    MatrixXi rxn_mat(n_rxns, 7 + n_colors + n_states + n_resources);	
    rxn_mat.setZero();
    k = 0;
    for (int i = 0; i < n_rxns; i++){
        for (int j = 0; j < (7 + n_colors + n_states + n_resources); j++){
            rxn_mat(i, j) = rxns[k] ;
            k++;
        }
    }    

// the list of reaction ids, maps rows in the reaction matrix to propensites
    VectorXi rxn_ids(n_rxns);
    k = 0;
    for (int i = 0; i < n_rxns; i++){
        rxn_ids[i] = rxnids[k];
        k++;
    }    


// default stepping rate matrix for all 3 frames of the mRNA
    MatrixXd elong_mat(3, L);	
    elong_mat.setZero();
    k = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < L; j++){
            elong_mat(i, j) = elongs[k] ;
            k++;
        }
    }    


// Probe location matrix, 0 if no probe on this codon, 1-N for different colors
    MatrixXi probe_mat(3, L);	
    probe_mat.setZero();
    k = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < L; j++){
            probe_mat(i, j) = probes[k] ;
            k++;
        }
    }    

    // SETUP RNG
	std::mt19937_64 rng; // initalize twister RNG
    rng.seed(seed);
	std::uniform_real_distribution<double> unif(0, 1); // define the uniform dist

    // doubles for calculating next reaction and time
	double a0, r1, r2;

	// map the vectors to fill and return to python
	Eigen::Map<Eigen::MatrixXi> Particle_array(particle_array, Nt, (4 + n_colors + n_rxns)*max_rib);
    Eigen::Map<Eigen::MatrixXi> State_array(state_array,n_states,Nt);
    Eigen::Map<Eigen::MatrixXi> Resource_array(resource_array,n_resources,Nt);

    Eigen::VectorXd wn(n_constant_reactions + n_ribosome_reactions*max_rib); // initalize the largest possible propensity
    wn.setZero();

	while (tc < tf){
		
        // CHECK IF ANYTHING WENT WRONG IN THE SIMULATION AND ERROR OUT WITH A CODE!
        if ((rib_arr.array() < 0).any() || (state_arr.array() < 0).any() 	|| (resource_arr.array() < 0).any() ){
            std::cout << "ERROR NEGATIVE PARTICLES, STATES, OR RESOURCES!!! " << std::endl;
            error_code = -1;
            return;
        }

        // GENERATE NEW PROPENSITIES
        BLANKPROP(wn, parameters,
                                 tc, 
                                 rib_arr,
                                 elong_mat,
                                 occupied,
                                 lattice_arr,
                                 probe_mat,
                                 state_arr,
                                 resource_arr,
                                 NR,
                                 n_constant_reactions, n_ribosome_reactions, L, footprint, 
                                 max_rib,
                                 rib_arr_size,
                                 n_states,
                                 n_resources,
                                 npars,
                                 tc);
                               

        // CHECK IF ANYTHING WENT WRONG IN THE SIMULATION AND ERROR OUT WITH A CODE!
        if ((wn.array() < 0).any() ){
            std::cout << "ERROR NEGATIVE PROPENSITY RATES!!! Returning..." << std::endl;
            std::cout << seed << std::endl;
            std::cout << "rib_arr" << std::endl;
            std::cout << "______________" << std::endl;
            std::cout << rib_arr << std::endl;
            std::cout << "running propensity..." << std::endl;
            std::cout << "wn" << std::endl;
            std::cout << "______________" << std::endl;
            std::cout << wn << std::endl;
            error_code = -2;
            return;
        }

		//check particle locations and error out 

        // Generate 2 random numbers.
        r1 =  unif(rng);

        // MAKE SURE r1 IS NOT ZERO OR THIS WILL CRASH
		while((r1==0)){
			r1 =  unif(rng);			
		}

        // sum of the propensities
		a0 = wn.sum();
		if (a0 == 0){ // special case if all reactions are zero, end the simulation, fast_rxn 
                	 // should be added as an option in the future!! TODO
    		tc = time_vector[Nt-1] + 1; // increment the time to be past the recording
    		// this will fill the remainder of the arrays with the previous state
    		}
		
		else{
            // Update the time vector, what time did the next reaction happen?
            tc -= log(r1)/a0;
        }

        r2 =  unif(rng);

        // fill up recording matrixes if time passed current time index
        while( (tindex < Nt) && (tc > time_vector[tindex])) {
            Eigen::Map<VectorXi> v(rib_arr.data(),rib_arr.size());
            Particle_array.row(tindex) = v;
            State_array.col(tindex) << state_arr;
            Resource_array.col(tindex) << resource_arr;
            tindex +=1;
            if (tindex == Nt){  // manually end the while loop if tc > time_vector[-1]
                error_code = 0; // 0 means ran successfully
                return;
            }
         }	

        // figure out which event happened based on unif r2
        event = 0;
        while (wn.head(event).sum() < r2*a0)
        {	
            event +=1;
        }
        event -=1;

    
        /*
        if (event >= n_constant_reactions){ //if its a ribosome reaction, on which ribosome did it occur
            rib_ind = (event-n_constant_reactions)%NR ; // WHICH RIBOSOME IS THIS REACTION HAPPENING TOO
            event = ((event-n_constant_reactions)/NR) + n_constant_reactions; // WHICH REACTION IS HAPPENING
        }
        */
        if (event >= n_constant_reactions){ //if its a ribosome reaction, on which ribosome did it occur
            rib_ind = ((event-n_constant_reactions)/n_ribosome_reactions); // WHICH REACTION IS HAPPENING
            event = (event-n_constant_reactions)%n_ribosome_reactions + n_constant_reactions ; // WHICH RIBOSOME IS THIS REACTION HAPPENING TOO

        }        


        event = rxn_ids[event]; // edit event to match rxn matrix since it can be in any order
        


        

        
        ribosome_moved = 0; // reset moved/left booleans
        ribosome_left = 0;

        if (rxn_mat(event,0) == 2){ // lattice reaction

                //fr = rxn_mat[event][2]
                //loc = rxn_mat[event][3]
                //rib_ind = np.where(rib_arr[:,3] == loc)[0][0]

            // find which ribosome the reaction is happening to based on location
            rib_ind = 0;
            for (int i = 0; NR < 0; i++){
                if ( (rib_arr(i, 2) == rxn_mat(event, 2)) && (rib_arr(i, 3) == rxn_mat(event, 3)) ){
                    rib_ind = i;
                }
            }

            if (rxn_mat(event, 4) == 1){ // RIBOSOME ARRIVING
                // set the ribosome id to 1+ current ribosome count
                rib_ind = NR;
                rib_ind += 1; // fill in the particle array
                rib_arr(NR, 0) = rib_ind;
                rib_arr(NR, 1) = 1; // exists? flag
                rib_arr(NR, 2) += rxn_mat(event, 2);
                rib_arr(NR, 3) += rxn_mat(event, 3);
                rib_arr(NR, 4 + n_colors + event) += 1;
                NR += 1;
                ribosome_moved = 1; // set boolean to 1, ribosomes changed in someway
            } 
            else if (rxn_mat(event, 4) == -1){ // RIBOSOME LEAVING

                rib_arr.block(rib_ind, 0, NR - rib_ind, rib_arr_size) << rib_arr.block(rib_ind+1, 0, NR - rib_ind, rib_arr_size); // shift all rows below this up
                rib_arr.row(NR).setZero();
                NR -=1;
                ribosome_left = 1;
            }
            else{ // RIBOSOME MOVING
                rib_arr.row(rib_ind).block(0,1,1,n_colors+4) << rib_arr.row(rib_ind).block(0,1,1,n_colors+4) + rxn_mat.row(event).block(0, 4, 1, n_colors+3);
                rib_arr(rib_ind,4+n_colors+event) += 1;

            }
            // if there are states, update states based on reaction taken
            if (rxn_mat(event,1) == 1){
                if (n_states > 0){
                    state_arr = state_arr +  rxn_mat.row(event).block(0, n_colors+7, 1, n_states).transpose();
                }
                if (n_resources > 0){// if there are resources, update states based on reaction taken
                    resource_arr = resource_arr +  rxn_mat.row(event).block(0, n_colors+7+n_states, 1, n_resources).transpose();
                }
            }
            if (rxn_mat(event,6) !=0){ // MOVING ALONG LATTICE
                ribosome_moved = 1;
            } 
            if (rxn_mat(event,5) !=0){ // MOVING FRAME
                ribosome_moved = 1;
            }

        }

        if (rxn_mat(event,0) == 0){ // ribosome reaction

            // place the dexist, dframe, dloc, dcolors
            rib_arr.row(rib_ind).block(0,1,1,n_colors+4) << rib_arr.row(rib_ind).block(0,1,1,n_colors+3) + rxn_mat.row(event).block(0, 4, 1, n_colors+3); // place the dexist, dframe, dloc, dcolors
            rib_arr(rib_ind,4+n_colors+event) += 1; // update the reaction that happened counter
            
            if (rxn_mat(event,1) == 1){
                // update states if needed
                if (n_states > 0){
                    state_arr = state_arr +  rxn_mat.row(event).block(0, n_colors+7, 1, n_states).transpose();
                }
    
                // update resources if needed
                if (n_resources > 0){
                    resource_arr = resource_arr +  rxn_mat.row(event).block(0, n_colors+7+n_states, 1, n_resources).transpose();
                    }
            }
            
            
            ribosome_moved = 0;
            if (rxn_mat(event,6) != 0){ // ribosome moved in lattice (dlocation !=0 )
                ribosome_moved = 1;
            }
            if (rxn_mat(event,5) != 0){ // ribosome moved in frames (dframe !=0 )
                ribosome_moved = 1;
            }
            ribosome_left = 0;
            if (rxn_mat(event,4) == -1){ // Ribosome left (dexist = -1)
                rib_arr.block(rib_ind, 0, NR - rib_ind, rib_arr_size) << rib_arr.block(rib_ind+1, 0, NR - rib_ind, rib_arr_size); // shift all rows below this up
                rib_arr.row(NR).setZero();
                NR -=1;
            }
            
        
        }        

        if (rxn_mat(event,0) == 1){ // state reaction
            state_arr = state_arr +  rxn_mat.row(event).block(0, n_colors+7, 1, n_states).transpose();
        }

        
        if (rxn_mat(event,0) == 3){ // resource reaction
            resource_arr = resource_arr +  rxn_mat.row(event).block(0, n_colors+7+n_states, 1, n_resources).transpose();
        }

        if (ribosome_left == 1){
            // a ribosome left so we have to recalculate the occupation and lattice arrays
            occupied << rib_arr.col(3); // occupation vector is the third column of the agent array
            lattice_arr.setZero();
            for (int i = 0; i < NR; i++){ // for each non zero occupation value, set the lattice arr location to 1
                if (occupied[i] > 0){
                    lattice_arr[occupied[i]] = 1;
                }
            }
        }

        if (ribosome_moved== 1){ // if the ribosome moved update the occupation vector
            // TODO: THIS CAN BE SPED UP AND MADE MORE EFFICIENT IN THE FUTURE**
            occupied << rib_arr.col(3); // occupation vector is the third column of the agent array
            lattice_arr.setZero();
            for (int i = 0; i < NR; i++){ // for each non zero occupation value, set the lattice arr location to 1
                if (occupied[i] > 0){
                    lattice_arr[occupied[i]] = 1;
                }
            }
            if (probe_mat(rib_arr(rib_ind,2), rib_arr(rib_ind,3)) != 0){
                if (pf){
                    rib_arr(rib_ind, 3 +  probe_mat(rib_arr(rib_ind,2), rib_arr(rib_ind,3) )) +=1;
                }
            }
        }
    }

}
    