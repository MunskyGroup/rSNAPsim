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

void translationSSA_lowmem(double* kelong, double* t_array, int Nt, double kbind, double kcompl, int* SSA_intensity, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double* SSA_ribtimes, int* watched_ribs, int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp, int* x0, int r_footprint, int* SSA_probe, int Ncolor, int* flags, double kon, double koff, double* k_probe,int* probe_loc,int* n_probes, int rib_max)
{
    // Declare the variables
	
	/*
	Inputs explaination: 
	
	kelong - a double 1xGeneLength vector of each codon at postition rate
	t_array - a time array of when to record the SSA_result
	Nt - the size of the t_array int
	kbind - ki the rate of initation
	kcompl - ktermination the rate of unbinding / release
	SSA_result - the result of the ssa, preallocated from python (N time points x 200 ribosomes x N trajectories)
	
	N - number trajectories
	FRAP - true false to do frap at inhibit time
	
	Inhibitor - True false to inhibit k_in at inhibit time
	inhibit_time - time of inhibition
	
	
	seed - RNG seed passed from python to ensure it doesnt use the same seeds for groups of trajectories - preallocated by python random vector
	rib_times - an array to be filled with the times of a watched_ribs to watch how long they take to go through the simulation
	watched_ribs - the number of ribosomes to watch, this is to get around a dynamic array; basically the code will just watch "X" amount of ribosomes and record their data
	        instead of trying to do a dynamic array
	ribtimesize - the size of the ribosome watching array
	
	so to do the frap the simulation just returns the section that underwent frap, then the python wrapper will do the logic to bleach that section, as it gets complex 
	with ribosomes that were partially bleached
	
	fNt - the number of frap times, basically we need to record the section of the simulation that undergoes frap to manually "bleach it later" as the logic is complex
	frap_result - the area of the ssa result that undergoes frap 
	
	cNt - the size of the collision recording array 
	col_result - records collisions
	
	colNp - number of collisions to record as well - vectorized
	col_x - posistions of collisions
	col_t - time of said collision
	

	
	
	
	
	*/
	
    int R = r_footprint; // ribosome exclusion.
    int N_rib = rib_max; // maximum number of ribosomes. 

	
	std::mt19937_64 rng; 
	//INCLUDE <#chrono> to seed from computer clock
    // initialize the random number generator with time-dependent seed
    //uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    //std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(seed);
	std::uniform_real_distribution<double> unif(0, 1);
	
	
	
	
	MatrixXi leaky_probe_matrix(Ncolor*N_rib,N);	// setup stuff for flags may not be used
	leaky_probe_matrix.setZero(Ncolor*N_rib,N);

	
	MatrixXi col(1,N_rib);
	col.setZero(1,N_rib);
	
	MatrixXi T(1,N_rib); //ribosome travel time array
	T.setZero(1,N_rib);
	
	//VectorXd T_array(200);
	int t_counter = 0;
	int col_counter = 0;
	
	// Create an eigen matrix that stores the results. 
	Eigen::Map<Eigen::VectorXd> T_array(SSA_ribtimes,ribtimesize);  // this map function will fill the python allocated array
	Eigen::Map<Eigen::VectorXi> col_array(col_result,cNt);
	Eigen::Map<Eigen::VectorXd> colarrayt(col_t,colNp);
	Eigen::Map<Eigen::VectorXi> colarrayx(col_x,colNp);
	Eigen::Map<Eigen::VectorXi> n_ribs(watched_ribs,1);
	int tsize = T_array.size();
	int col_size = colarrayt.size();		

    int it = 0;
	int number_ribs = 0;
	int fit = 0;
	
	int bursting = flags[0];
	int leaky = flags[1];
	int stats = flags[2];
	
    // int N = 10; // gene length
    //bool Inhibitor = 0;
    //bool FRAP = 0;
    double Inhibit_condition=1;
    bool inhibit_off=1;
    //double inhibit_time = 0;
	
	bool inhibitor_pres=0;
	bool FRAP_pres =0;
	
	if (Inhibitor==1){
		inhibitor_pres =1;
	}
	if (FRAP==1){
		FRAP_pres =1;
	}	
    int NR = 0;
	int old_NR = 0;
//	N = 10;
    double a0, r1, r2, t, tf;
    t=t_array[0];
    tf = t_array[Nt-1];
    int ind = 0;
	
	int burst = 1;
	if (bursting){
		
		r1 = unif(rng);
		int burst = 0; // is the system bursting on or off_type with proprotion to kon / koff
		if (r1 < kon / (kon+koff)){
			int burst = 1;
		}
	}
	// predefine stuff that may remain blank
	int total_probes = 0;
	for (int i = 0; i < Ncolor; i++){
		total_probes += n_probes[i];
	}
	

	
    MatrixXi X(1,N_rib);	
    X.setZero(1,N_rib);

	for(int i=0; i < N_rib; i++)	{
		X(0,i) = x0[i];

	}
	
	
	Eigen::Map<Eigen::VectorXi> X_array(SSA_intensity,Nt,Ncolor);
	Eigen::Map<Eigen::MatrixXi> frap_array(frap_result,fNt,N_rib);
	
	
	
	MatrixXi probe(Ncolor,N);
	probe.setZero(Ncolor,N);
	for (int i=0; i < Ncolor; i++){
		for (int j=0; j < N; j++){
			probe(i,j) = SSA_probe[i*N + j ]; 
			//std::cout << i << "," << j << " " << i*N + j << " " << SSA_probe[i*N + j ]  << std::endl;
		}
	}
	



    while( t < tf)
    {

        // Determine inhibitor stuff
        if (inhibitor_pres) {
            if (t>=inhibit_time){
                inhibit_off = 0;
                Inhibit_condition = 0;
            } else {
                inhibit_off = 1; 
                Inhibit_condition=1;
            }}
        else { 
            Inhibit_condition=1;
        }


		
        // Update the number of ribosomes, always making sure to have
        // a zero on the right side of X (space for a new ribosome) 
		
        int NR = 0;
        while (X(0,NR)>0){
            NR+=1;
        }
		
		if (stats){
			if (NR > old_NR){

				T(0,NR-1) = t;
				number_ribs +=1;
				
			}
		}
		
		
		if (leaky){
			if (NR> old_NR){
				for (int k=0; k < total_probes*2; k+=2){
					
					r1 = unif(rng);
					if ( r1 > 1-k_probe[probe_loc[k]]  ){
						leaky_probe_matrix(Ncolor*(NR-1) + probe_loc[k] ,probe_loc[k+1]) = 1;
					}
				}
			}
			//std::cout << "adding probes: " << leaky_probe_matrix.block(0,0,NR+3,25)  << std::endl;
		}


        old_NR = 0;
        while (X(0,old_NR)>0){
            old_NR+=1;
			
        }
			
		
				
        //std::cout << "X: " << X << std::endl;
        MatrixXi Sn(NR+1,NR+1);
        Sn.setIdentity(NR+1,NR+1);
        //std::cout << "Stoichiometry Matrix: \n" << Sn << std::endl;
        Eigen::VectorXd wn;
		if (bursting){
			wn.setZero(NR+2);
		}
		else{
			wn.setZero(NR+1);
		}
		// for loop instead of "where"       
        // also, account for the exclusion here.
        for(int i=0; i <= NR; i++)
        {
            if( X(0,i) > 0 ){
                wn(i) = kelong[X(0,i)-1];
            }
            if( i>0){
                if( X(0,i-1)-X(0,i)<=R ){
					
                    wn(i) = 0;
                }
            }
        }
		
		if (bursting){
			if (burst == 1){	
				wn(NR+1) = koff;	
			}
			else{
				wn(NR+1) = kon;
			}
		}

        // If a nascent protein reaches full length
        // add in the completion rate
//       std::cout  << X << std::endl;
        if( X(0,0) == N ){
            // shift everyone if completing
            for(int i=0; i <= NR-1; i++)
            {
                Sn(0,i) = X(0,i+1) - X(0,i); 
            }
            Sn(0,NR-1) = -X(0,NR-1);
            wn(0) = kcompl;
//            std::cout << "Updated propensity: \n" << wn << std::endl;
//            std::cout << "Updated stoichiometry: \n " << Sn << std::endl;
        }

        // include the initiation condition
        if( (NR==0) || ( X(0,NR-1)>R) ){
            wn(NR) = kbind*Inhibit_condition*burst;
        }
        // Update the propensity
        a0 = wn.sum();

        // Generate some random numbers.
		
        r1 =  unif(rng);
        r2 =  unif(rng);
		
		// if rand() gets a 0 resample, since ln(0) = -inf 
		if((r1==0)){
			r1 = unif(rng);			
		}
		
		
        t -= log(r1)/a0; // time of next reaction

		

		
        while( (it<=Nt-1) && (t>t_array[it])) {

			for(int j = 0; j < Ncolor; j++){

				int intensity = 0;
				for(int i=0; i < NR; i++){

					if (leaky){
						intensity += leaky_probe_matrix.block( Ncolor*i  +j,0,1,X(0,i)-1).sum();
						//std::cout << "intensity: " << intensity << std::endl;
						//std::cout << "current matrix: " << leaky_probe_matrix.block(0,0,NR+3,25)  << std::endl;
					}
					else{
						//std::cout << "intensity: " << intensity << std::endl;
						intensity += probe(j,X(0,i)-1);
					}

				}

				X_array(it,j) = intensity;
				
				
				if (FRAP_pres){
					
					if ( (t>= inhibit_time) && (t< inhibit_time+20)) {
						frap_array.row(fit) = X.row(0);
						fit +=1;
					}
				}			
			}
            it+=1;
			//std::cout << it << std::endl;
			
         }
		 

		
        // update the state
        ind = 1;
        while (wn.head(ind).sum() < r2*a0)
        {	
            ind +=1;
        }
		
		if (bursting){
			if (ind == NR+2){
				burst+=1; // flip the bursting state
				burst = burst%2;
			}
			else{
				X.topLeftCorner(1,NR+1) = X.topLeftCorner(1,NR+1) + Sn.row(ind-1);
			}
		}
		else{
			X.topLeftCorner(1,NR+1) = X.topLeftCorner(1,NR+1) + Sn.row(ind-1);
		}
		
		if (leaky){
			if ( (Sn.row(ind-1).sum() < 0) && (ind != NR+2)){
				//std::cout << "removing row: " << leaky_probe_matrix.block(0,0,NR+3,25)  << std::endl;
				leaky_probe_matrix.block(0,0,NR*Ncolor,N) = leaky_probe_matrix.block(Ncolor,0,Ncolor*NR+Ncolor,N);
				//std::cout << "removing row: " << leaky_probe_matrix.block(0,0,NR+3,25)  << std::endl;
			}
		}
		
		if (stats){
			if (Sn.row(ind-1).sum() < 0){

				if (t_counter < tsize){
					T_array(t_counter) = t - T(0,0);
							
					T.block(0,0,1,NR) = T.block(0,1,1,NR);
					T(0,NR) = 0;	
					
					col_array(t_counter) = col(0,0);	
					col.block(0,0,1,NR) = col.block(0,1,1,NR);	
					col(0,NR)= 0;
					t_counter +=1;		
					
				}
			}
			else {
				
				
				if (X(0,ind-2) == X(0,ind-1) + R){
					col(0,ind-1) +=1;
					if (col_counter < col_size){
						colarrayt(col_counter) = t;
						colarrayx(col_counter) = X(0,ind-1);
						col_counter+=1;
					}
					
				}
				
			}
		}
		//std::cout << "oneiter" << std::endl;

    }
	if (stats){
		n_ribs(0) = number_ribs;
	}
	
        
}
 
/*
int main()
{
    //double k[2][5] = {{1, 2, 3, 4, 5},{6,7,8,9,10}};
    //double k[1][10] = {{3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5}};
    double k[10] = {3.1,3.2,3.3,3.4,3.5,3.1,3.2,3.3,3.4,3.5};
    double t_array[3] = {0, 50,100};
    double kbind = 9.0; 
    double kcompl = 15.0; 
	double inhibit_time = 10.0;
	int FRAP = 0;
	int Inhibitor = 1;
    int Nt;
    int result[10];
	double ribtimes[50];
	int seed;
    //seed = srand (time(NULL));
    Nt = (sizeof(t_array)/sizeof(*t_array));
    for(int kk=0; kk<3; kk++){
        translationSSA(k, t_array, 3, kbind, kcompl, result, 10, FRAP, Inhibitor, inhibit_time,15,ribtimes);
    }
}   

*/
