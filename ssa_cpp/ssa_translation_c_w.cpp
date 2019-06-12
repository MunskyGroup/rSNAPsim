#include <iostream>
#include <ctime>
#include <math.h>

#include <eigen/Eigen/Dense>
#include <cstdlib>
#include <random>
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

void translationSSA(double* kelong, double* t_array, int Nt, double kbind, double kcompl, int* SSA_result, int N, int FRAP, int Inhibitor, double inhibit_time, int seed, double* SSA_ribtimes, int* nribs, int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result, double* col_t, int* col_x, int colNp )
{
    // Declare the variables
    int R = 9; // ribosome exclusion.
    int N_rib = 200; // maximum number of ribosomes. 
	srand(seed);
    int it = 0;
	int number_ribs = 0;
	int fit = 0;
	std::cout << fNt << std::endl;
	
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
    //srand (1537478928);  
	//std::cout << time(NULL) << std::endl;
    // print a test random number
    // Define the state vector
    MatrixXi X(1,N_rib);
    X.setZero(1,N_rib);
	
	MatrixXi col(1,N_rib);
	col.setZero(1,N_rib);
	
    MatrixXi T(1,N_rib); //ribosome travel time array
    T.setZero(1,N_rib);
	
	//VectorXd T_array(200);
	int t_counter = 0;
	int col_counter = 0;
	
    // Create an eigen matrix that stores the results. 
    Eigen::Map<Eigen::MatrixXi> X_array(SSA_result,Nt,N_rib);
	
	Eigen::Map<Eigen::MatrixXi> frap_array(frap_result,fNt,N_rib);
	
	
	
	
	Eigen::Map<Eigen::VectorXd> T_array(SSA_ribtimes,ribtimesize);
	Eigen::Map<Eigen::VectorXi> col_array(col_result,cNt);
	
	
	Eigen::Map<Eigen::VectorXd> colarrayt(col_t,colNp);
	Eigen::Map<Eigen::VectorXi> colarrayx(col_x,colNp);
	
	Eigen::Map<Eigen::VectorXi> n_ribs(nribs,1);

	int tsize = T_array.size();
	int col_size = colarrayt.size();

    while( t < tf)
    {
        //std::cout << "-------t=" << t << "-------" << std::endl;
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
		

		
		
		if (NR > old_NR){

			T(0,NR-1) = t;
			number_ribs +=1;
			
			
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
        wn.setZero(NR+1);
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
            wn(NR) = kbind*Inhibit_condition;
        }
        // Update the propensity
        a0 = wn.sum();

        // Generate some random numbers.
		
        r1 =  ((double) rand()/ (RAND_MAX));
        r2 =  ((double) rand() / (RAND_MAX));
		
		
		// if rand() gets a 0 resample, since ln(0) = -inf 
		if((r1==0)){
			//std::cout << r1 << " " << r1a << " " << t <<  std::endl;
			
			r1 =  ((double) rand() / (RAND_MAX));			
		}
		
		

        // Update the time vector
//        std::cout << "TIME " << t << std::endl;
        t -= log(r1)/a0;
		//std::cout << t << std::endl;
/* 		if((t > inhibit_time) && (t < inhibit_time+30)){
			std::cout << "TIME " << t << std::endl;
			std::cout << wn << std::endl;
			std::cout << X << std::endl;
		}
		 */
		
        //std::cout << "TIME " << t << std::endl;
  //      std::cout << "-------------" << std::endl;
        // Store the X vector
        //while( ((it<=Nt-1) || (t>t_array[it])) ){
        while( (it<=Nt-1) && (t>t_array[it])) {
			//std::cout << it << std::endl;
            X_array.row(it) = X.row(0);
			
			if (FRAP_pres){
				
				if ( (t>= inhibit_time) && (t< inhibit_time+20)) {
					frap_array.row(fit) = X.row(0);
					fit +=1;
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
		
        //std::cout << "Stoichiometry of reaction: " << Sn.row(ind-1) << std::endl;
        
		
		
        X.topLeftCorner(1,NR+1) = X.topLeftCorner(1,NR+1) + Sn.row(ind-1);
		//std::cout << "rxn "<< Sn.row(ind-1).sum() << std::endl;
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
		//std::cout << "oneiter" << std::endl;

    }
	n_ribs(0) = number_ribs;
	
        
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
