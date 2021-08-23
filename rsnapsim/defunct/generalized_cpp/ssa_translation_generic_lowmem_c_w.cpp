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

//void translationSSA_general(double* kelong, double* t_array, double kbind, double kcompl, int* SSA_result, int N, int* stop_sites, int FRAP, int Inhibitor, double inhibit_time, int seed, double* SSA_ribtimes, int* nribs, int ribtimesize, int fNt, int* frap_result, int cNt, int* col_result,int ASL, double ki_alt)

void translationSSA_generic(double* kelong, double* t_array, int* SSA_result, int N, int Nt, double* inhibitors, int seed, int fNt, int* frap_result, double* k_add, int n_enters,int n_pauses,int n_stops, int n_jumps,int* SSA_probe, int Ncolor, int Nlocs, int watched_ribs )
{
    // Declare the variables
	
	//std::cout << "init"  << "-------" << std::endl;
	
    int R = 9; // ribosome exclusion.
    int N_rib = 200; // maximum number of ribosomes. 
	srand(seed);
	
    int it = 0;
	//int number_ribs = 0;
	int fit = 0;
	
	int nrxns = 0;
	//N is gene length
	int additional_rxns_cnt = 0;
	int binds = 0;
	int add_index = 0;
	int add_index_2 = 0;
	bool extra_reactions = false;
	int removed_rib = 0;
	//int Nt = sizeof(t_array)/sizeof(*t_array); // get number of time points
	//std::cout << Nlocs << std::endl;
	MatrixXi probe(Ncolor,Nlocs);
	probe.setZero(Ncolor,Nlocs);
	//std::cout << "seting up probe"  << Nlocs << std::endl;
	//std::cout << "seting up probe"  << Ncolor << std::endl;
	
	for (int i=0; i < Ncolor; i++){
		for (int j=0; j < Nlocs; j++){
			probe(i,j) = SSA_probe[i*Nlocs + j ]; 
			//std::cout << i << "," << j << " " << i*Nlocs + j << " " << SSA_probe[i*Nlocs + j ]  << std::endl;
		}
	}
	//std::cout << "finished probe"  << Ncolor << std::endl;
	int nprobes = 0;
	for (int i=0; i < Ncolor; i++){
		nprobes += probe.row(i).sum();
	}
	//std::cout << "probe " << probe << "-------" << std::endl;
	MatrixXi probe_locs(2,nprobes);
	probe_locs.setZero(2,nprobes);
	//std::cout << "probe " << probe_locs << "-------" << std::endl;
	
	int nn = 0;
	for (int i=0; i < Ncolor; i++){
		for (int j=0; j < Nlocs; j++){
			if (probe(i,j) > 0){
				probe_locs(1,nn) = j; 
				probe_locs(0,nn) = i; 
				//std::cout << "probe " << probe_locs << "-------" << std::endl;
				nn+=1;
			}
			//std::cout << i << "," << j << " " << i*Nlocs + j << " " << SSA_probe[i*Nlocs + j ]  << std::endl;
		}
	}
	
	MatrixXi intensity_mat(Ncolor,watched_ribs);
	intensity_mat.setZero(Ncolor,watched_ribs);
	
	
	//std::cout << "probe " << probe_locs << "-------" << std::endl;
	////std::cout << "im " << intensity_mat << "-------" << std::endl;
	//std::cout << "n_pauses" << n_pauses << "-------" << std::endl;
	//std::cout << "n_stops" << n_stops << "-------" << std::endl;
	//std::cout << "n_enters" << n_enters << "-------" << std::endl;
	//std::cout << "n_enters" << n_jumps << "-------" << std::endl;
	
	Eigen::MatrixXd jumps(n_jumps,3);  // set up matrices of these reactions 
	Eigen::MatrixXd pauses(n_pauses,2);
	Eigen::MatrixXd stops(n_stops,2);
	Eigen::MatrixXd enters(n_enters,2);
	
	
	jumps.setZero();
	pauses.setZero();
	stops.setZero();
	enters.setZero();

/* 	Eigen::VectorXd k_enters;
	Eigen::VectorXd k_pauses;
	Eigen::VectorXd k_stops;
	Eigen::VectorXd k_jumps;
 */
	
	Eigen::VectorXd wn_bind;
	Eigen::VectorXd wn_pauses;
	Eigen::VectorXd wn_jumps;
	Eigen::VectorXd wn_stops;
	
	
	Eigen::MatrixXi jump_events;
	Eigen::MatrixXi stop_events;
	Eigen::MatrixXi pause_events;
	Eigen::MatrixXi bind_events;
	
	
	int n_rxns = n_jumps+n_enters+n_pauses+n_stops;
	
	Eigen::VectorXd important_locs;
	important_locs.setZero(n_jumps+n_pauses+n_stops);
	
	// first map all the passed new reactions to matrices for easier handling
	
	int k = 0;
	int j = 0;
	int m = 0;
	
	

	for(int i=0; i < n_enters; i++){ 
	
		enters(j,0) = k_add[m];
		m+= 1;
		enters(j,1) = k_add[m];
		m+= 1;
		j += 1;
		
	}
	

	
	
	j = 0;
	for(int i=0; i < n_pauses; i++){
		
		pauses(j,0) = k_add[m];
		m+= 1;
		pauses(j,1) = k_add[m];	
		m+= 1;
		important_locs(k) = pauses(j,0);	
		j += 1;		
		k += 1;
	}	

	j = 0;
	for(int i=0; i < n_stops; i++){
				
		stops(j,0) = k_add[m];
		m+= 1;
		stops(j,1) = k_add[m];
		m+= 1;
		important_locs(k) = stops(j,0);
		j += 1;		
		k += 1;
	}	

	
	j = 0;
	for(int i=0; i < n_jumps; i++){
		
		jumps(j,0) = k_add[m];
		m+= 1;
		jumps(j,1) = k_add[m];
		m+= 1;
		jumps(j,2) = k_add[m];
		m+= 1;
		important_locs(k) = jumps(j,0);

		j += 1;		
		k += 1;
		
	}

  	//std::cout << "important locs  " << important_locs << std::endl;
	//std::cout << "enters  " << enters << std::endl;
	//std::cout << "pauses  " << pauses << std::endl;
	//std::cout << "stops  " << stops  << std::endl;
	//std::cout << "jumps  " << jumps  << std::endl; 
 
	
    // int N = 10; // gene length
    //bool Inhibitor = 0;
    //bool FRAP = 0;
    double Inhibit_condition=1;
    bool inhibit_off=1;
	//double inhibit_time
    double inhibit_time = 0;
	
	bool inhibitor_pres=false;
	bool FRAP_pres =false;
	
	// inhibitors = [ harringtonine-bool, FRAP-bool, time]
	//std::cout << "inhibs  " << inhibitors[0]  << std::endl;
	if (inhibitors[0]==1){
		inhibitor_pres =true;
	}
	if (inhibitors[1]==1){
		FRAP_pres =true;
	}	
	
	if (inhibitor_pres || FRAP_pres){
		inhibit_time = inhibitors[2];
	}		
	
	
    int NR = 0;
	int old_NR = 0;
//	N = 10;
    double a0, r1, r2, t, tf;
    t=t_array[0];
    tf = t_array[Nt-1];
    int ind = 0;
	//int transferX;
    //srand (1537478928);  
	//std::cout << time(NULL) << std::endl;
    // print a test random number
    // Define the state vector
    MatrixXi X(1,N_rib);
    X.setZero(1,N_rib);
	
/* 	MatrixXi col(1,N_rib);
	col.setZero(1,N_rib);
	
    MatrixXi T(1,N_rib); //ribosome travel time array
    T.setZero(1,N_rib);
	 */
	//VectorXd T_array(200);
	int t_counter = 0;
	
    // Create an eigen matrix that stores the results. 
    Eigen::Map<Eigen::MatrixXi> X_array(SSA_result,Nt,Ncolor);
	
	Eigen::Map<Eigen::MatrixXi> frap_array(frap_result,fNt,N_rib);
	
	
	int asl_free = 0;
	
/* 	Eigen::Map<Eigen::VectorXd> T_array(SSA_ribtimes,ribtimesize);
	Eigen::Map<Eigen::VectorXi> col_array(col_result,cNt);
	Eigen::Map<Eigen::VectorXi> n_ribs(nribs,1); */
	
	MatrixXi spatial_X(1,N_rib);

	//int tsize = T_array.size();
	
	//std::cout << Nt << std::endl;
	//std::cout << t_array << std::endl;
	//std::cout << sizeof(t_array) << std::endl;
    while( t < tf)
    {   //std::cout << "new iteration" << "-------" << std::endl;
        //std::cout << "-------t=" << t << "-------" << std::endl;
		//std::cout << "-------tf=" << tf << "-------" << std::endl;
		//std::cout << "x=" << X.topLeftCorner(1,5) << "-------" << std::endl;
		
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
		
		
/* 		if (NR > old_NR){

			T(0,NR-1) = t;
			number_ribs +=1;
			
			
		} */


        old_NR = 0;
        while (X(0,old_NR)>0){
            old_NR+=1;
			
        }
			
		nrxns = 0;
		binds = 0;
		
		spatial_X.setZero();
		for (int i = 0; i <= NR; i++){
			spatial_X(0,i) = (X(0,i) + (X(0,i)>N) + (X(0,i) > (2*N-1)) )%N; 
			
		}
		//std::cout << "spatial_X " << spatial_X << "-------" << std::endl;
		//spatial_X.topLeftCorner(1,NR+1) == (X.array() - (X.array()>N)-(X.array() > (2*N-1)) )%N; // convert to spatial locs
		
		// parse the binds first because they change the stoich column count
		additional_rxns_cnt = 0;
		
		
		

		wn_bind.resize(0);
		wn_jumps.resize(0);
		wn_stops.resize(0);
		wn_pauses.resize(0);
		
        
		jump_events.resize(0,0);
		stop_events.resize(0,0);
		pause_events.resize(0,0);
		bind_events.resize(0,0);
	
		for (int i = 0; i <n_enters; i++){
			
			//bool loc_free = spatial_X.head(NR).array()- enters(i,0)
			int j = enters(i,0);
			int enter_loc = (j + (j>N) + (j > (2*N-1)) )%N;
			//std::cout  << enter_loc << "----" <<std::endl;
			//std::cout  << j << "----"<<std::endl;
			
			bool loc_not_free = (((spatial_X.block(0,0,1,NR).array())-enter_loc).abs() < R).any();
			
			//std::cout  << " check "<< ((spatial_X.block(0,0,1,NR).array())-enter_loc).abs() << "----" <<std::endl;
			
			
			//std::cout << (((spatial_X.block(0,0,1,NR).array())-enter_loc).abs() < R) << std::endl;
			//std::cout  << loc_not_free << "----"<<std::endl;
			
			//std::cout << "spatial x " << ((spatial_X.block(0,0,1,NR).array())-j).abs() << "----"<<std::endl;
			//std::cout <<"loc free flag " << loc_not_free << "----"<<std::endl;
			

			// is this location free?
			if (loc_not_free == false){
				//std::cout  << binds << "----"<<std::endl;
				
				if (binds == 0){
					//Eigen::MatrixXd pause_events(1,NR+1);
					bind_events.resize(1,NR+1);
					bind_events.setZero();
					bind_events(0,NR) = enters(i,0);
					
					additional_rxns_cnt+=1;
					wn_bind.conservativeResize(1);
					wn_bind(0) = enters(i,1);
					binds +=1;
					
				}
				else{
					bind_events.conservativeResize(bind_events.rows()+1,bind_events.cols());
					bind_events.row(binds).setZero();
					
					bind_events(binds,NR) = enters(i,0);
					additional_rxns_cnt+=1;
					wn_bind.conservativeResize(wn_bind.size()+1);
					wn_bind(wn_bind.size()-1) = enters(i,1);
					binds +=1;
					//std::cout  << bind_events << "----"<<std::endl;
					
				}
			}
		}
		
		//std::cout  << bind_events << "----"<<std::endl;
		
		extra_reactions = false;
		for (int i = 0; i <= n_jumps+n_pauses+n_stops; i++){
			if ((X.block(0,0,1,NR).array()== important_locs(i)).any()){
				extra_reactions = true;
				//std::cout << "adding reaction "  << extra_reactions << "----"<<std::endl;
			}			
		}
		
		
		bool rib_on_pause = false;
		if (extra_reactions == true){ 
			//handle the pauses now
			
			
			for (int i = 0; i <n_pauses; i++){
				bool pause_check = (X.block(0,0,1,NR).array()== pauses(i,0)).any(); 
				if (pause_check){
					rib_on_pause = true;
				}
				//std::cout  <<"pause check "<< pause_check << "----"<<std::endl;	
				//std::cout  <<"on pause "<< rib_on_pause << "----"<<std::endl;	
			}
			
			
			
/* 			add_index = 0;
			add_index_2 = 0;
			
			
		
			
			for (int i = 0; i <n_pauses; i++){
				bool rib_on_pause = (X.block(0,0,1,NR).array()== pauses(i,0)).any();   //if there is a ribosome in the pause position
				//std::cout  <<"on pause "<< rib_on_pause << "----"<<std::endl;
				//std::cout  <<"pause loc "<< pauses(i,0) << "----"<<std::endl;
				
				if (rib_on_pause== true){
					
					
					//for (int j = 0; j <=n_pauses; j++){
						//if ( abs(X(j)-pauses(i,0)) > R){
						//	add_index = j
						//}
					//}
					
					if (add_index == 0){
						//Eigen::MatrixXd pause_events(1,NR+binds);
						pause_events.resize(1,NR+(binds>0));
						pause_events.setZero();
						add_index += 1;
						additional_rxns_cnt+=1;
						wn_pauses.conservativeResize(1);
						wn_pauses(0) = pauses(i,1);
						//std::cout  <<"wn p "<< wn_pauses << "----"<<std::endl;
						
						//std::cout  <<"S p "<< pause_events << "----"<<std::endl;

					}
					else{
						pause_events.conservativeResize(pause_events.rows()+1,pause_events.cols());
						pause_events.row(add_index).setZero();
						add_index += 1;
						additional_rxns_cnt+=1;
						wn_pauses.conservativeResize(wn_pauses.size()+1);
						wn_pauses(wn_pauses.size()-1) = pauses(i,1);						
					}
				}
				
				
			}
 */

			// handle jumps
			add_index = 0;
			add_index_2 = 0;
			
			
			
			for (int i = 0; i <n_jumps; i++){
				bool on_jump = (X.block(0,0,1,NR).array()== jumps(i,0)).any(); 
				
				add_index = 0;
				
				
				if (on_jump == true){
					
					
					for (int j = 0; j <=NR; j++){
						if ( X(0,j) == jumps(i,0) ){
							add_index_2 = j;
						}
					}
					
					
					int jump_loc = (X(0,add_index_2) + (X(0,add_index_2)>N) + (X(0,add_index_2) > (2*N-1)) )%N;
					
					bool jump_ex_1 = ((spatial_X.block(0,0,1,add_index_2).array() - jump_loc).abs() < R).any();
					bool jump_ex_2 = false;
					
					if (add_index_2 < NR){
						jump_ex_2 = ((spatial_X.block(0,add_index_2+1,1,NR).array() - jump_loc).abs() < R ).any();
					}
					
					
					if (jump_ex_1==false && jump_ex_1==false){ //not excluded, this thing can jump there
					
						if (add_index == 0){
							//Eigen::MatrixXd jump_events(1,NR+binds);
							jump_events.resize(1,NR+(binds>0));
							jump_events.setZero();
							jump_events(add_index,add_index_2) = jumps(i,1)-jumps(i,0);
							add_index +=1;
							additional_rxns_cnt+=1;
							wn_jumps.resize(1);
							wn_jumps(0) = jumps(i,2);			
							
							//std::cout  <<"wn_jumps "<< wn_jumps << "----"<<std::endl;
							//std::cout  <<"add_index_2 "<< add_index_2 << "----"<<std::endl;
							//std::cout  <<"jump_events "<< jump_events << "----"<<std::endl;	
							
						}
						else{
							jump_events.conservativeResize(jump_events.rows()+1,jump_events.cols());
							jump_events.row(add_index).setZero();
							jump_events(add_index,add_index_2) = jumps(i,1)-jumps(i,0);
							add_index +=1;
							additional_rxns_cnt+=1;
							wn_jumps.conservativeResize(wn_jumps.size()+1);
							wn_jumps(wn_jumps.size()-1) = jumps(i,2);							
						}
					}
				}
				
				
			}

			//handle stops
			add_index = 0;
			add_index_2 = 0;
			
			
			for (int i = 0; i < n_stops; i++){
				bool on_stop = (X.block(0,0,1,NR).array()== stops(i,0)).any(); 
				//std::cout << "wn stops i " << i << "----"<<std::endl;
				add_index = 0;
				if (on_stop == true){
					
					
					for (int j = 0; j < NR; j++){
						if ( X(0,j) == stops(i,0)){
							add_index_2 = j;
						}
					}
					
					if (add_index == 0){
						//Eigen::MatrixXd stop_events(1,NR+binds);
						stop_events.resize(1,NR+ (binds>0)  );
						stop_events.setZero();
						stop_events(add_index,add_index_2) = -stops(i,0);
						add_index +=1;
						additional_rxns_cnt+=1;
						wn_stops.resize(1);
						
						wn_stops(0) = stops(i,1);
						//std::cout << "add_index " << add_index << "----"<<std::endl;
						//std::cout << "wn stops " << wn_stops << "----"<<std::endl;
						//std::cout << "stops  " << stop_events << "----"<<std::endl;
						//std::cout << "stop prop " << stops(i,1) << "----"<<std::endl;
						//std::cout << "add_index_2  " << add_index_2 << "----"<<std::endl;

						
					}
					else{
						stop_events.conservativeResize(stop_events.rows()+1,stop_events.cols());
						stop_events.row(add_index).setZero();
						stop_events(add_index,add_index_2) = -stops(i,0);
						add_index +=1;
						additional_rxns_cnt+=1;
						wn_stops.conservativeResize(wn_stops.size()+1);
						wn_stops(wn_stops.size()-1) = stops(i,1);	
						//std::cout << "add_index " << add_index << "----"<<std::endl;
						//std::cout << "wn stops " << wn_stops << "----"<<std::endl;
						//std::cout << "stops  " << stop_events << "----"<<std::endl;		
						//std::cout << "add_index " << add_index << "----"<<std::endl;
						//std::cout << "wn stops " << wn_stops << "----"<<std::endl;
						//std::cout << "stops  " << stop_events << "----"<<std::endl;
						//std::cout << "stop prop " << stops(i,1) << "----"<<std::endl;
						//std::cout << "add_index_2  " << add_index_2 << "----"<<std::endl;

												
						
					}
				}
				
				
			}
		}
		//std::cout << "wn stops " << wn_stops << "----"<<std::endl;
				
        //std::cout << "X: " << X << std::endl;
		
		//std::cout << "making sn_elong " << "----"<<std::endl;
        MatrixXi Sn_elong(NR,NR+(binds>0));
		Sn_elong.setZero();
        Sn_elong.topLeftCorner(NR,NR).setIdentity();
        
        Eigen::VectorXd wn;
        wn.setZero(NR);
		
		//std::cout << "sn_elong" << Sn_elong << "----"<<std::endl;
		
		//std::cout << "making sn_elong done " << "----"<<std::endl;
		
        // for loop instead of "where"       
        // also, account for the exclusion here.
		
		//std::cout << "NR " << NR << "----"<<std::endl;

		
        for(int i=0; i < NR; i++)
        {
			
			
            if( X(0,i) > 0 ){
                wn(i) = kelong[X(0,i)-1];
            }
            if( i>0){
                if( X(0,i-1)%N-X(0,i)%N<R ){
                    wn(i) = 0;
                }
                
            }
			
			if (rib_on_pause){  // if previously detected that a ribosome is on a pause, check location, if its a pause multiply by rate
				
				//std::cout << "rib_on_pause " << rib_on_pause << "----"<<std::endl;
				for (int j = 0; j <n_pauses; j++){
					if (X(0,i) == pauses(j,0)){
						//std::cout << " old wn(i) " << wn(i) << "----"<<std::endl;
						wn(i) = kelong[X(0,i)-1]*pauses(j,1);
						//std::cout << "wn(i) " << wn(i) << "----"<<std::endl;
					}
				}
			}
			
			//std::cout << "wn(i) " << wn(i) << "----"<<std::endl;
        }
		
		Eigen::VectorXd wn_all_rxns(wn.size() + wn_bind.size() + wn_pauses.size()+ wn_jumps.size()+ wn_stops.size());
		wn_all_rxns << wn, wn_bind, wn_pauses,wn_jumps,wn_stops;
		
		
		Eigen::MatrixXi Sn(Sn_elong.rows()+additional_rxns_cnt,Sn_elong.cols());
		Sn.setZero();
		//std::cout  <<"all b4 "<< Sn << "----"<<std::endl;
		Sn << Sn_elong,bind_events,pause_events,jump_events,stop_events;
		
		//std::cout  <<"X "<< X << "----"<<std::endl;
		//std::cout  <<"t "<< t << "----"<<std::endl;
		//std::cout  <<"spatial X "<< spatial_X << "----"<<std::endl;
		
		//std::cout  <<"all "<< Sn << "----"<<std::endl;
		//std::cout << "seed" << seed << "-------" << std::endl;
		//std::cout << "Elong S" << Sn_elong << "----"<<std::endl;
		
		
		//std::cout << "bindsS " << bind_events << "----"<<std::endl;
		
		//std::cout << "stopsS " << stop_events << "----"<<std::endl;
		
		//std::cout << "jumpS  " << jump_events << "----"<<std::endl;
		//std::cout << "PauseS  " << pause_events << "----"<<std::endl;
		
		//std::cout << "wn " << wn_all_rxns << "----"<<std::endl;
		
		//std::cout << "wn stops " << wn_stops << "----"<<std::endl;
		

		
		
        // If a nascent protein reaches full length
        // add in the completion rate
//       std::cout  << X << std::endl;
/*         if( ((X(0,0)%N) == 0) && (X(0,0) !=0) ){
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
        if( (NR==0) || ( X(0,NR-1) > R) ){
            wn(NR) = kbind*Inhibit_condition;
        }
		 */
		

        // Update the propensity
        a0 = wn_all_rxns.sum();

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
			//std::cout << t << std::endl;
			//std::cout << "starting intensity" << std::endl;
           
			
			
			for(int j=0; j < Ncolor; j++){
/*     			int intensity = 0;
    			for(int i=0; i < NR; i++){
    
    				//std::cout << " i, j=" << i << "," <<j << std::endl;

    				intensity += probe(j,X(0,i)-1);
                    //std::cout << "-------sum=" << leaky_probe_matrix.block( Ncolor*i  +j,0,1,X(0,i)-1).sum()<< "-------" << std::endl;
    				
    			}
    			 */
				int intensity = 0;
				//std::cout << "-------sum=" << intensity_mat.topLeftCorner(3,5)<< "-------" << std::endl;
				//std::cout << "----- Ncolor" << Ncolor << std::endl;
				//std::cout << "-------sum=" << intensity_mat.row(j).sum()<< "-------" << std::endl;

				intensity = intensity_mat.row(j).sum();
                X_array(it,j) = intensity;
                //std::cout << " int = " << t << " " <<intensity << std::endl;
            }			

			
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
        while (wn_all_rxns.head(ind).sum() < r2*a0)
        {	
            ind +=1;
        }
		
        //std::cout << "Stoichiometry of reaction: " << Sn.row(ind-1) << std::endl;
        //std::cout << "Stoichiometry: " << Sn << std::endl;
		//std::cout << "X: " << X.topLeftCorner(1,5) << std::endl;


		//Sn.row(ind-1).size();
		
        X.topLeftCorner(1,Sn.row(ind-1).size()) = X.topLeftCorner(1,Sn.row(ind-1).size()) + Sn.row(ind-1);
		
		//std::cout << "rxn S added "<< Sn.row(ind-1) << std::endl;


		
		//std::cout << "rxn "<< Sn.row(ind-1).sum() << std::endl;
		//std::cout  <<"ind "<< ind << "----"<<std::endl;
		if (ind >= NR+1){  // alternative reaction
			//std::cout  <<"alt reaction selected "<< ind << "----"<<std::endl;
			
			if (Sn.row(ind-1).sum() < 0){ //something left the system
			
				
			
				//std::cout  <<"negative reaction " << "----"<<std::endl;
/* 				if (t_counter < tsize){
					T_array(t_counter) = t - T(0,0);
							
					T.block(0,0,1,NR) = T.block(0,1,1,NR);
					T(0,NR) = 0;	
					
					col_array(t_counter) = col(0,0);	
					col.block(0,0,1,NR) = col.block(0,1,1,NR);	
					col(0,NR)= 0;
					t_counter +=1;		
					
				} */
				int removed_rib = 0;
				//std::cout  <<"old X " << X<< "----"<<std::endl;
				for(int i =0; i<=NR;i++){
					if (X(0,i) == 0){
						removed_rib = i;
						break;
						
					}
				}
				
				//std::cout  <<"removed_rib  " << removed_rib << "----"<<std::endl;
				//std::cout  <<"num_ribs  " << NR << "----"<<std::endl;
				
				intensity_mat.block(0,removed_rib,Ncolor,NR) = intensity_mat.block(0,removed_rib+1,Ncolor,NR);
				for(int i=removed_rib; i <= NR; i++){// shift everyone if completing
					
					X(0,i) = X(0,i+1); 
					//std::cout  <<"removed_im before  " << intensity_mat << "----"<<std::endl;
					//intensity_mat.block(0,i,Ncolor,NR) = intensity_mat.block(0,i+1,Ncolor,NR);
					//std::cout  <<"removed_im  " << intensity_mat << "----"<<std::endl;
				}
				//std::cout  <<"new X " << X<< "----"<<std::endl;
			}
		}
		else{           // something moved forward, add probes
			//std::cout << "checking probes "  << "-------" << std::endl;
			int x_changed = ind-1;
/* 			for (int i = 0; i <= Sn.row(ind-1).size(); i++){
				if (Sn.row(ind-1)(i) == 1){
					x_changed = i;
					
				} 
				//std::cout << "NR "  << NR << "-------" << std::endl;
			} */
			for (int j = 0; j < nprobes; j++){ 
				//std::cout << "j "  << j << "-------" << std::endl;
				//std::cout << "probespot "  << X(0, x_changed) << " " << probe_locs(1,j) << "-------" << std::endl;
				
				if (X(0, x_changed) == probe_locs(1,j)){
					
					//std::cout << "xchange " << x_changed << "-------" << std::endl;
					//std::cout << "ind " << ind << "-------" << std::endl;
					
					//std::cout << "stoich " << Sn.row(ind-1) << "-------" << std::endl;
					//std::cout << "pl " << probe_locs(1,j) << "-------" << std::endl;
					//std::cout << "im before " << X(0,x_changed) << "-------" << std::endl;
					//std::cout << "im before " << intensity_mat << "-------" << std::endl;
					intensity_mat(probe_locs(0,j),x_changed) += 1;
					//std::cout << "im after  " << intensity_mat << "-------" << std::endl;
				}
					
			
			}
			//std::cout << "done with probe "  << "-------" << std::endl;
			
		}
		
		
/* 		if((ind <= NR+1)){  // if something elongated and hit another ribosome, record it.
			
			if (X(0,ind-2) == X(0,ind-1) + R){
				col(0,ind-1) +=1;
			}
			
		} */
/* 			if (Sn.row(ind-1).sum() < 0){  // something left the system with this reaction

				if (t_counter < tsize){
					T_array(t_counter) = t - T(0,0);
							
					T.block(0,0,1,NR) = T.block(0,1,1,NR);
					T(0,NR) = 0;	
					
					col_array(t_counter) = col(0,0);	
					col.block(0,0,1,NR) = col.block(0,1,1,NR);	
					col(0,NR)= 0;
					t_counter +=1;		
					
				}
				
				
				for(int i =0; i<=NR;i++){
					if (X(0,i) == 0){
						int removed_rib = i;
						
					}
				for(int i=removed_rib; i <= NR; i++){// shift everyone if completing
					
					X(0,i) = X(0,i+1); 
						
					
					}
				
                if( ((X(0,0)%N) == 0) && (X(0,0) !=0) ){
            // shift everyone if completing
				for(int i=0; i <= NR-1; i++)
				{
					Sn(0,i) = X(0,i+1) - X(0,i); 
				}				
					 */
					
/* 				}
				else {
					if (X(0,ind-2) == X(0,ind-1) + R){
						col(0,ind-1) +=1;
					}


			}
			//std::cout << "oneiter" << std::endl;
		} */
/* 		transferX = -1;
		if((ind == NR+2)){
			//std::cout << "altstart" << std::endl;
			//std::cout << X << std::endl;
			if((NR != 0)){
				transferX = -1;
				for(int i=0; i <= NR; i++){
					
					
					
					
					if (transferX == -1){
						//std::cout << "transferx" << std::endl;
						if((X(0,i)) < ASL && (X(0,i))!=0 ){
							transferX = 1;
							
							
							for(int j=NR+1; j >= i; j--){
								
								//std::cout << X << std::endl;
								X(0,j+1) = X(0,j);
								
							}
							
							//std::cout << X << std::endl;
							X(0,i) = ASL+N+1;
							X(0,NR+1) = 0;
							//std::cout << "adding rib" << std::endl;
							//std::cout << X << std::endl;
							
							
							
						}
						
					}
				}
			}
			else{
				X(0,0) = ASL+N+1;
				
			} */

					
			
				
					
			//std::cout << X << std::endl;
		
	
		
		
    
	
		

		} 
	//n_ribs(0) = number_ribs;		
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
