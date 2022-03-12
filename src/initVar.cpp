

double getVar(std::string file, std::string idVar){


	/*if you input two strings it gives the the 
	variable value insie the file*/

	std::string id;
	double var;

	std::ifstream root_dir("app_root_dir", std::ifstream::in);
	std::string app_root_dir; 
	root_dir >> app_root_dir;

	std::ifstream theFile(app_root_dir + "input_data\\controlVars\\" + file);

    if(!theFile.is_open()){ 
        std::cout << "ERROR Opening: input_data\\controlVars\\" << file;
    }

    while(theFile >> id >> var){

    	if (id == idVar){
    		theFile.close();
    		return var;
    	}
    }

}


void init_n_print(){

	std::cout << std::endl << 
	"=================================================="
	"==============================" << std::endl;
	std::cout << "Rui Vieira, 2019" << std::endl;
	std::cout << "HP_HT2Dt_TE1Dt" << std::endl;
	std::cout << "version 10" << std::endl;
	std::cout << "2D transient simulation of heat transfer" 
	" with thermoelectric effects" << std::endl;

    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "Started computation on " << std::ctime(&end_time);

	std::cout <<
	"=================================================="
	"==============================" << std::endl << std::endl;
}



void ETA(auto start_time, double time_par, double final_time_par, 
	double error_par, int it, double write_increment){

    std::chrono::duration<double> elapsed_time = 
        std::chrono::system_clock::now()-timeElapsed1_old;

	std::cout << std::endl << std::endl << 
	"=================================================="
	"==============================" << std::endl << std::endl;
	std::cout << "Writing step: " << time_par << "/" 
	<< final_time_par << " [s]" << std::endl;
	std::cout << "Error= " << error_par;

	if (write_increment == time_par+1){


		std::cout << " in " << it << " iterations; " <<
		-((elapsed_time.count())*( 1/write_increment)*(time_par-final_time_par))/(60*60)
		<< "hrs remaining (estimated)" << std::endl;

	}
	else{

		std::cout << " in " << it << " iterations for " <<
		elapsed_time.count() << "s; ETA " <<
		-( ( elapsed_time.count() )*( 1/write_increment)*
			(time_par-final_time_par))/(60*60) << " hrs" << std::endl;

	}

    timeElapsed1_old = std::chrono::system_clock::now();
}


