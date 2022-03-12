


void writing(double t, std::string outTimeStamp){


    std::vector<std::vector<double>> u(Nx+1, 
        std::vector<double> (N_HP + Ncool + 1));
    
    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < N_HP + 1; j++){
            u[i][j] = u_HP[i][j];
        }
        for (int j = 0; j < Ncool + 1; j++){
            u[i][j + N_HP] = u_cool[i][j];
        }
    }


	std::ifstream root_dir("app_root_dir", std::ifstream::in);
	std::string app_root_dir; 
	root_dir >> app_root_dir;
    std::string n = std::to_string(t*dt);
    std::string folder = "results\\results_";
    std::string folderName = app_root_dir + folder + outTimeStamp;
    _mkdir( folderName.c_str() );


    /* ====================== _T2Dt.csv ====================== */
	std::ofstream myfile; 
    

    myfile.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_T2Dt.csv");
    
    if(!myfile.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_T2Dt.csv";
    }

    for (int i = 0; i < Nx + 1; i++){
        for (int j = 0; j < N_HP + Ncool + 1; j++){
            myfile << u[i][j] << "; ";
        }
    myfile << "\n";
    }

    myfile.close();

    /* ====================== _u_TEG_N.csv ====================== */
	std::ofstream myfile_1; 


    myfile_1.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_u_TEG_N.csv");
    
    if(!myfile_1.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_u_TEG_N.csv";
    }
    for (int i = 0; i < Nx + 1; i++){
        for (int j = 0; j < Nx_TEG + 1; j++){
            myfile_1 << u_TEG_N_main[i][j] << "; ";
        }
    myfile_1 << "\n";
    }

    myfile_1.close();

    /* ====================== _u_TEG_P.csv ====================== */
	std::ofstream myfile_2; 

    myfile_2.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_u_TEG_P.csv");
    
    if(!myfile_2.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_u_TEG_P.csv";
    }

    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < Nx_TEG + 1; j++){
            myfile_2 << u_TEG_P_main[i][j] << "; ";
        }
    myfile_2 << "\n";
    }

    myfile_2.close();

    /* ====================== _uLoo.csv ====================== */
	std::ofstream myfile_3; 


    myfile_3.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_uLoo.csv");
    
    if(!myfile_3.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_uLoo.csv";
    }

    for (int i = 0; i < Nx+1; i++){
        myfile_3 << uoo[i] << "; ";
    }

    myfile_3.close();

    /* ====================== _uLcc.csv ====================== */
	std::ofstream myfile_4; 


    myfile_4.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_uLcc.csv");
    
    if(!myfile_4.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_uLcc.csv";
    }

    for (int i = 0; i < Nx+1; i++){
        myfile_4 << ucc[i] << "; ";
    }

    myfile_4.close();

    /* ====================== _Q_exc.csv ====================== */
	std::ofstream myfile_5; 


    myfile_5.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Q_exc.csv");
    
    if(!myfile_5.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Q_exc.csv";
    }
    
    for (int i = 0; i < Nx+1; i++){
        myfile_5 << ucc[i] << "; ";
    }

    myfile_5.close();

    /* ====================== _Qoo.csv ====================== */
	std::ofstream myfile_6; 


    myfile_6.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Qoo.csv");
    
    if(!myfile_6.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Qoo.csv";
    }
    
    for (int i = 0; i < Nx+1; i++){
        myfile_6 << Qoo[i] << "; ";
    }

    myfile_6.close();
            
    /* ====================== _Qcc.csv ====================== */
	std::ofstream myfile_7; 


    myfile_7.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Qcc.csv");
    
    if(!myfile_7.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_Qcc.csv";
    }
    
    for (int i = 0; i < Nx+1; i++){
        myfile_7 << Qcc[i] << "; ";
    }

    myfile_7.close();

    /* ====================== _I.csv ====================== */
	std::ofstream myfile_8; 


    myfile_8.open (app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_I.csv");
    
    if(!myfile_8.is_open()){ 
        std::cout << std::endl << "Error writing: " + app_root_dir + "results\\results_" + outTimeStamp + "\\" + n.substr(0, n.size()-5) + "_I.csv";
    }
    
    myfile_8 << I_current;

    myfile_8.close();

}


