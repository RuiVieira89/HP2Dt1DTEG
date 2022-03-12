
/* ======= SOLVER ======= */

void HP_2Dtsolve(){ 

    std::cout << "Delta t = " << dt << " [s]" << "\n";
    std::cout << "\n" << "initializing... " << std::endl << std::endl;

    /*======= Results folder timestamp =======*/

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H-%M-%S",timeinfo);
    std::string outTimeStamp(buffer);

    /*======= Input cycle.csv initial conditions =======*/

	std::ifstream root_dir("app_root_dir", std::ifstream::in);
	std::string app_root_dir; 
	root_dir >> app_root_dir;

	std::ifstream theFile(app_root_dir + "input_data\\cycle\\cycle.csv");

    if(!theFile.is_open()){ 
        std::cout << std::endl << "Error opening: input_data\\cycle\\cycle.csv";
    }

	double a;
	double b;
    double c;
    char d;


    int size=0;

    while(theFile >> a >> d >> b >> d >> c){  
        
        size += 1;
    }

    theFile.clear();
    theFile.seekg(0, theFile.beg);


	double pExhaustPower;
	double pExhaustMass;
    double pExhaustTemperature;
    char comma;

	std::vector<double> ExhaustPower(size);
	std::vector<double> ExhaustMass(size);
    std::vector<double> ExhaustTemperature(size);

    int i = 0;

    while(theFile >> pExhaustPower >> comma >> pExhaustMass >> 
        comma >> pExhaustTemperature){

        ExhaustPower[i] = pExhaustPower;
        ExhaustMass[i] = pExhaustMass;
        ExhaustTemperature[i] = pExhaustTemperature;

        i += 1;
    }

    theFile.close();

    /* ======= Stability Conditions ======= */

    double max_HP = *std::max_element(Fy_HP.begin(), Fy_HP.end()) + 
        *std::max_element(Fx_HP.begin(), Fx_HP.end());

    double max_cool = *std::max_element(Fy_cool.begin(), Fy_cool.end()) + 
        *std::max_element(Fx_cool.begin(), Fx_cool.end());

    std::cout << "Stabillity coef, Fy.max() + Fx.max() < 0.5 = " << 
        std::max(max_HP,max_cool) << std::endl << std::endl;

    if  (std::max(max_HP,max_cool)> 0.5){
        std::cout << "WARNING not stable, proced? [1]" << std::endl;
        std::cin >> a;
        if (a!=1){
            std::cout << "Stability conditions not met! " << std::endl;
            exit(0);
        }
    }

    /* ======= Cooler matrix  ======= */
    int j = 0;
    for (int i = 0; i < Nx + 1; i++){
        int p = m(i,j); A_cool[p][p] = 1;
    }
    // Loop over all internal mesh points in y diretion
    // and all mesh points in x direction
    for (int j = 1; j < Ncool; j++){
        i = 0;  int p = m(i,j);  A_cool[p][p] = 1; // boundary

        for (int i = 1; i < Nx; i++){
            p = m(i,j);
            A_cool[p][m(i,j-1)] = - Fy_cool[j];
            A_cool[p][m(i-1,j)] = - Fx_cool[j];
            A_cool[p][p]        = 1 + 2*(Fx_cool[j]+Fy_cool[j]);
            A_cool[p][m(i+1,j)] = - Fx_cool[j];
            A_cool[p][m(i,j+1)] = - Fy_cool[j];
        }
        i = Nx;  p = m(i,j);  A_cool[p][p] = 1;  // boundary
    }
    // Equations corresponding to j=Ny, i=0,1,... (u known)
    j = Ncool;
    for (int i = 0; i < Nx + 1; i++){
        int p = m(i,j);  A_cool[p][p] = 1;
    }

    /* ======= HP matrix  ======= */
    j = 0;
    for (int i = 0; i < Nx + 1; i++){
        int p = m(i,j); A_HP[p][p] = 1;
    }
    // Loop over all internal mesh points in y diretion
    // and all mesh points in x direction
    for (int j = 1; j < N_HP; j++){
        i = 0;  int p = m(i,j);  A_HP[p][p] = 1; // boundary

        for (int i = 1; i < Nx; i++){
            p = m(i,j);
            A_HP[p][m(i,j-1)] = - Fy_HP[j];
            A_HP[p][m(i-1,j)] = - Fx_HP[j];
            A_HP[p][p]        = 1 + 2*(Fx_HP[j]+Fy_HP[j]);
            A_HP[p][m(i+1,j)] = - Fx_HP[j];
            A_HP[p][m(i,j+1)] = - Fy_HP[j];
        }
        i = Nx;  p = m(i,j);  A_HP[p][p] = 1;  // boundary
    }
    // Equations corresponding to j=Ny, i=0,1,... (u known)
    j = N_HP;
    for (int i = 0; i < Nx + 1; i++){
        int p = m(i,j);  A_HP[p][p] = 1;
    }

    std::cout << std::endl << "Nx = " << Nx << " N_HP = " <<
        N_HP << " Ncool = " << Ncool << 
        " Nx_TEG = " << Nx_TEG << std::endl << std::endl;

        std::cout << std::endl << "(Nx+1)*(Ncool+1) = " 
            << (Nx+1)*(Ncool+1)<< std::endl;

        std::cout << std::endl << "(Nx+1)*(N_HP+1) = " 
            << (Nx+1)*(N_HP+1)<< std::endl; 


    std::cout << "========================================================================="
        "============" << std::endl;
    std::cout << "========================================================================="
        "============";
    std::cout << std::endl << std::endl << "solving" << std::endl;

    double m_air;
    double sumExcess;
    double Q_reject;
    double uoo_out;
    double ucc_out;
    int theta = 1;
    double Q_qout;
    double Q_in;
    double R_cont = 1/(10000*L_lat*dx);
    double Avoid = voidPercentage_TEG*(Area_TEG_N + Area_TEG_P);

    double I_current_old = 0;
    double e_HP_error = 0;
    double e_cool_error = 0;
    double e_N_error = 0;
    double e_P_error = 0;
    double e_I = 0;
    double Erro_tot = 0;
    double Error = 0;
    double e1 = 0;
    double e2 = 0;
    double e3 = 0;

    double air_visc;
    double air_density;
    double k_air;
    double Cp_air;

    double water_visc;
    double water_density;
    double k_water;
    double Cp_water;

    double write = write_increment;

    Eigen::VectorXd bv_HP((Nx+1)*(N_HP+1));
    Eigen::VectorXd xv_HP((Nx+1)*(N_HP+1));
    Eigen::SparseMatrix<double> Av_HP((Nx+1)*(N_HP+1), (Nx+1)*(N_HP+1));

    Eigen::VectorXd bv_cool((Nx+1)*(Ncool+1));
    Eigen::VectorXd xv_cool((Nx+1)*(Ncool+1));
    Eigen::SparseMatrix<double> Av_cool((Nx+1)*(Ncool+1), (Nx+1)*(Ncool+1));
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_B;

    std::vector<std::vector<double>> u_HP_old(Nx+1, 
        std::vector<double> (N_HP+1, initTemperatureField));
    std::vector<std::vector<double>> u_cool_old(Nx+1, 
        std::vector<double> (Ncool+1, initTemperatureField));
    std::vector<std::vector<double>> u_TEG_N_old(Nx+1, 
        std::vector<double> (Nx_TEG+1, initTemperatureField));
    std::vector<std::vector<double>> u_TEG_P_old(Nx+1, 
        std::vector<double> (Nx_TEG+1, initTemperatureField));


    std::vector<std::vector<double>> e_HP(Nx+1, 
        std::vector<double> (N_HP+1));
    std::vector<std::vector<double>> e_cool(Nx+1, 
        std::vector<double> (Ncool+1));
    std::vector<std::vector<double>> e_N(Nx+1, 
        std::vector<double> (Nx_TEG+1));
    std::vector<std::vector<double>> e_P(Nx+1, 
        std::vector<double> (Nx_TEG+1));


    for (int t = 0; t < 1 + ExhaustTemperature.size()/dt; t++){


        double error = 1;
        int it = 0;

        // cycle information
        if (t*dt > ExhaustTemperature.size()){ 
            uoo[0] = ExhaustTemperature[int(ceil(t*dt))];
            m_air = (1e-3)*ExhaustMass[int(ceil(t*dt))]/2;
        } else{ 
            uoo[0] = (int(ceil(t*dt)) - t*dt)*ExhaustTemperature[int(ceil(t*dt))] + 
                (1-(int(ceil(t*dt)) - t*dt))*ExhaustTemperature[int(ceil(t*dt))+1];

            m_air = (1e-3)*( (int(ceil(t*dt)) - t*dt)*ExhaustMass[int(ceil(t*dt))] + 
                (1-(int(ceil(t*dt)) - t*dt))*ExhaustMass[int(ceil(t*dt))+1] )/2;
        }

        /* ======= heat spreading  ======= */
        for (int i; i < Nx + 1; i++){

            sumExcess = std::accumulate(Q_excess.begin(), Q_excess.end(), 0.0f);

            if (Q_excess[i] == 0){
                Q_HP[i][N_HP] = sumExcess;
                
            } else {
                if (sumExcess > 0 & i == Nx + 1){
                    Q_reject = sumExcess;
                }
            }
        }

        Error = 1;

        while ( Error > error_tol){

            /*for ( int i = N_slices_cooler-1; i < Nx + 1; i = i + N_slices_cooler){

                ucc[i] = Twater_in; // cooler inlet
            }*/
            ucc[Nx] = Twater_in;

            for (int i = 0; i < Nx; i++){
                
                air_visc = prop_air(uoo[i])[0];
                //air_density = prop_air(uoo[i])[1];
                k_air = prop_air(uoo[i])[2];
                Cp_air = prop_air(uoo[i])[3];

                // ======= HE conv =======
                R_air[i] = r_Exhaust(air_visc,k_air,
                    Cp_air,m_air);
                uoo[i+1] = ((uoo[i]-u_HP[i][0])/R_air[i] - m_air*Cp_air*uoo[i])/
                    (- m_air*Cp_air);
                Qoo[i] = (uoo[i]-uoo[i+1])*m_air*Cp_air; 
                uLoo[i] = u_HP_n[i][0] + (Qoo[i]*dt)/(density_HP*dx*dy*L_lat*cp_HP) - 
                    ( ( (L_lat*dx)*k_HP*(u_HP[i][0] - u_HP[i][1])/dy ) *dt ) / 
                    (density_HP*dx*dy*L_lat*cp_HP);

                
                // cooler conv
                water_visc = water_properties(ucc[Nx-i])[0];
                water_density = water_properties(ucc[Nx-i])[1];
                k_water = water_properties(ucc[Nx-i])[2];
                Cp_water = water_properties(ucc[Nx-i])[3];

                R_water[Nx-i] = r_cooler_water(water_density,
                    k_water, k_cool, Cp_water, water_visc);


                bool condition = (i != (ceil(i/N_slices_cooler))*N_slices_cooler);
                if ((condition == 1) || (i == 0)){

                    R_water[0] = r_cooler_water(water_density,
                        k_water, k_cool, Cp_water, water_visc);

                    ucc[Nx-1-i] = ((u_cool[Nx-i][Ncool] - ucc[Nx-i])/R_water[Nx-i] + 
                        2*m_cooler*Cp_water*ucc[Nx-i])/
                        (2*m_cooler*Cp_water);
                    Qcc[Nx-i] = (ucc[Nx-1-i]-ucc[Nx-i])*2*m_cooler*Cp_water;
                    uLcc[Nx-i] = u_cool_n[Nx-i][Ncool] - (Qcc[Nx-i]*dt)/(water_density*
                        dx*dy*L_lat*Cp_water) + 
                        ( ( (L_lat*dx)*k_cool*(u_cool[Nx-i][Ncool-1] - u_cool[Nx-i][Ncool])/dy ) *dt ) / 
                        (water_density*dx*dy*L_lat*Cp_water);

                } else {
                    ucc[Nx-1-i] = Twater_in;

                water_visc = water_properties(ucc[Nx-i])[0];
                water_density = water_properties(ucc[Nx-i])[1];
                k_water = water_properties(ucc[Nx-i])[2];
                Cp_water = water_properties(ucc[Nx-i])[3];


                R_water[(Nx-i)] = r_cooler_water(water_density, k_water,k_cool, Cp_water,
                    water_visc);
                ucc_out = ((u_cool[Nx-i][Ncool] - ucc[Nx-i])/R_water[Nx-i] + 
                    2*m_cooler*Cp_water*ucc[Nx-i])/
                    (2*m_cooler*Cp_water);
                Qcc[Nx-i] = (ucc_out-ucc[Nx-i])*2*m_cooler*Cp_water;
                uLcc[Nx-i] = u_cool_n[Nx-i][Ncool] - 
                    (Qcc[Nx-i]*dt)/(density_cool*dx*dy*L_lat*cp_cool) + 
                    ( ( (L_lat*dx)*k_cool*(u_cool[Nx-i][Ncool-1] - u_cool[Nx-i][Ncool])/dy ) *dt) / 
                    (density_cool*dx*dy*L_lat*cp_cool);    


                }
                //std::cout << "cond=" << ((condition == 1) || (i == 0)) << " Nx-1-i=" << Nx-1-i << " " << uLcc[Nx-i] << std::endl;
            }
                //getchar();
            // ======= outflows =======
            i = Nx;

            air_visc = prop_air(uoo[i-1])[0];
            //air_density = prop_air(uoo[i-1])[1];
            k_air = prop_air(uoo[i-1])[2];
            Cp_air = prop_air(uoo[i-1])[3];


            R_air[i] = r_Exhaust( air_visc,k_air,
                    Cp_air,m_air );
            uoo_out = ((uoo[i]-u_HP[i][0])/R_air[i] - m_air*Cp_air*uoo[i])/
                (- m_air*Cp_air);
            Qoo[i] = (uoo[i]-uoo_out)*m_air*Cp_air;
            uLoo[i] = u_HP_n[i][0] + (Qoo[i]*dt)/(density_HP*dx*dy*L_lat*cp_HP) - 
                ( ( (L_lat*dx)*k_HP*(u_HP[i][0] - u_HP[i][1])/dy ) *dt ) / 
                (density_HP*dx*dy*L_lat*cp_HP);

            water_visc = water_properties(ucc[Nx-i])[0];
            water_density = water_properties(ucc[Nx-i])[1];
            k_water = water_properties(ucc[Nx-i])[2];
            Cp_water = water_properties(ucc[Nx-i])[3];

            R_water[Nx-i] = r_cooler_water(water_density, 
                k_water,k_cool, Cp_water,
                water_visc);
            ucc_out = ((u_cool[Nx-i][Ncool] - ucc[Nx-i])/R_water[Nx-i] + 
                2*m_cooler*Cp_water*ucc[Nx-i])/
                (2*m_cooler*Cp_water);
            Qcc[Nx-i] = (ucc_out-ucc[Nx-i])*2*m_cooler*Cp_water;
            uLcc[Nx-i] = u_cool_n[Nx-i][Ncool] - (Qcc[Nx-i]*dt)/(density_cool*dx*dy*L_lat*cp_cool) + 
                ( ( (L_lat*dx)*k_cool*(u_cool[Nx-i][Ncool-1] - u_cool[Nx-i][Ncool])/dy ) *dt) / 
                (density_cool*dx*dy*L_lat*cp_cool);


            /*for (int i = N_slices_cooler-1; i < number_of_coolers*(N_slices_cooler-1) + 1; i = i + N_slices_cooler-1){

                water_visc = water_properties(ucc[Nx-i])[0];
                water_density = water_properties(ucc[Nx-i])[1];
                k_water = water_properties(ucc[Nx-i])[2];
                Cp_water = water_properties(ucc[Nx-i])[3];


                R_water[(Nx-i)] = r_cooler_water(water_density, k_water,k_cool, Cp_water,
                    water_visc);
                ucc_out = ((u_cool[Nx-i][Ncool] - ucc[Nx-i])/R_water[Nx-i] + 
                    2*m_cooler*Cp_water*ucc[Nx-i])/
                    (2*m_cooler*Cp_water);
                Qcc[Nx-i] = (ucc_out-ucc[Nx-i])*2*m_cooler*Cp_water;
                uLcc[Nx-i] = u_cool_n[Nx-i][Ncool] - 
                    (Qcc[Nx-i]*dt)/(density_cool*dx*dy*L_lat*cp_cool) + 
                    ( ( (L_lat*dx)*k_cool*(u_cool[Nx-i][Ncool-1] - u_cool[Nx-i][Ncool])/dy ) *dt) / 
                    (density_cool*dx*dy*L_lat*cp_cool);    
            }*/
            // ======= End BC convection =======

            if (it > max_iter or std::isnan(error) == 1 or std::isinf(error) == 1){

                std::cout << std::endl << "========================================================================="
                    "============" << std::endl;
                std::cout << std::endl << "Does not converge!" << std::endl;
                std::cout << "========================================================================="
                    "============" << std::endl;
                std::cout << std::endl << "time" << t*dt << "time step: " << t << 
                "[s], Error: " << error << " with " << it << " iterations "  << std::endl;
                exit(0);
            }

            /* ======= HP section  ======= */

            // Compute b_HP
            j = 0;
            for (int i = 0; i < Nx + 1; i++){

                b_HP[m(i,j)] = uLoo[i]; // convective boundary
            }
            for (int j = 1; j < N_HP; j++){
                i = 0;
                b_HP[m(i,j)] = u_HP[i+1][j] + (1-theta)*
                (Fx_HP[j]*(u_HP_n[i+1][j] - u_HP_n[i][j] ) +
                Fy_HP[j]*(u_HP_n[i][j+1] - 2*u_HP_n[i][j] + u_HP_n[i][j-1])) +
                theta*dt*Q_HP[i][j]/C_HP[j] + 
                (1-theta)*dt*Q_HP[i][j]/C_HP[j]; // adiabatic boundary

                for (int i = 1; i < Nx; i ++){ 
                    int p = m(i,j);                    // interior
                    b_HP[p] = u_HP_n[i][j] + 
                    (1-theta)*(Fx_HP[j]*(u_HP_n[i+1][j] - 2*u_HP_n[i][j] + u_HP_n[i-1][j]) +
                    Fy_HP[j]*(u_HP_n[i][j+1] - 2*u_HP_n[i][j] + u_HP_n[i][j-1])) +
                    theta*dt*Q_HP[i][j]/C_HP[j] + 
                    (1-theta)*dt*Q_HP[i][j]/C_HP[j];
                }    
                i = Nx;
                b_HP[m(i,j)] = u_HP[i-1][j] + (1-theta)*(dx*Fx_HP[j]*(u_HP_n[i][j] + u_HP_n[i-1][j]) +
                dy*Fy_HP[j]*(u_HP_n[i][j+1] - 2*u_HP_n[i][j] + u_HP_n[i][j-1])) +
                theta*dt*Q_HP[i][j]/C_HP[j] + 
                (1-theta)*dt*Q_HP[i][j]/C_HP[j]; // adiabatic boundary
            // boundary
            }
            j = N_HP;
            for (int i = 0; i < Nx + 1; i++){ 
                Q_qout = Q1[i];
                Q_in = (L_lat*dx)*k_HP*(u_HP_n[i][N_HP-1]-u_HP_n[i][N_HP])/dy;
                b_HP[m(i,j)] = u_HP_n[i][j] + dt*Q_HP[i][j]/C_HP[j] + 
                    ((Q_in-Q_qout)*dt)/(density_HP*dx*dy*L_lat*cp_HP);
            }

            //SOLVER
            
            for (int i = 0; i < Nx +1; i++){ //temperature control
                if (u_HP[i][N_HP] > T_TEG_MAX){ 
                    b_HP[m(i,j)] = T_TEG_MAX;            
                    u_HP[i][j] = T_TEG_MAX;
                }
            }

            for (int row = 0; row < (Nx+1)*(N_HP+1); ++row){

                bv_HP(row) = b_HP[row];

                for (int col = 0; col < (Nx+1)*(N_HP+1); ++col){
                    Av_HP.coeffRef(row,col) = A_HP[row][col];
                }
            }

            solver.compute(Av_HP);
            xv_HP = solver.solve(bv_HP);

            // Fill u with vector c
            for (int i = 0; i < Nx + 1; i ++){ 
                for (int j = 0; j < N_HP + 1; j++){
                    u_HP[i][j] = xv_HP[m(i,j)];
                }
            }
                    
            // ================================= HP vapor =================================            
            for (int i = 0; i < Nx + 1; i++){
                if (u_HP[i][N_HP] > T_TEG_MAX){
                    
                    Q_excess[i] = ( (u_HP[i][N_HP] - T_TEG_MAX)*density_HP*
                        (dx*dy*L_lat)*cp_HP )/dt;

                    b_HP[m(i,j)] = T_TEG_MAX;           
                    u_HP[i][j] = T_TEG_MAX;
                }
            }

            //================================= Cooler section =================================
            
            // Compute b
            j = 0;
            for (int i = 0; i < Nx + 1; i++){ 
                Q_qout = (L_lat*dx)*k_cool*
                    (u_cool_n[i][0]-u_cool_n[i][1])/dy;
                Q_in = Q2[i];
                b_cool[m(i,j - (0))] = u_cool_n[i][j] + 
                    ((Q_in - Q_qout)*dt)/(density_cool*dx*dy*L_lat*cp_cool);
            }
            
            for (int j = 1; j < Ncool; j ++){ 
                i = 0;
                b_cool[m(i,j - (0))] = u_cool[i+1][j] + 
                    (1-theta)*(Fx_cool[j]*(u_cool_n[i+1][j] - 
                    u_cool_n[i][j] ) +
                Fy_cool[j]*(u_cool_n[i][j+1] - 2*u_cool_n[i][j] + 
                u_cool_n[i][j-1])) +
                theta*dt*Q_cool[i][j]/C_cool[j] + 
                (1-theta)*dt*Q_cool[i][j]/C_cool[j]; // adiabatic boundary

                for (int i = 1; i < Nx; i++){ 
                    int p = m(i,j - (0));   // interior
                    b_cool[p] = u_cool_n[i][j] + 
                    (1-theta)*(Fx_cool[j]*(u_cool_n[i+1][j] - 
                    2*u_cool_n[i][j] + u_cool_n[i-1][j]) +
                    Fy_cool[j]*(u_cool_n[i][j+1] - 2*u_cool_n[i][j] + 
                    u_cool_n[i][j-1])) +
                    theta*dt*Q_cool[i][j]/C_cool[j] + 
                    (1-theta)*dt*Q_cool[i][j]/C_cool[j];
                }
                    
                i = Nx;
                b_cool[m(i,j - (0))] = u_cool[i-1][j]; + 
                (1-theta)*(dx*Fx_cool[j]*(u_cool_n[i][j] + u_cool_n[i-1][j]) +
                dy*Fy_cool[j]*(u_cool_n[i][j+1] - 2*u_cool_n[i][j] + 
                u_cool_n[i][j-1])) +
                theta*dt*Q_cool[i][j]/C_cool[j] + 
                (1-theta)*dt*Q_cool[i][j]/C_cool[j]; // adiabatic boundary
            }    
            // boundary
            j = Ncool;
            for (int i = 0; i < Nx + 1; i++){
                b_cool[m(i,j - (0))] = uLcc[i]; // convective boundary

            }

            //SOLVER 
            for (int row = 0; row < (Nx+1)*(Ncool+1); ++row){

                bv_cool(row) = b_cool[row];

                for (int col = 0; col < (Nx+1)*(Ncool+1); ++col){
                    Av_cool.coeffRef(row,col) = A_cool[row][col];
                }
            }

            solver_B.compute(Av_cool);
            xv_cool = solver_B.solve(bv_cool);
            // Fill u with vector c
            for (int i = 0; i < (Nx+1); i++){
                for (int j = 0; j < (Ncool+1); j++){

                    u_cool[i][j] = xv_cool[m(i,j-(0))];
                }
            }

            //================================= TEG section =================================
            for (int i = 0; i < Nx + 1; i++){ 
                
                Q1[i] = ( u_HP[i][N_HP] - (u_TEG_P_main[i][0] + u_TEG_N_main[i][0])/2)/R_cont;
                Q2[i] = ((u_TEG_N_main[i][Nx_TEG] + u_TEG_P_main[i][Nx_TEG])/2 - u_cool[i][0])/R_cont;
                    
                //n - type
                u_TEG_P_main[i] = funTEG(Q1[i], Q2[i], I_current, u_TEG_n_P[i], u_TEG_P_main[i], "p");
                u_TEG_N_main[i] = funTEG(Q1[i], Q2[i], I_current, u_TEG_n_N[i], u_TEG_N_main[i], "n");

                V = std::accumulate(Vy.begin() + N_cer, Vy.end() - N_cer, 0.0f);
                R = std::accumulate(Ry.begin() + N_cer, Ry.end() - N_cer, 0.0f);

                Rtot[i] = ( R * L_lat * dx ) / ( Avoid + Area_TEG_N + Area_TEG_P );
                U[i] = ( ( V * L_lat * dx ) / ( Avoid + Area_TEG_N + Area_TEG_P ) )/2; 
                // vemf - > OCV for one pair
            }

            //std::cout << std::endl;
            
            I_current = (std::accumulate(U.begin(), U.end(), 0.0f))/
            (2*std::accumulate(Rtot.begin(), Rtot.end(), 0.0f)); 
            // 2*Rtot -> Rload = Rinterno

            for (int i = 0; i < Nx; i++){ 
                U_matchedLoad[i] = U[i]-I_current*Rtot[i];
            }

            /* ======= Error calculation ======= */
            e_HP_error = 0;
            e_cool_error = 0;
            e_N_error = 0;
            e_P_error = 0;

            for (int i = 0; i < Nx + 1; i++){
                for (int j = 0; j < N_HP + 1; j ++){
                    e_HP[i][j] = abs(u_HP[i][j] - u_HP_old[i][j]);
                    if (std::to_string(e_HP[i][j]) > std::to_string(e_HP_error)){ 
                        e_HP_error = e_HP[i][j];
                    }
                }
                for (int j = 0; j < Ncool + 1; j ++){
                    e_cool[i][j] = abs(u_cool[i][j] - u_cool_old[i][j]);
                    if (e_cool[i][j] > e_cool_error){ 
                        e_cool_error = e_cool[i][j];
                    }
                }
                for (int j = 0; j < Nx_TEG + 1; j ++){
                    e_N[i][j] = abs(u_TEG_N_main[i][j] - u_TEG_N_old[i][j]);
                    if (e_N[i][j] > e_N_error){ 
                        e_N_error = e_N[i][j];
                    }
                }
                for (int j = 0; j < Nx_TEG + 1; j ++){
                    e_P[i][j] = abs(u_TEG_P_main[i][j] - u_TEG_P_old[i][j]);
                    if (e_P[i][j] > e_P_error){ 
                        e_P_error = e_P[i][j];
                    }
                }
            }

            e_I = abs(I_current_old - I_current);

            e1 = std::max(e_HP_error, e_cool_error);
            e2 = std::max(e_N_error, e_P_error);
            e3 = std::max(e1, e_I);
            Error = std::max(e3, e2);


            Erro_tot = e_HP_error + e_cool_error + e_N_error + e_P_error + e_I;

            I_current_old = I_current;

            for (int i = 0; i < Nx +1; i++){
                for (int j = 0; j < N_HP +1; j++){
                    u_HP_old[i][j] = u_HP[i][j];
                }
            }

            for (int i = 0; i < Nx +1; i++){
                for (int j = 0; j < Ncool +1; j++){
                    u_cool_old[i][j] = u_cool[i][j];
                }
            }

            for (int i = 0; i < Nx +1; i++){
                for (int j = 0; j < Nx_TEG +1; j++){
                    u_TEG_N_old[i][j] = u_TEG_N_main[i][j];
                }
            }
            for (int i = 0; i < Nx +1; i++){
                for (int j = 0; j < Nx_TEG +1; j++){
                    u_TEG_P_old[i][j] = u_TEG_P_main[i][j];
                }
            }

            it += 1;
        }

        /* ======= write cycle ======= */
        if ( ( std::to_string(t*dt) == std::to_string(write) ) || ( t == (ExhaustTemperature.size()/dt) ) ){ 
            

            ETA(start,t*dt,ExhaustTemperature.size(),Error,it, write_increment);        

            std::cout << t*dt << " Current: " << I_current << " U_matchedLoad: " << 
                U_matchedLoad[0] << " U0: " << U[0] << " Rtot: " << Rtot[0] << 
                " Pelect[W]: " << 49*U_matchedLoad[0]*I_current << " Thot: " <<
                u_TEG_N_main[0][0] - 273 << " Tcold: " << u_TEG_N_main[0][Nx_TEG] - 273
                << std::endl; 

            writing(t, outTimeStamp);
                            
            write = write + write_increment;
        }


        for (int i = 0; i < Nx +1; i++){
            for (int j = 0; j < N_HP +1; j++){
                u_HP_n[i][j] = u_HP[i][j];
            }
        }

        for (int i = 0; i < Nx +1; i++){
            for (int j = 0; j < Ncool +1; j++){
                u_cool_n[i][j] = u_cool[i][j];
            }
        }

        for (int i = 0; i < Nx +1; i++){
            for (int j = 0; j < Nx_TEG +1; j++){
                u_TEG_n_N[i][j] = u_TEG_N_main[i][j];
            }
        }
        for (int i = 0; i < Nx +1; i++){
            for (int j = 0; j < Nx_TEG +1; j++){
                u_TEG_n_P[i][j] = u_TEG_P_main[i][j];
            }
        }

    }

    std::cout << std::endl << "Solver ran successfully \n" << "End \n";
}

