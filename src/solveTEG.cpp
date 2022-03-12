

std::vector<double> funTEG(double Q1TEG, double Q2TEG, double I_current,std::vector<double> u_TEG_n,
	std::vector<double> u_TEG, std::string n_or_p){
	
	
		/* TEG leg */
	for (int y = N_cer + N_connector; y < N_cer + N_connector +N_legs; y++){

		std::vector<double> PropTEG = select_propTEG(n_or_p, u_TEG[y]);

		alfa_TEG[y] = PropTEG[0];
		k_TEG[y] = PropTEG[2];
		cp_TEG[y] = prop_TEG_Cp(u_TEG[y]);
		rho_TEG[y] = PropTEG[1];
		density_TEG[y] = rhoDense_leg_TEG;
	}

	
    double Qin = k_TEG[0]*select_areaTEG(n_or_p)*(-u_TEG[Nx_TEG] + u_TEG[Nx_TEG-1])
		/delta_x;

	double TQ1 = (Q1TEG*dt)/(density_TEG[0]*VsliceTEG*cp_TEG[0]);
	double TQ2 = ((Qin - Q2TEG) *dt)/(density_TEG[Nx_TEG-1]*VsliceTEG*cp_TEG[Nx_TEG-1]);

	/* Matrix Assembly */

	for (int y = 1; y < Nx_TEG; y++){

		double E_T_Teg = (select_propTEGtau(n_or_p) * I_current * dt) / 
		(cp_TEG[y] * density_TEG[y] * select_areaTEG(n_or_p) * delta_x);

        double E_Teg = (rho_TEG[y]*pow(I_current,2)*dt)/
		(cp_TEG[y]*density_TEG[y]*pow(select_areaTEG(n_or_p),2));

		double F_Teg = ( (k_TEG[y]/(density_TEG[y]*cp_TEG[y]))*dt )/pow(delta_x,2);

        A_TEG[y][y-1] = - E_T_Teg - F_Teg;
        A_TEG[y][y+1] = E_T_Teg - F_Teg;
        A_TEG[y][y] = 1 + 2*F_Teg;
        
        b_TEG[y] = u_TEG_n[y] + E_Teg;

        
	}

	double E_Teg1 = (select_propTEGtau(n_or_p) * I_current * dt) / 
		(cp_TEG[0] * density_TEG[0] * select_areaTEG(n_or_p) * delta_x);
    A_TEG[0][0] = 1;
    A_TEG[Nx_TEG][Nx_TEG] = 1;
    b_TEG[0] = u_TEG_n[0] + E_Teg1 + TQ1;

	double E_Teg2 = (select_propTEGtau(n_or_p) * I_current * dt) / 
		(cp_TEG[Nx_TEG-1] * density_TEG[Nx_TEG-1] * select_areaTEG(n_or_p) * delta_x);
    b_TEG[Nx_TEG] = u_TEG_n[Nx_TEG] + E_Teg2 + TQ2;


	/* ============ TEG Solver ============ */

    Eigen::VectorXd b(Nx_TEG +1);
    Eigen::VectorXd x(Nx_TEG +1);

    Eigen::SparseMatrix<double> A(Nx_TEG +1, Nx_TEG +1);

    for (int row = 0; row < Nx_TEG +1; ++row){

        b(row) = b_TEG[row];

        for (int col = 0; col < Nx_TEG +1; ++col){
            A.coeffRef(row,col) = A_TEG[row][col];
        }
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);

	for (int y = 1; y < Nx_TEG + 1; y++){

        Vy[y-1] = std::abs( alfa_TEG[y-1]*(u_TEG[y-1] - u_TEG[y]) ) + 
			std::abs( alfa_TEG[y-1]*(u_TEG[y-1] - u_TEG[y]) );

            //std::cout << Vy[y-1] << " ";

        Ry[y-1] = (rho_TEG[y-1]*delta_x)/select_areaTEG(n_or_p) + 
			(rho_TEG[y-1]*delta_x)/select_areaTEG(n_or_p);


            //std::cout << Ry[y-1] << " ";

	}

    //std::cout << std::endl;
	
    std::vector<double> Out(Nx_TEG+1);

    for (int row = 0; row < Nx_TEG +1; ++row){
        Out[row] = x(row);
    }

    return Out;
}

int m( int i, int j){

    return j*(Nx+1) + i;

}

