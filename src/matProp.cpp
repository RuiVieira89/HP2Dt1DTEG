


std::vector<double> water_properties(double T){

	// [Pa.s] LIQUID water viscosity (ZGGRAFOS,1987)
	double visc_liq = 3.8208*pow(10,-2)*pow( (T-252.33), -1); 

	// [kg/m3] LIQUID water volumic mass 1atm (ZGGRAFOS,1987)
	double rho_liq = -3.0115*pow(10,-6)*pow(T,3)+9.6272*pow(10,-4)*pow(T,2)
	-0.11052*T+1022.4;

	// [W/mK] water thermal conductivity
	double k_water = 4.2365*pow(10,-9)*pow(T,3)-1.144*pow(10,-5)*pow(T,2)
	+7.1959*pow(10,-3)*T-0.63262;

	// [J/kgK] water cp (Zografos, 1986)
	double Cp_water = 1000*(1.785*pow(10,-7)*pow(T,3)-1.9149*pow(10,-4)*pow(T,2)
		+6.7953*pow(10,-2)*T-3.7559);


	properties[0] = visc_liq; 
	properties[1] = rho_liq;
	properties[2] = k_water;
	properties[3] = Cp_water;

	return properties;

}




double r_cooler_water(double rho_water, double k_water, double k_Al, double Cp_water, 
	double visc_water){
	
    double fin_area_coller = 2 * fin_heigt * channel_length * Nc;  
	// [m2] Af finned area
    double fin_perimeter_coller = 2 * (fin_widt + channel_length);  
	// [m] fin tip perimeter
    double non_fin_area = Nc * fin_widt * channel_length;  
	// [m2] Ab area of the base with no fin
    double total_HT_area = fin_area_coller + non_fin_area;  
	// [m2] At=Ab+Af total heat tranfer area
    double flow_area = fin_widt * fin_heigt;  
	// [m2] single channel flow area
    double Dh_cooler = (4 * flow_area) / (2 * (fin_widt + fin_heigt));
	// [m] hidraulic diameter 4A/P
    double D_charact_cooler = pow(flow_area,0.5);
	// [m] characteristic lenght A**1/2 

	/*Cooler convection model constants, check tese Rui Vieira 2017*/

	/*FLOW Conditions*/
    double Pr_water = (Cp_water *visc_water )/k_water;

    double C1=3.24; // Uniform heat power
    double C2=1.5;
    double C3=0.409; // Uniform heat power
    double C4=2;
    double Par_forma=0.1; // Par de forma (>=90)
    double ee=fin_widt/fin_heigt; // channel aspect ratio e=a/b of a rectangule

    double fPr =0.564/( pow( (1+ pow( (1.664* pow(Pr_water, (1/6))), 9/2) ), 2/9)); 
	// Uniform heat power
    double m_coef_cool =2.27+1.65*pow(Pr_water,(1/3)); // m coefficient 

    // [m3/s] cooling water volumic flow per channel
    double Q_water_channel =(m_cooler/rho_water )/Nc;

    // [m/s] water velocity in the channels
    double u_water =Q_water_channel /flow_area;


    //water Reynolds number
    double Re_water =(rho_water *u_water *D_charact_cooler)/visc_water;

    double L_coef =0.058*D_charact_cooler*Re_water;  // [m] L' entry region
    double Z_plus =pow((L_coef /channel_length),-1)/(Re_water ); // Z+ coefficient
    double Z_star =((channel_length/D_charact_cooler)/(Re_water *Pr_water ));
	// Z* coefficient

    // fRe(A**0.5) coefficient
    double fRe = pow(( pow((12/((pow(ee,0.5))*(1+ee)*(1-((192*ee)/(pow(M_PI,5)))*
	tanh(M_PI/(2*ee))))),2)+pow( (3.44/pow(Z_plus, 0.5)), 2)), 0.5);




    //Nusselt number
    double Nu_water = pow( pow( (C4*fPr )/pow(Z_star, 0.5), m_coef_cool) + 
	( pow( (( pow((C2*C3* pow(fRe /Z_star, 1/3)), 5)) + 
	pow( (C1*(fRe /(8*pow(M_PI,0.5)*pow(ee,Par_forma) ))),5) ), m_coef_cool /5)), 
	1/m_coef_cool );
    // [Wm2/K] heat transfer coefficient
    double h_water =(Nu_water *k_water )/D_charact_cooler;

    // m fin efficiency coeficient
    double m_coef_water = pow((h_water *fin_perimeter_coller)/(k_Al*non_fin_area), 0.5);
    // fin efficiency
    double fin_eff_water =(tanh(m_coef_water *fin_heigt))/(m_coef_water *fin_heigt);
    // group fin efficiency
    double fin_eff_water_group =1-(fin_area_coller/total_HT_area)*(1-fin_eff_water );


    //Cooler thermal resistance [K/W]
    double R_cooler = (1/(h_water *(total_HT_area/(channel_length/dx))*fin_eff_water_group ))/2;

	return R_cooler;
}



std::vector<double> prop_air(double T){

    //AIR PROPERTIES IN S.I. (Zografos, 1986) (equações na tese Rui Vieira, 2017 
	//no ultimo anexo)
    //Calculates the air resistance for Th 

    //viscosity as a function of hot inlet temperature [Pa.s]
    double air_viscosity = ((2.5914*pow(10,-15))*pow(T,3) ) - 
	(1.4346*pow(10,-11)*pow(T,2)) + 
	((5.0523*pow(10,-8))*T)+4.113*pow(10,-6); //[Pa.s]
            
    //density as a function of hot inlet temperature [kg/m3]
    double air_density = 101325/((8314.4598/28.97)*T); //[kg/m3]
            
    //air conductivity [W/mK]
    double air_conductivity = ((0.000000000015207)*pow(T,3) ) 
	- ((0.000000048574)*pow(T,2) )+( (0.00010184)*T)- 0.00039333;
            
    //Specific Heat [J/kg.K]
    double air_specific_heat=(((1.3864*pow(10,-13))*pow(T,4)) 
	- ((6.4747*pow(10,-10))*pow(T,3))
	+ ((1.0234*pow(10,-6))*pow(T,2))-(4.3282*pow(10,-4))*T+1.0613)*1000;
    
    properties_air[0] = air_viscosity;
	properties_air[1] = air_density;
	properties_air[2] = air_conductivity;
	properties_air[3] = air_specific_heat;
    
    return properties_air;
}
    


double r_Exhaust(double air_viscosity,double air_conductivity,double air_specific_heat,
	double m_air){
	
    //Dynamic and themal parameters 
    //Reynolds Number
    double Re=(4/M_PI)*(((m_air/number_of_channels))/(air_viscosity*D_hidraulic));

    //Colburn factor (j)
    double j=(0.6522*pow(Re,-0.5403)*pow(alfa_geo,-0.1541)*pow(sigma_geo,0.1499)*
	pow(gama_geo,-0.0678))*pow( (1+(5.269*pow(10,-5))*pow(Re,1.34)*pow(alfa_geo,0.504)*
	pow(sigma_geo,0.456)*pow(gama_geo,-1.055)),0.1);

    //Prandtl number
    double Pr=(air_specific_heat*air_viscosity)/air_conductivity;

    //Nusselt number
    double Nu=j*Re*pow(Pr,1/3);

    //Heat tranfer coefficient [W/m2K]
    double h_air=(air_conductivity*Nu)/D_hidraulic;


    //Fin efficiency coefficient

    double m_fin = pow((h_air*2*(tick+l))/(k_HE*tick*l),0.5);

    //Fin efficiency 
    double fin_eff=(tanh(m_fin*h/2))/(m_fin*h/2);

    //Fin group efficiency  
    double fin_eff_group = 1-((2*h/2)/(2*h/2+s))*(1-fin_eff);

    //Exhaust heat exchanger thermal resistance [K/W]
    double R_air = (1/(((total_area/2)/N_slices)*h_air*fin_eff_group));
    
    return R_air;
}



double prop_TEG_Cp(double T){

    return (36.858*pow(10,-3)*T + 117.02 - 161744*pow(T,-2) )/(800.761*pow(10,-3));
}
    


std::vector<double> prop_TEG_N(double T){


    // all temperatures must be in [K]

	if(select == 1){

        //// Hi-z data
        
        // Seebeck coeficient [V/K]
        double alfa_Bi2Te3N = 0.00007423215 - 0.0000015018 * T + 
		0.0000000029361 * pow(T,2) - 0.000000000002499 * pow(T,3) + 
		0.000000000000001361 * pow(T,4); // [V/K]
        
        // electrical resistivity [ohm.m]
        double rho_Bi2Te3N = -0.00195922 + 0.00001791526 * T - 
		0.00000003818 * pow(T,2) + 0.000000000049186 * pow(T,3) - 
		0.0000000000000298 * pow(T,4);// ohm.cm
        rho_Bi2Te3N = rho_Bi2Te3N/100; // ohm.m
        
        // thermal conductivity [W/mK]
        double k_Bi2Te3N = 1.425785 + 0.006514882 * T - 0.00005162 * pow(T,2) + 
		0.00000011246 * pow(T,3) - 0.000000000076 * pow(T,4); // W/mK

		Prop_N[0] = alfa_Bi2Te3N;
		Prop_N[1] = rho_Bi2Te3N;
		Prop_N[2] = k_Bi2Te3N;
		
		return Prop_N;

	} else{

        // Cyprus data
        
        double alfa_Bi2Te3N = -(3.8273*pow(10,-11)*pow(T,5) - 
		9.3841*pow(10,-8)*pow(T,4) + 9.1048*pow(10,-5)*pow(T,3) - 
		4.3783*pow(10,-2)*pow(T,2) + 1.0608*pow(10,+1)*T - 
		8.7621*pow(10,2))*pow(10,-6); //R2 = 9.981710**-01
        
        double rho_Bi2Te3N = 1.3292*pow(10,-15)*pow(T,4) - 2.8350*pow(10,-12)*pow(T,3) +
		 2.3015*pow(10,-9)*pow(T,2) - 8.3722*pow(10,-7)*T + 1.7220*pow(10,-4);
		 //R2 = 9.9562E-01
        
        double k_Bi2Te3N = -7.3633*pow(10,-16)*(T,6) + 2.4093*pow(10,-12)*pow(T,5) - 
		3.1602*pow(10,-9)*pow(T,4) + 2.1370*pow(10,-6)*pow(T,3) - 
		7.8786*pow(10,-4)*pow(T,2) + 1.5046*pow(10,-1)*T - 1.1021*pow(10,+1);
		//R2 = 9.9508E-01

		Prop_N[0] = alfa_Bi2Te3N;
		Prop_N[1] = rho_Bi2Te3N;
		Prop_N[2] = k_Bi2Te3N;
		
		return Prop_N;

	}        

}



std::vector<double> prop_TEG_P(double T){

    // all temperatures must be in [K]

	if(select == 1){

        //// Hi-z data
        // Seebeck coeficient [V/K]
        double alfa_Bi2Te3P = -0.0002559037 + 0.0000023184 * T - 
		0.000000003181 * pow(T,2) + 0.0000000000009173 * pow(T,3) - 
		0.000000000000000488 * pow(T,4); // [V/K]
        
        // electrical resistivity [ohm.m]
        double rho_Bi2Te3P = -0.002849603 + 0.00001967684 * T - 0.00000003317 * pow(T,2) + 
		0.000000000034733 * pow(T,3) - 0.000000000000019 * pow(T,4); // ohm.cm
        rho_Bi2Te3P = rho_Bi2Te3P/100; // ohm.m
        
        // thermal conductivity [W/mK]
        double k_Bi2Te3P = 6.9245746 - 0.05118914 * T + 0.000199588 * pow(T,2) - 
		0.0000003891 * pow(T,3) + 0.00000000030382 * pow(T,4); // W/mK

		Prop_P[0] = alfa_Bi2Te3P;
		Prop_P[1] = rho_Bi2Te3P;
		Prop_P[2] = k_Bi2Te3P;
		
		return Prop_P;

	} else{

        // Cyprus data
        
        double alfa_Bi2Te3P = (3.8273*pow(10,-11)*pow(T,5) - 9.3841*pow(10,-8)*pow(T,4) + 
		9.1048*pow(10,-5)*pow(T,3) - 4.3783*pow(10,-2)*pow(T,2) + 1.0608*pow(10,+1)*T - 
		8.7621*pow(10,+2))*pow(10,-6); //R² = 9.981710,-01
        
        double rho_Bi2Te3P = 1.3292*pow(10,-15)*pow(T,4) - 2.8350*pow(10,-12)*pow(T,3) + 
		2.3015*pow(10,-9)*pow(T,2) - 8.3722*pow(10,-7)*T + 1.7220*pow(10,-4);
		 //R² = 9.9562E-01
        
        
        double k_Bi2Te3P = -7.3633*pow(10,-16)*pow(T,6) + 2.4093*pow(10,-12)*pow(T,5) - 
		3.1602*pow(10,-9)*pow(T,4) + 2.1370*pow(10,-6)*pow(T,3) - 
		7.8786*pow(10,-4)*pow(T,2) + 1.5046*pow(10,-1)*T - 1.1021*pow(10,+1);
		//R² = 9.9508E-01

		Prop_P[0] = alfa_Bi2Te3P;
		Prop_P[1] = rho_Bi2Te3P;
		Prop_P[2] = k_Bi2Te3P;
		
		return Prop_P;

	}        

}


double select_areaTEG(std::string n_or_p){

	if (n_or_p == "n"){
		return Area_TEG_N;
	} else {
		return Area_TEG_P;
	}

}

std::vector<double> select_propTEG(std::string n_or_p, double u){

	if (n_or_p == "n"){
		return prop_TEG_N(u);
	} else {
		return prop_TEG_P(u);
	}

}

double select_propTEGtau(std::string n_or_p){

	if (n_or_p == "n"){
		return Tau_N;
	} else {
		return Tau_P;
	}

}



void initPropTEG(){

    double rho_cer = 0.001;

	/* 1st Ceramic part */
	for (int y = 0; y < N_cer; y++){

		alfa_TEG[y] = 0;
		k_TEG[y] = k_cer_TEG;
		cp_TEG[y] = cp_cer_TEG;
		rho_TEG[y] = rho_cer;
		density_TEG[y] = density_cer_TEG;
	}

	/* 1st Connector part */

	for (int y = N_cer; y < N_cer + N_connector; y++){

		alfa_TEG[y] = alfa_connector_TEG;
		k_TEG[y] = k_connector_TEG;
		cp_TEG[y] = cp_connector_TEG;
		rho_TEG[y] = rho_connector_TEG;
		density_TEG[y] = density_connector_TEG;
	}

	/* 2nd Connector part */

	for (int y = N_cer + N_connector +N_legs; 
		y < N_cer + N_connector + N_legs + N_connector; y++){

		alfa_TEG[y] = alfa_connector_TEG;
		k_TEG[y] = k_connector_TEG;
		cp_TEG[y] = cp_connector_TEG;
		rho_TEG[y] = rho_connector_TEG;
		density_TEG[y] = density_connector_TEG;
	}

	/* 2nd Ceramic part */
	for (int y = N_cer + N_connector +N_legs + N_connector; 
		y < Nx_TEG; y++){

		alfa_TEG[y] = 0;
		k_TEG[y] = k_cer_TEG;
		cp_TEG[y] = cp_cer_TEG;
		rho_TEG[y] = rho_cer;
		density_TEG[y] = density_cer_TEG;
	}

}
