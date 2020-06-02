
#include <iostream>
#include <chrono>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>

#include <conio.h> 
#include <io.h> 
#include <process.h> 
#include <stdio.h> 

#include <C:\dev\THERMOSS\HP2Dt\cpp\src\eigen\eigen\Sparse>
#include <C:\dev\THERMOSS\HP2Dt\cpp\src\eigen\eigen\Dense>
#include <C:\dev\THERMOSS\HP2Dt\cpp\src\eigen\Eigen\IterativeLinearSolvers>

auto start = std::chrono::system_clock::now();
auto timeElapsed1_old = std::chrono::system_clock::now();



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



std::string fileHE = "HEDim";
std::string fileMat = "matProp";
std::string fileCool = "coolerDim";
std::string fileNum = "numVal";



/*======= COOLER Geometry =======*/

double m_cooler = getVar(fileCool, "m_cooler"); 
//[kg/s] cooling water mass flow rate per plate
double Nc = getVar(fileCool, "Nc"); //number of channels
double channel_widt = getVar(fileCool, "channel_widt"); //channel width [m]
double fin_widt = getVar(fileCool, "fin_widt"); // fin widt [m]
double fin_heigt = getVar(fileCool, "fin_heigt");  // fin heigt [m]
double channel_length = getVar(fileCool, "channel_length"); // channel length [m]
double lateral_length = getVar(fileCool, "lateral_length"); 
// lateral length [m] comprimento transversal
double Metal_tick = getVar(fileCool, "Metal_tick"); 
// [m] metal thickness joined by the fins*/

double k_cool  = getVar(fileMat, "k_cool"); 
double cp_cool  = getVar(fileMat, "cp_cool"); 
double density_cool  = getVar(fileMat, "density_cool"); 

double number_of_coolers = getVar(fileCool, "number_of_coolers"); 

double Twater_in = getVar(fileCool, "Twater_in"); 


/*======= HE Geometry =======*/

/*geometry definition*/

double s = getVar(fileHE, "s"); // [m] fin spacing
double h = getVar(fileHE, "h"); // [m] fin heigth
double tick = getVar(fileHE, "tick"); // [m] fin width
double l = getVar(fileHE, "l"); // [m] offset length
double number_of_channels = getVar(fileHE, "number_of_channels"); 
// [-] number of channels 
double number_of_offsets = getVar(fileHE, "number_of_offsets"); 
// [-] number of offsets
double k_HE = getVar(fileMat, "k_HE"); // [W/mK] thermal conductivity for the offset fins
//20 W/mK average value for the thermal conductivity of the stainless steel
double L_lat = number_of_channels*(s+tick);


//estas areas são para o PC todo, para os calculos das alhetas,ao
// dividir o PC a meio, é necessário dividir por 2!
double total_area = 2*(number_of_channels*number_of_offsets)*(s*l+h*l);
//[m2] total heat transfer area 

//hidraulic diameter for offset fins

double alfa_geo = s/h;
double sigma_geo = tick/l;
double gama_geo = tick/s;

double D_hidraulic = (4*s*h*l)/(2*(s*l+h*l+tick*h)+tick*s); 
// [m] hydrauslic dyameter

double k_HP = getVar(fileMat, "k_HP"); 
double cp_HP = getVar(fileMat, "cp_HP"); 
double density_HE = getVar(fileMat, "density_HE"); 

double density_HP = getVar(fileMat, "density_HP"); 

/*======= TEG Geometry =======*/

double select = getVar("matProp", "selectTEG"); // 1 - Hi-z data

double leg_length = getVar("TEGDim", "leg_length"); 
// [m] TEG length in y direction
double L_Al_connector = getVar("TEGDim", "L_Al_connector"); 
// [m] connector thickkness
double Lcer = getVar("TEGDim", "Lcer"); 
// [m] ceramic thickness
double delta_x = getVar("numVal", "delta_x");

double leg_side_N = getVar("TEGDim", "leg_side_N");
double leg_side_P = getVar("TEGDim", "leg_side_P");


//double L_total = leg_length + 2*L_Al_connector + 2*Lcer;

int N_legs = round(leg_length/delta_x); // number of slices in the legs
int N_connector = round(L_Al_connector/delta_x);
int N_cer = round(Lcer/delta_x);
// number of slices in each connector
int Nx_TEG = N_legs + 2*N_connector + 2*N_cer; // total number of slices


double rhoDense_leg_TEG = getVar("matProp", "rhoDense_leg_TEG"); // kg/m3

double k_cer_TEG = getVar("matProp", "k_cer_TEG");
double cp_cer_TEG = getVar("matProp", "cp_cer_TEG");
double density_cer_TEG = getVar("matProp", "density_cer_TEG");


double Tamb = getVar("matProp", "Tamb");

double Area_TEG_N = leg_side_N*leg_side_N;
double Area_TEG_P = leg_side_P*leg_side_P;

std::vector<std::vector<double>> A_TEG(
	(Nx_TEG + 1), std::vector<double> (Nx_TEG + 1));
std::vector<double> b_TEG(Nx_TEG + 1);
std::vector<double> alfa_TEG(Nx_TEG);
std::vector<double> k_TEG(Nx_TEG);
std::vector<double> cp_TEG(Nx_TEG);
std::vector<double> rho_TEG(Nx_TEG); // electrical resistivity [ohm.m]
std::vector<double> density_TEG(Nx_TEG); // [kg/m3]

std::vector<double> Ry(Nx_TEG);
std::vector<double> Vy(Nx_TEG);

double Tau_P = getVar("matProp", "Tau_P"); // [V/K] thomson coefficient
double Tau_N = getVar("matProp", "Tau_N"); // [V/K] thomson coefficient

double k_connector_TEG = getVar("matProp", "k_connector_TEG");
double cp_connector_TEG = getVar("matProp", "cp_connector_TEG");
double density_connector_TEG = getVar("matProp", "density_connector_TEG");
double rho_connector_TEG = getVar("matProp", "rho_connector_TEG");
double alfa_connector_TEG = getVar("matProp", "alfa_connector_TEG");

double V;
double R;

double I_current = 0;

double voidPercentage_TEG = getVar("TEGDim", "voidPercentage_TEG");

/*======= Numerical Parameters =======*/

double dt = getVar(fileNum, "dt");
double dx = getVar(fileNum, "dx");
double dy = getVar(fileNum, "dy");
double write_increment = getVar(fileNum, "write_increment");

int Nx = round(l*number_of_offsets/dx); // int Nx 
double L_HP = getVar("HEDim", "L_HP");
int N_HP = round(L_HP/dy);
int Ncool = round(Metal_tick/dy);

float initTemperatureField = getVar(fileNum, "initTemperatureField");

int N_slices_cooler = round(Nx / number_of_coolers);
double N_slices = round(l*number_of_offsets/dx); 
//1+HE_length/dx // number of slices 
double VsliceTEG = delta_x*dx*L_lat; // volume of a TEG slice

float error_tol = getVar(fileNum, "error_tol");
float max_iter = getVar(fileNum, "max_iter");

/*======= Simulation conditions =======*/

double T_TEG_MAX = getVar(fileHE, "T_TEG_MAX");

std::vector<double> U(Nx+1);

std::vector<std::vector<double>> A_HP(
	(Nx+1)*(N_HP+1), std::vector<double> ((Nx+1)*(N_HP+1)));
std::vector<double> b_HP((Nx+1)*(N_HP+1));

std::vector<std::vector<double>> A_cool(
	(Nx+1)*(Ncool+1), std::vector<double> ((Nx+1)*(Ncool+1), 0.0));
std::vector<double> b_cool((Nx+1)*(Ncool+1), 0.0);

// u init
std::vector<std::vector<double>> u_HP(Nx+1, 
    std::vector<double> (N_HP+1, initTemperatureField));
std::vector<std::vector<double>> u_cool(Nx+1, 
    std::vector<double> (Ncool+1, initTemperatureField));
std::vector<std::vector<double>> u_HP_n(Nx+1, 
    std::vector<double> (N_HP+1, initTemperatureField));
std::vector<std::vector<double>> u_cool_n(Nx+1, 
    std::vector<double> (Ncool+1, initTemperatureField));

std::vector<std::vector<double>> u_TEG_P_main(Nx+1, 
    std::vector<double> (Nx_TEG+1, initTemperatureField));
std::vector<std::vector<double>> u_TEG_N_main(Nx+1, 
    std::vector<double> (Nx_TEG+1, initTemperatureField));
std::vector<std::vector<double>> u_TEG_n_P(Nx+1, 
    std::vector<double> (Nx_TEG+1, initTemperatureField));
std::vector<std::vector<double>> u_TEG_n_N(Nx+1, 
    std::vector<double> (Nx_TEG+1, initTemperatureField));

    

std::vector<double> alpha_HP(N_HP+1, k_HP/(cp_HP * density_HP));
std::vector<double> C_HP(N_HP+1, cp_HP * density_HP);
std::vector<double> Fy_HP(N_HP+1, 
    ( (k_HP/(cp_HP * density_HP))*dt )/pow(dy,2) );
std::vector<double> Fx_HP(N_HP+1, 
    ( (k_HP/(cp_HP * density_HP))*dt )/pow(dx,2) );


std::vector<double> alpha_cool(Ncool+1, k_cool/(cp_cool * density_cool));
std::vector<double> C_cool(Ncool+1, cp_cool * density_cool);
std::vector<double> Fy_cool(Ncool+1, 
    ( (k_cool/(cp_cool * density_cool))*dt )/pow(dy,2) );
std::vector<double> Fx_cool(Ncool+1, 
    ( (k_cool/(cp_cool * density_cool))*dt )/pow(dx,2) );

std::vector<double> uoo((Nx+1));
std::vector<double> ucc((Nx+1));
std::vector<double> Qoo((Nx+1));
std::vector<double> Qcc((Nx+1));
std::vector<double> uLoo((Nx+1));
std::vector<double> uLcc((Nx+1));
std::vector<double> R_air((Nx+1));
std::vector<double> R_water((Nx+1));
std::vector<double> Rtot((Nx+1));
std::vector<double> Q_excess((Nx+1));
std::vector<double> Q1 ((Nx+1));
std::vector<double> Q2 ((Nx+1));
std::vector<double> U_matchedLoad ((Nx+1));

std::vector<std::vector<double>> Q_HP(Nx+1, 
    std::vector<double> (N_HP+1, 0));
std::vector<std::vector<double>> Q_cool(Nx+1, 
    std::vector<double> (Ncool+1, 0));


std::vector<double> properties (4);
std::vector<double> Prop_N (3);
std::vector<double> Prop_P (3);
std::vector<double> properties_air(4);

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


void disp(double print){
	std::cout << std::endl << print << std::endl;
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
                
                std::cout << Q_HP[i][N_HP] << " ";

            } else {
                if (sumExcess > 0 & i == Nx + 1){
                    Q_reject = sumExcess;
                }
            }
        }

        Error = 1;

        while ( Error > error_tol){


            for ( int i = 0; i < Nx + 1; i = i + N_slices_cooler){

                ucc[i] = Twater_in; // cooler inlet
            }
            
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

                if ( i != N_slices_cooler-1 & i != (2*N_slices_cooler-1) & 
                    i != (3*N_slices_cooler-1) & i != (4*N_slices_cooler-1) & 
                    i != (5*N_slices_cooler-1) ){

                    ucc[Nx-1-i] = ((u_cool[Nx-i][Ncool] - ucc[Nx-i])/R_water[Nx-i] + 
                        2*m_cooler*Cp_water*ucc[Nx-i])/
                        (2*m_cooler*Cp_water);
                    Qcc[Nx-i] = (ucc[Nx-1-i]-ucc[Nx-i])*2*m_cooler*Cp_water;
                    uLcc[Nx-i] = u_cool_n[Nx-i][Ncool] - (Qcc[Nx-i]*dt)/(water_density*
                        dx*dy*L_lat*Cp_water) + 
                        ( ( (L_lat*dx)*k_cool*(u_cool[Nx-i][Ncool-1] - u_cool[Nx-i][Ncool])/dy ) *dt ) / 
                        (water_density*dx*dy*L_lat*Cp_water);
                }
            }

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


            for (int i = N_slices_cooler-1; i < number_of_coolers*(N_slices_cooler-1) + 1; i = i + N_slices_cooler-1){

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


int main() {

    init_n_print();
	initPropTEG();
    HP_2Dtsolve();


	std::cin.get();

    return 0;
}
