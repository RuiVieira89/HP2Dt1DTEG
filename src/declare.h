#ifndef DECLARE
#define DECLARE


auto start = std::chrono::system_clock::now();
auto timeElapsed1_old = std::chrono::system_clock::now();


// initialization 
double getVar(std::string file, std::string idVar);
void ETA(auto start_time, double time_par, double final_time_par, 
	double error_par, int it, double write_increment);
void init_n_print();


// Material Properties
std::vector<double> water_properties(double T);
double r_cooler_water(double rho_water, double k_water, double k_Al, double Cp_water, 
	double visc_water);
std::vector<double> prop_air(double T);  
double r_Exhaust(double air_viscosity,double air_conductivity,double air_specific_heat,
	double m_air);
double prop_TEG_Cp(double T);
std::vector<double> prop_TEG_N(double T);
std::vector<double> prop_TEG_P(double T);
double select_areaTEG(std::string n_or_p);
std::vector<double> select_propTEG(std::string n_or_p, double u);
double select_propTEGtau(std::string n_or_p);
void initPropTEG();


// write to csv
void writing(double t, std::string outTimeStamp);


// Solvers
std::vector<double> funTEG(double Q1TEG, double Q2TEG, double I_current,std::vector<double> u_TEG_n,
	std::vector<double> u_TEG, std::string n_or_p);
void HP_2Dtsolve();


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
std::vector<double> Q_excess((Nx+1), 0);
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

#endif
