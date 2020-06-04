# HP2Dt1DTEG
A 2D transient heat transfer model to predict the electrical power output of a TEG system on a car. This would model the waste heat recovery potential of a Heat Pipe assisted TEG system.

It employs Eigen, a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms. Check their website: http://eigen.tuxfamily.org/index.php?title=Main_Page
The Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > method was used to solve the equations.

In the input_data folder there are all the variables that are relevant for the simulation. 
In the cycle folder there's the driving cycle used for this simulation, a WLTP template is provided. The engine conditions necessary to fulfil the driving cycle are predicted using a steady-state engine model map, from which the exhaust mass flow and temperature are extracted. The file cycle.csv is an example of the thermal output of a driving cycle.

In the folder input_data\controlVars there are the main variables the can be modified.

# coolerDim - Cooler variables
A finned cooler was simulated
m_cooler - Water cooler mass flow [kg/s]
Nc - Number of flow channels in the cooler
channel_widt - Channel widt [m]
fin_widt - fin widt [m]
fin_heigt fin heigt [m]
channel_length - channel_length [m]
lateral_length - transversal length [m]
Metal_tick - cooler thickness from the base of the fins to the TEG module [m]
number_of_coolers - number of coolers [m]
Twater_in - inlet water temperature [K]

# HEDim - Heat exchanger dimensions
An offset fin Heat exchanger was simulated in this example.
s - fin spacing [m]
h - fin height [m]
tick - fin thickness [m] 
l - offset length [m] 
number_of_channels - number of channels [m]
number_of_offsets - number of offsets [m]
L_HP - Heat Pipe length [m]
T_TEG_MAX - Maximum TEG temperature [K]

# matProp - material properties
selectTEG 1
k_HE 
k_HP 
cp_HP 
density_HP 
density_HE 
k_cool 
cp_cool 
density_cool 
Tamb 
rhoDense_leg_TEG 
Tau_P 
Tau_N 
alfa_connector_TEG 
k_connector_TEG 
cp_connector_TEG 
density_connector_TEG 
rho_connector_TEG
k_cer_TEG 
cp_cer_TEG 
density_cer_TEG 

# numVal - Numerical variables
dx - x discretization [m]
dy - y discretization [m]
dt - time discretization [s]
delta_x - TEG discretization [m]
write_increment - time step to write results [s]
initTemperatureField - initial temperature of the system [K]
error_tol - maximum error allowed
max_iter - maximmum number of iterations

# TEGDim - TEG dimensions
n_leg - number of TEG legs
leg_side_N - N type TEG side length
leg_side_P  P type TEG side length
leg_length - leg length
L_Al_connector connector thickness
Lcer ceramic thickness
voidPercentage_TEG - empty space percentage between modules
