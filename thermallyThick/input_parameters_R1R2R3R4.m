function [ T_end, geometry, delta_i, A_rectangle, L_cylinder, ...
           T_g, u_g, G, Y_g_O2, pres_g, k_g0, cp_g0, nu_g0, Pr, MW_g, ...          
           temp_i, x_ws_i, x_ds_i, x_c_i, x_a_i, Y_O2_i, pres_i, ...
           rho_ws, rho_ds, rho_c, rho_a, ...
           k_ws, k_ds, k_c, k_a, ...
           c_ws, c_ds, c_c,c_a, ...
           eps_ws, eps_ds, eps_c, eps_a, ...
           psi_ws, psi_ds, psi_c, psi_a, ...
           Kperm_ws, Kperm_ds, Kperm_c, Kperm_a, ...
           R, sigma, ...
           A_R1, Ta_R1, n_R1, DeltaH_R1, eta_ds_R1, ...
           A_R2, Ta_R2, n_R2, DeltaH_R2, eta_c_R2, ...
           A_R3, Ta_R3, n_R3, n_O2_R3, DeltaH_R3, eta_c_R3, eta_O2_R3, ...
           A_R4, Ta_R4, n_R4, n_O2_R4, DeltaH_R4, eta_a_R4, eta_O2_R4 ] ...
           = input_parameters
       
% Duration of the simulation [s]
T_end = 10000;

% Particle shape (available shapes: rectangle - cylinder - sphere)
geometry = "rectangle";  

% Characteristic length scales of the particle gemoetry
%  - delta_i: Initial half-thickness (geometry = "rectangle") or
%             intial radius (geometry = "cylinder" or "sphere") [m]
%  - A_rectangle: Area of exposed surface (geometry = "rectangle") [m2]
%  - L_cylinder: Axial length of the cylinder (geometry = "cylinder") [m]
delta_i     = 1.25e-2;   
A_rectangle = 1e-0;
L_cylinder  = 1e-0;

% Conditions in the ambient external gas
T_g    = 300;     % Temperature [K]
u_g    = 0.0;     % Flow velocity [m/s]
%%AT G      = 40e+3;   % Averaged irradiation [W/m2]
u_g = 10.0;
G   = 50e+3;
%%AT
Y_g_O2 = 0.233;   % Oxygen mass fraction [-]
pres_g = 101325;  % Absolute pressure [Pa]

k_g0   = 0.026;   % Thermal conductivity in air at normal temp. [W/m/K]
cp_g0  = 1100;    % Heat capacity (at constant pressure) in air
                  % at normal temperature [J/kg/K]
nu_g0  = 1.6d-5;  % Kinematic viscosity in air at normal temp./pres. [m2/s]
Pr     = 0.7;     % Prandtl number [-]
MW_g   = 0.029;   % Molecular weight [kg/mole]

% Initial conditions inside the particle (t = 0)
temp_i = 300;     % Initial temperature [K] (solid/gas phases)
x_ws_i = 1.0;     % Initial volume fraction of wet solid [-] (solid phase)
x_ds_i = 0.0;     % Initial volume fraction of dry solid [-] (solid phase)
x_c_i  = 0.0;     % Initial volume fraction of char [-] (solid phase)
x_a_i  = 0.0;     % Initial volume fraction of ash [-] (solid phase)
Y_O2_i = Y_g_O2;  % Initial mass fraction of oxygen [-] (gas phase)
pres_i = 0;       % Initial gauge pressure [Pa] (gas phase)

%  - Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
rho_ws = 390;     % Mass density of wet solid [kg/m3]
                  % NB: at constant volume and porosity in reaction (R1),
                  %     rho_ws = rho_ds*(1+FMC) or FMC = (rho_ws/rho_ds)-1
rho_ds = 390;     % Mass density of dry solid [kg/m3]
rho_c  = 390;     % Mass density of char [kg/m3]
                  % NB: at constant volume and porosity in reaction (R2),
                  %     eta_c_R2 = (rho_c/rho_ds)
rho_a  = 390;     % Mass density of ash [kg/m3]

k_ws   = 0.186;   % Conductivity of wet solid [W/m/K]
k_ds   = 0.176;   % Conductivity of dry solid [W/m/K]
k_c    = 0.065;   % Conductivity of char [W/m/K]
k_a    = 0.058;   % Conductivity of ash [W/m/K]

c_ws   = 1764;    % Heat capacity of wet solid [J/kg/K]
c_ds   = 1664;    % Heat capacity of dry solid [J/kg/K]
c_c    = 1219;    % Heat capacity of char [J/kg/K]
c_a    = 1244;    % Heat capacity of ash [J/kg/K]

eps_ws = 0.757;   % Surface emissivity of wet solid [-]
eps_ds = 0.759;   % Surface emissivity of dry solid [-]
eps_c  = 0.957;   % Surface emissivity of char [-]
eps_a  = 0.955;   % Surface emissivity of ash [-]

psi_ws = 1-(380/390); % Porosity of the particle when pure wet solid [-]
psi_ds = 1-(361/390); % Porosity of the particle when pure solid [-]
psi_c  = 1- (73/390); % Porosity of the particle when pure char [-]
psi_a  = 1-(5.7/390); % Porosity of the particle when pure ash [-]

Kperm_ws = 1e-10;   % Permeability of the particle when pure wet solid [m2]
Kperm_ds = 1e-10;   % Permeability of the particle when pure solid [m2]
Kperm_c  = 1e-10;   % Permeability of the particle when pure char [m2]
Kperm_a  = 1e-10;   % Permeability of the particle when pure ash [m2]

% Miscellaneous
R     = 8.3145;    % Ideal gas constant [J/K/mol]
sigma = 5.67e-8;   % Stefan-Boltzmann constant [W/m2/K4]

% Moisture evaporation model (reaction R1)
%  - Source: Shen et al. (2007) Fire Safety J. 42:210-217
A_R1      = 4.29e+3;    % Pre-exponential factor [1/s]
E_R1      = 43.8e+3;    % Activation energy [J/mol]
Ta_R1     = E_R1/R;     % Activation temperature [K]
n_R1      = 1;          % Solid-phase reactant exponent [-]
DeltaH_R1 =-2410e+3;    % Heat of evaporation [J/kg] (< 0, endothermic)
eta_ds_R1 = (361/380);  % Mass yield of dry solid in reaction (R1) [-]
                        % NB: at constant volume/porosity in reaction (R1),
                        %     eta_ds_R1 = 1/(1+FMC) = (rho_ds/rho_ws)

% Thermal pyrolysis model (reaction R2)
%  - Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
A_R2      = 3.29e+9;    % Pre-exponential factor [1/s]
E_R2      = 135e+3;     % Activation energy [J/mol]
Ta_R2     = E_R2/R;     % Activation temperature [K]
n_R2      = 1;          % Solid-phase reactant exponent [-]
DeltaH_R2 =-533e+3;     % Heat of pyrolysis [J/kg] (< 0, endothermic)
eta_c_R2  = (73/361);   % Mass yield of char in reaction (R2) [-]
                        % NB: at constant volume and porosity in (R2),
                        %     eta_c_R2 = (rho_c/rho_ds)

% Oxidative pyrolysis model (reaction R3)
%  - Source: Lautenberger & Fernandez-Pello (2009) Combust. Flame
%            156:1503-1513
A_R3      = 6.00e+9;    % Pre-exponential factor [1/s]
E_R3      = 124.2e+3;   % Activation energy [J/mol]
Ta_R3     = E_R3/R;     % Activation temperature [K]
n_R3      = 1;          % Solid-phase reactant exponent [-]
n_O2_R3   = 1;          % Gas-phase oxygen exponent [-]
DeltaH_R3 =+994e+3;     % Heat of pyrolysis [J/kg] (> 0, exothermic)
eta_c_R3  = (73/361);   % Mass yield of char in reaction (R3) [-]
eta_O2_R3 = 0.1*(1-eta_c_R3); % Oxygen-to-dry-solid mass ratio in (R3) [-]

% Char oxidation model (reaction R4)
%  - Source: Lautenberger & Fernandez-Pello (2009) Combust. Flame
%            156:1503-1513
A_R4      = 9.79e+13;   % Pre-exponential factor [1/s]
E_R4      = 192.4e+3;   % Activation energy [J/mol]
Ta_R4     = E_R4/R;     % Activation temperature [K]
n_R4      = 1;          % Solid-phase reactant exponent [-]
n_O2_R4   = 1;          % Gas-phase oxygen exponent [-]
DeltaH_R4 =+37700e+3;   % Heat of pyrolysis [J/kg] (> 0, exothermic)
eta_a_R4  = (5.7/73);   % Mass yield of ash in reaction (R4) [-]
eta_O2_R4 = 2.0*(1-eta_a_R4); % Oxygen-to-char mass ratio in (R4) [-]

%%AT Uncomment the lines below for tests with fixed volume dV
%%{
rho_ws    = 380;
rho_ds    = 361;
rho_c     = 73;
rho_a     = 5.7;
psi_ws    = 0.01;
psi_ds    = 0.01;
psi_c     = 0.01;
psi_a     = 0.01;
eta_ds_R1 = (rho_ds/rho_ws);
eta_c_R2  = (rho_c/rho_ds);
eta_c_R3  = (rho_c/rho_ds);
eta_O2_R3 = 0.1*(1-eta_c_R3);
eta_a_R4  = (rho_a/rho_c);
eta_O2_R4 = 2.0*(1-eta_a_R4);
%}
%%AT


end 