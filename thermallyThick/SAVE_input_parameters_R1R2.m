function [ T_end, geometry, delta_i, A_rectangle, L_cylinder, ...
           T_g, u_g, G, Y_g_O2, pres_g, k_g0, cp_g0, nu_g0, Pr, MW_g, ...          
           temp_i, x_ws_i, x_ds_i, x_c_i, x_a_i, Y_O2_i, pres_i, ...
           rho_ws, rho_ds, rho_c, rho_a, ...
           k0_ws, k0_ds, k0_c, k0_a, nk_ws, nk_ds, nk_c, nk_a, ...
           gamma_ws, gamma_ds, gamma_c, gamma_a, ...
           c0_ws, c0_ds, c0_c, c0_a, nc_ws, nc_ds, nc_c, nc_a, ...
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
T_end = 600;

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
u_g    = 1.0;     % Flow velocity [m/s]
G      = 40e+3;   % Averaged irradiation [W/m2]
x_g_O2 = 0.21;    % Oxygen mole fraction [-]
Y_g_O2 = x_g_O2*32/(x_g_O2*32+(1-x_g_O2)*28); % Oxygen mass fraction [-]
pres_g = 101325;  % Absolute pressure [Pa]

k_g0   = 0.026;   % Thermal conductivity in air at normal temp. [W/m/K]
cp_g0  = 1000;    % Heat capacity (at constant pressure) in air
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
rho_ws = 729.3;   % Mass density of wet solid [kg/m3]
                  % NB: at constant volume and porosity in reaction (R1),
                  %     rho_ws = rho_ds*(1+FMC) or FMC = (rho_ws/rho_ds)-1
rho_ds = 663;     % Mass density of dry solid [kg/m3]
rho_c  = 132.6;   % Mass density of char [kg/m3]
                  % NB: at constant volume and porosity in reaction (R2),
                  %     eta_c_R2 = (rho_c/rho_ds)
rho_a  = 132.6;   % Mass density of ash [kg/m3]

k0_ws  = 0.126;   % Conductivity of wet solid [W/m/K]
k0_ds  = 0.126;   % Conductivity of dry solid [W/m/K]
k0_c   = 0.126;   % Conductivity of char [W/m/K]
k0_a   = 0.126;   % Conductivity of ash [W/m/K]
nk_ws  = 0.0;     % Temperature exponent of conductivity of wet solid
nk_ds  = 0.0;     % Temperature exponent of conductivity of dry solid
nk_c   = 0.0;     % Temperature exponent of conductivity of char
nk_a   = 0.0;     % Temperature exponent of conductivity of ash

gamma_ws = 0;      % Effective radiation conductivity of wet solid [m]
gamma_ds = 0;      % Effective radiation conductivity of dry solid [m]
gamma_c  = 0;      % Effective radiation conductivity of char [m]
gamma_a  = 0;      % Effective radiation conductivity of ash [m]

c0_ws  = 2.52e+3; % Heat capacity of wet solid [J/kg/K]
c0_ds  = 2.52e+3; % Heat capacity of dry solid [J/kg/K]
c0_c   = 2.52e+3; % Heat capacity of char [J/kg/K]
c0_a   = 2.52e+3; % Heat capacity of ash [J/kg/K]
nc_ws  = 0.0;     % Temperature exponent of heat capacity of wet solid
nc_ds  = 0.0;     % Temperature exponent of heat capacity of dry solid
nc_c   = 0.0;     % Temperature exponent of heat capacity of char
nc_a   = 0.0;     % Temperature exponent of heat capacity of ash

eps_ws = 0.9;     % Surface emissivity of wet solid [-]
eps_ds = 0.9;     % Surface emissivity of dry solid [-]
eps_c  = 0.9;     % Surface emissivity of char [-]
eps_a  = 0.9;     % Surface emissivity of ash [-]

psi_ws = 0.001;   % Porosity of the particle when pure wet solid [-]
psi_ds = 0.001;   % Porosity of the particle when pure solid [-]
psi_c  = 0.001;   % Porosity of the particle when pure char [-]
psi_a  = 0.001;   % Porosity of the particle when pure ash [-]

Kperm_ws = 1e-10;   % Permeability of the particle when pure wet solid [m2]
Kperm_ds = 1e-10;   % Permeability of the particle when pure solid [m2]
Kperm_c  = 1e-10;   % Permeability of the particle when pure char [m2]
Kperm_a  = 1e-10;   % Permeability of the particle when pure ash [m2]

% Miscellaneous
R     = 8.3145;    % Ideal gas constant [J/K/mol]
sigma = 5.67e-8;   % Stefan-Boltzmann constant [W/m2/K4]

% Moisture evaporation model (reaction R1)
%  - Source: Shen et al. (2007) Fire Safety J. 42:210-217
A_R1      = 5.13e+10;   % Pre-exponential factor [1/s]
E_R1      = 88e+3;      % Activation energy [J/mol]
Ta_R1     = E_R1/R;     % Activation temperature [K]
n_R1      = 1;          % Solid-phase reactant exponent [-]
DeltaH_R1 =-2440e+3;    % Heat of evaporation [J/kg] (< 0, endothermic)
eta_ds_R1 = 0.90909091; % Mass yield of dry solid in reaction (R1) [-]
                        % NB: at constant volume/porosity in reaction (R1),
                        %     eta_ds_R1 = 1/(1+FMC) = (rho_ds/rho_ws)

% Thermal pyrolysis model (reaction R2)
%  - Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
A_R2      = 5.25e+7;    % Pre-exponential factor [1/s]
E_R2      = 1.256e+5;   % Activation energy [J/mol]
Ta_R2     = E_R2/R;     % Activation temperature [K]
n_R2      = 1;          % Solid-phase reactant exponent [-]
DeltaH_R2 =-0;          % Heat of pyrolysis [J/kg] (< 0, endothermic)
eta_c_R2  = 0.20;       % Mass yield of char in reaction (R2) [-]
                        % NB: at constant volume and porosity in (R2),
                        %     eta_c_R2 = (rho_c/rho_ds)
% Oxidative pyrolysis model (reaction R3)
%  - Source: Lautenberger & Fernandez-Pello (2009) Combust. Flame
%            156:1503-1513
A_R3      = 0;          % Pre-exponential factor [1/s]
E_R3      = 0;          % Activation energy [J/mol]
Ta_R3     = E_R3/R;     % Activation temperature [K]
n_R3      = 1;          % Solid-phase reactant exponent [-]
n_O2_R3   = 1;          % Gas-phase oxygen exponent [-]
DeltaH_R3 =+0;          % Heat of pyrolysis [J/kg] (< 0, endothermic)
eta_c_R3  = 0;          % Mass yield of char in reaction (R3) [-]
eta_O2_R3 = 0;          % Oxygen-to-dry-solid mass ratio in react. (R3) [-]

% Char oxidation model (reaction R4)
%  - Source: Lautenberger & Fernandez-Pello (2009) Combust. Flame
%            156:1503-1513
A_R4      = 0;          % Pre-exponential factor [1/s]
E_R4      = 0;          % Activation energy [J/mol]
Ta_R4     = E_R4/R;     % Activation temperature [K]
n_R4      = 1;          % Solid-phase reactant exponent [-]
n_O2_R4   = 1;          % Gas-phase oxygen exponent [-]
DeltaH_R4 =+0;          % Heat of pyrolysis [J/kg] (< 0, endothermic)
eta_a_R4  = 0;          % Mass yield of ash in reaction (R4) [-]
eta_O2_R4 = 0;          % Oxygen-to-char mass ratio in reaction (R4) [-]

end 