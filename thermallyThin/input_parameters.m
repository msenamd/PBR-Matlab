function [T_g, u_g, G, x_O2_g, delta_i, FMC, ...
          rho_m, rho_vs, rho_c, k_m, k_vs, k_c, c_m, c_vs, c_c, ...
          eps_m, eps_vs, eps_c, A_R1, Ta_R1, DeltaH_R1, ...
          A_R2, Ta_R2, DeltaH_R2, eta_c, A_R3, Ta_R3, DeltaH_R3,...
          temp_surf_i,T_end,geometry,A_rectangle,L_cylinder] ...
          = input_parameters

% Particle type (available types: rectangle - cylinder - sphere)
geometry="rectangle"; 

% Conditions of the particle
delta_i     = 1e-3;    % Initial half-thickness or radius of particle [m]
% NOTE: delta_i = Initial half-thickness for rectagular particles 
% or intial radius for cylindrical and spherical particles
FMC         = 0.0;     % Initial fuel moisture content [-]
temp_surf_i = 300;     % Initial surface temperature [K]
A_rectangle = 1e-0;    % Area of exposed surface (for rectangular particles only) [m2]
L_cylinder  = 1e-3;    % Cylinder length (for cylindrical particles only) [m]

% Duration of simulation [s]
T_end = 500;              

% Local gas conditions
T_g    = 700;     % Local ambient gas temperature [K]
u_g    = 1;       % Local gas velocity [m/s]
G      = 0.0;   % Averaged irradiation [W/m2]
x_O2_g = 0;       % Local oxygen mole fraction in ambient gas [-]

%  - Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
rho_m  = 1000;    % Mass density of liquid water [kg/m3]
rho_vs = 663;     % Mass density of virgin solid [kg/m3]
rho_c  = 132.6;   % Mass density of solid char [kg/m3]

k_m    = 0.6;     % Conductivity of liquid water [W/m/K]
k_vs   = 0.126;   % Conductivity of virgin solid [W/m/K]
k_c    = 0.126;   % Conductivity of solid char [W/m/K]

c_m    = 4.18e+3;   % Heat capacity of liquid water [J/kg/K]
c_vs   = 2.52e+3;   % Heat capacity of virgin solid [J/kg/K]
c_c    = 2.52e+3;   % Heat capacity of solid char [J/kg/K]

eps_m  = 0.9;   % Surface emissivity of liquid water [-]
eps_vs = 0.9;   % Surface emissivity of virgin solid [-]
eps_c  = 0.9;   % Surface emissivity of solid char [-]

% Ideal gas constant
R = 8.3145;   % [J/K/mol]

% Moisture evaporation model (reaction R1)
%  - Source: Shen et al. (2007) Fire Safety J. 42:210-217
A_R1      = 5.13e+10;  % Pre-exponential factor [1/s]
E_R1      = 88e+3;     % Activation energy [J/mol]
Ta_R1     = E_R1/R;    % Activation temperature [K]
DeltaH_R1 =-2440e+3;   % Heat of evaporation [J/kg] (< 0, endothermic)

% Pyrolysis model (reaction R2)
%  - Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
A_R2      = 5.25e+7;    % Pre-exponential factor [1/s]
E_R2      = 1.256e+5;   % Activation energy [J/mol]
Ta_R2     = E_R2/R;     % Activation temperature [K]
DeltaH_R2 =-0;          % Heat of pyrolysis [J/kg] (< 0, endothermic)
eta_c     = 0.0;        % Char yield [-]

% Char oxidation model (reaction R3)
A_R3      = 0;        % Pre-exponential factor [1/s]
E_R3      = 0;        % Activation energy [J/mol]
Ta_R3     = E_R3/R;   % Activation temperature [K]
DeltaH_R3 =+0;        % Heat of char oxidation [J/kg] (> 0, exothermic)


end 