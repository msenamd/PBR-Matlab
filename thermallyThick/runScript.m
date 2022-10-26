% Particle Burning Rate Model - PBR_model.m

clc;
clear;
close all;

global geometry A_rectangle L_cylinder ...
       T_g u_g G Y_g_O2 pres_g k_g0 cp_g0 nu_g0 MW_g R sigma ...
       rho_ws rho_ds rho_c rho_a ...
       k0_ws k0_ds k0_c k0_a nk_ws nk_ds nk_c nk_a ...
       gamma_ws gamma_ds gamma_c gamma_a ...
       c0_ws c0_ds c0_c c0_a nc_ws nc_ds nc_c nc_a ...
       eps_ws eps_ds eps_c eps_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 DeltaH_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 DeltaH_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 DeltaH_R3 eta_c_R3 eta_O2_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 DeltaH_R4 eta_a_R4 eta_O2_R4 ...
       dx_i IFilter nFilter Tmax_R4


% Read input parameters
[ T_end, geometry, delta_i, A_rectangle, L_cylinder, ...
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
         = input_parameters;

% NOTE: delta_i = Initial half-thickness (geometry = "rectangle") or
%                 initial radius (geometry = "cylinder" or "sphere")
%       Select spatial resolution ~ 100 microns
%AT dx0  = 100e-6;               % Spatial resolution [m]
%AT dx0  = 50e-6;
%AT dx0  = 25e-6;
dx0  = 10e-6;
nx_i = round(delta_i/dx0);   % Initial number of grid cells
dx_i = delta_i/nx_i;              % Initial grid cell size
fprintf(' dx_i = %g \n',dx_i);
fprintf(' nx_i = %g \n',nx_i);
fprintf(' \n');

% Calculate initial time step
dt_i = timeStep(rho_ws, rho_ds, rho_c, rho_a, ...
                k0_ws, k0_ds, k0_c, k0_a, nk_ws, nk_ds, nk_c, nk_a, ...
                gamma_ws, gamma_ds, gamma_c, gamma_a, ...
                c0_ws, c0_ds, c0_c, c0_a, nc_ws, nc_ds, nc_c, nc_a, ...
                psi_ws, psi_ds, psi_c, psi_a, ...
                Kperm_ws, Kperm_ds, Kperm_c, Kperm_a, ...
                A_R1, Ta_R1, ...
                A_R2, Ta_R2, ...
                A_R3, Ta_R3, n_O2_R3, ...
                A_R4, Ta_R4, n_O2_R4, ...
                Y_g_O2, nu_g0, MW_g, R, dx_i);

dt_i = min(dt_i,1e-3);
G_i  = G;

%n_f = round(T_end/dt_i);   % Total number of time steps
n_f = round(1000000);   % Total number of time steps (TBC)
fprintf(' n_f = %g \n',n_f);
fprintf(' \n');
pause;

% Parameter that controls the frequency at which quantities are saved
% to output files
%%AT i_output  = 50;
i_output  = 100;
%%AT i_output  = 500;

% Parameters that control the time step and solution accuracy
%%AT Threshold_Temp = 1.0;    % Max. value of temp variation during dt [K]
%%AT Threshold_Temp = 0.1;    % Max. value of temp variation during dt [K]
Threshold_Temp = 0.1;    % Max. value of temp variation during dt [K]
Threshold_xk   = 0.01;   % Max. value of x_k variation during dt [-]
Threshold_YO2  = 0.01;   % Max. value of Y_O2 variation during dt [-]
Threshold_pres = 0.1;    % Max. value of pres variation during dt [Pa]
lambda         = 0.1;    % Under-relaxation parameter for temp & Y_O2 [-]
lambda_xk      = 0.0;    % Under-relaxation parameter for x_k [-]
lambda_pres    = 0.0;    % Under-relaxation parameter for pressure [-]
%AT dt_max         = 0.1;    % Max. value of dt [s]
dt_max         = 0.01;    % Max. value of dt [s]
dt_max_i       = dt_max;

% Parameters that control filtering of reaction rates for R3 and R4
IFilter = 0; % IFilter =  1 to activate filtering
nFilter = 0; % nFilter >= 1 to activate filtering
%AT
Tmax_R4 = 2000; % Control of RR4
%AT

% Set Initial conditions (t = 0)
n         = 0;
time      = 0; 
n_output  = 0;

nx_new    = nx_i;

temp_surf = temp_i;    % Surface temperature [K] (t = 0)

% Calculation of convective heat transfer coefficient
T_film       = 0.5*(temp_surf+T_g);       % Film temperature [K]
nu_g         = nu_g0*(T_film/300)^(1.76); % Kin. visc. at film temp.
k_g          = k_g0*(T_film/300)^(0.76);  % Conductivity at film temp.
if geometry=="rectangle"
    D_eff    = 2*delta_i; % Effective diameter of particle [m]
    Re_D     = u_g*D_eff/nu_g;
    Nu_D     = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
elseif geometry=="cylinder"
    D_eff    = 2*delta_i;
    Re_D     = u_g*D_eff/nu_g;
    Nu_D     = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
elseif geometry=="sphere"
    D_eff    = 2*delta_i;
    k_g      = k_g0*(T_g/300)^(0.76);        % Conductivity at gas temp.
    nu_g     = nu_g0*(T_g/300)^(1.76);       % Kin. visc. at gas temp.
    nu_s     = nu_g0*(temp_surf/300)^(1.76); % Kin. visc. at surface temp.
    rho_g    = pres_g*MW_g/R/T_g;
    rho_surf = pres_g*MW_g/R/temp_surf;
    mu_g     = rho_g*nu_g;
    mu_s     = rho_surf*nu_s;
    Re_D     = u_g*D_eff/nu_g; % Calculated at gas temperature
    Nu_D     = 2+(0.4*Re_D^0.5+0.06*Re_D^(2/3))*Pr^0.4*(mu_g/mu_s)^0.25;
else
    fprintf(...
    ' Invalid geometry, valid geometries are: cylinder/sphere/rectangle ');
    return;
end
h_conv   = Nu_D*k_g/D_eff;   % Convective heat transfer coef. [W/m2/K]
%%AT
%%AT h_conv   = 10; % Model of Lautenberger & Fernandez-Pello (2009)
%%AT

q_surf_i = eps_ws*G - eps_ws*sigma*temp_surf^4 ...
         + h_conv*(T_g-temp_surf);  % Net surface heat flux [W/m2] (t = 0)
q_surf   = q_surf_i;

% Solid mass density (t = 0)
rho_s_i = rho_ws*x_ws_i + rho_ds*x_ds_i + rho_c*x_c_i + rho_a*x_a_i;

% Porosity (t = 0)
psi_i   = psi_ws*x_ws_i + psi_ds*x_ds_i + psi_c*x_c_i + psi_a*x_a_i;

% Initial mass of wet solid (t = 0) [kg]
if geometry=="rectangle"
    Mass_ws_i = rho_ws*(1-psi_ws)*delta_i*A_rectangle; % NB: no factor 2
elseif geometry=="cylinder"
    Mass_ws_i = rho_ws*(1-psi_ws)*pi*delta_i^2*L_cylinder;
elseif geometry=="sphere"
    Mass_ws_i = rho_ws*(1-psi_ws)*(4/3)*pi*delta_i^3;
end

dt     = dt_i;                 % Time increment
dx     = dx_i   *ones(nx_i,1); % Array containing grid cell sizes
temp   = temp_i *ones(nx_i,1); % Array containing solid temperatures
rho_s  = rho_s_i*ones(nx_i,1); % Array containing solid mass densities
psi_sg = psi_i  *ones(nx_i,1); % Array containing porosities
x_ws   = x_ws_i *ones(nx_i,1); % Array containing vol. frac. of wet solid
x_ds   = x_ds_i *ones(nx_i,1); % Array containing vol. frac. of dry solid
x_c    = x_c_i  *ones(nx_i,1); % Array containing vol. frac. of char
x_a    = x_a_i  *ones(nx_i,1); % Array containing vol. frac. of ash
Y_O2   = Y_O2_i *ones(nx_i,1); % Array containing mass fraction of O2 (gas)
pres   = pres_i *ones(nx_i,1); % Array containing pressure (gas)

% Allocate additional arrays for saved output quantities
time_save       =          zeros(round(n_f/i_output),1); % Time
MLR_save        =          zeros(round(n_f/i_output),1); % Mass loss rate
q_surf_save     = q_surf_i *ones(round(n_f/i_output),1); % Surf. heat flux
temp_surf_save  = temp_i   *ones(round(n_f/i_output),1); % Surf. temp.
temp_back_save  = temp_i   *ones(round(n_f/i_output),1); % Back temp.
delta_save      = delta_i  *ones(round(n_f/i_output),1); % Thickness
dt_save         = dt_i     *ones(round(n_f/i_output),1); % Time increment
nx_save         = nx_i     *ones(round(n_f/i_output),1); % Number of cells
DeltaQ_max_save =          zeros(4,round(n_f/i_output)); % Max. variation
YO2_surf_save   = Y_g_O2   *ones(round(n_f/i_output),1); % Surf. YO2
YO2_nx_save     = Y_g_O2   *ones(round(n_f/i_output),1); % YO2(nx)
RR1_peak_save   =          zeros(round(n_f/i_output),1); % Peak R1-RR value
RR2_peak_save   =          zeros(round(n_f/i_output),1); % Peak R2-RR value
RR3_peak_save   =          zeros(round(n_f/i_output),1); % Peak R3-RR value
RR4_peak_save   =          zeros(round(n_f/i_output),1); % Peak R4-RR value
xRR1_peak_save  =          zeros(round(n_f/i_output),1); % x-loc of peak R1
xRR2_peak_save  =          zeros(round(n_f/i_output),1); % x-loc of peak R2
xRR3_peak_save  =          zeros(round(n_f/i_output),1); % x-loc of peak R3
xRR4_peak_save  =          zeros(round(n_f/i_output),1); % x-loc of peak R4
h_conv_save     = h_conv   *ones(round(n_f/i_output),1); % Heat tran. coef.
%AT
MLR_R1_save     =          zeros(round(n_f/i_output),1); % MLR R1-RR
MLR_R2_save     =          zeros(round(n_f/i_output),1); % MLR R2-RR
MLR_R3_save     =          zeros(round(n_f/i_output),1); % MLR R3-RR
MLR_R4_save     =          zeros(round(n_f/i_output),1); % MLR R4-RR
temp_max_save   = temp_i   *ones(round(n_f/i_output),1); % Max. TEMP
YO2_back_save   = Y_g_O2   *ones(round(n_f/i_output),1); % Back YO2
YO2_min_save    = Y_g_O2   *ones(round(n_f/i_output),1); % Min. YO2
pres_surf_save  = pres_i   *ones(round(n_f/i_output),1); % Surf. pressure
pres_back_save  = pres_i   *ones(round(n_f/i_output),1); % Back pressure
pres_max_save   = pres_i   *ones(round(n_f/i_output),1); % Max. pressure
MLR_surf_save   =          zeros(round(n_f/i_output),1); % Surf. MLRPUA
vel_surf_save   =          zeros(round(n_f/i_output),1); % Surf. velocity
vel_back_save   =          zeros(round(n_f/i_output),1); % Back velocity
vel_max_save    =          zeros(round(n_f/i_output),1); % Max. velocity
temp_xRR1_peak_save  = temp_i*ones(round(n_f/i_output),1); %T @peak R1 loc
temp_xRR2_peak_save  = temp_i*ones(round(n_f/i_output),1); %T @peak R2 loc
temp_xRR3_peak_save  = temp_i*ones(round(n_f/i_output),1); %T @peak R3 loc
temp_xRR4_peak_save  = temp_i*ones(round(n_f/i_output),1); %T @peak R4 loc
RR1_surf_save   =          zeros(round(n_f/i_output),1); % Surf. R1-RR
RR1_back_save   =          zeros(round(n_f/i_output),1); % Back R1-RR
RR2_surf_save   =          zeros(round(n_f/i_output),1); % Surf. R2-RR
RR2_back_save   =          zeros(round(n_f/i_output),1); % Back R2-RR
RR3_surf_save   =          zeros(round(n_f/i_output),1); % Surf. R3-RR
RR3_back_save   =          zeros(round(n_f/i_output),1); % Back R3-RR
RR4_surf_save   =          zeros(round(n_f/i_output),1); % Surf. R4-RR
RR4_back_save   =          zeros(round(n_f/i_output),1); % Back R4-RR
mass_ws_save    = Mass_ws_i*ones(round(n_f/i_output),1); % Mass of wet s.
mass_ds_save    =          zeros(round(n_f/i_output),1); % Mass of dry s.
mass_c_save     =          zeros(round(n_f/i_output),1); % Mass of char
mass_a_save     =          zeros(round(n_f/i_output),1); % Mass of ash
%AT

XCells_save =         ones(nx_i,round(n_f/i_output));  % Coord. cell ctrs.
TEMP_save   = temp_i *ones(nx_i,round(n_f/i_output));  % Solid temperature
RHOS_save   = rho_s_i*ones(nx_i,round(n_f/i_output));  % Solid mass density
PSI_save    = psi_i  *ones(nx_i,round(n_f/i_output));  % Porosity
XWS_save    = x_ws_i *ones(nx_i,round(n_f/i_output));  % Vol. frac. wet s.
XDS_save    = x_ds_i *ones(nx_i,round(n_f/i_output));  % Vol. frac. dry s.
XC_save     = x_c_i  *ones(nx_i,round(n_f/i_output));  % Vol. frac. char
XA_save     = x_a_i  *ones(nx_i,round(n_f/i_output));  % Vol. frac. ash
RR1_save    =         zeros(nx_i,round(n_f/i_output)); % R1 reaction rate
RR2_save    =         zeros(nx_i,round(n_f/i_output)); % R2 reaction rate
RR3_save    =         zeros(nx_i,round(n_f/i_output)); % R3 reaction rate
RR4_save    =         zeros(nx_i,round(n_f/i_output)); % R4 reaction rate
Qdot_save   =         zeros(nx_i,round(n_f/i_output)); % Volumetric HRR
YO2G_save   = Y_O2_i *ones(nx_i,round(n_f/i_output));  % Mass frac. of O2
PRESG_save  = pres_i *ones(nx_i,round(n_f/i_output));  % Pressure
MLRPUV_save =         zeros(nx_i,round(n_f/i_output)); % Volumetric MLR
MFLUXG_save =         zeros(nx_i,round(n_f/i_output)); % Gas mass flux
VELG_save   =         zeros(nx_i,round(n_f/i_output)); % Gas velocity

% Set the computational grid
[xRight, xCenter, xLeft1, dV] = mesh(nx_i, dx);

% Arrays needed in exressions of R1-R2-R3-R4 reaction rates
sum_R1_i = rho_ws*(1-psi_ws)*x_ws_i;
sum_R1   = sum_R1_i.*dV;   % Array containing integral term in R1-RR
sum_R1   = max(sum_R1, 1e-14);

sum_R2_i = rho_ds*(1-psi_ds)*x_ds_i;
sum_R2   = sum_R2_i.*dV;   % Array containing integral term in R2-RR
sum_R2   = max(sum_R2, 1e-14);

sum_R3_i = rho_ds*(1-psi_ds)*x_ds_i;
sum_R3   = sum_R3_i.*dV;   % Array containing integral term in R3-RR
sum_R3   = max(sum_R3, 1e-14);

sum_R4_i = rho_c*(1-psi_c)*x_c_i;
sum_R4   = sum_R4_i.*dV;   % Array containing integral term in R4-RR
sum_R4   = max(sum_R4, 1e-14);

% Save output quantities at initial time
n_output            = n_output+1;
time_save(n_output) = time;
fprintf(' n_output = %g, time = %g \n',n_output,time_save(n_output));

% Update array containing coordinates of cell centers
XCells_save(1:nx_save(n_output),n_output) = xCenter(1:nx_save(n_output));

flag_burnout = 0;
% Calculate final value of half-thickness/radius of particle (at burn-out)
delta_f = 1000;
if( ( A_R3 == 0 ) && ( A_R4 == 0 ) )   % (R1)-(R2) reaction model
    delta_f = delta_i*(eta_c_R2*eta_ds_R1*rho_ws*(1-psi_ws) ...
                      /rho_c /(1-psi_c));
end
if( A_R4 ~= 0 )   % (R1)-(R4) reaction model
    % CAUTION: this expression ignores reaction (R3)
    delta_f = delta_i*(eta_a_R4*eta_c_R2*eta_ds_R1*rho_ws*(1-psi_ws) ...
                      /rho_a /(1-psi_a));
end

iter_max = 0;

fid = fopen('LogFile.txt','w');
%--- Begin Time loop ---
while ((n ~= n_f) && (flag_burnout ~= 2) && (time < T_end))
    
    n = n + 1;
    
    % Solution at previous time step (t = time-dt)
    nx_old   = nx_new;
    dV_old   = dV;
    dx_old   = dx;
    dt_old   = dt;    
    temp_old = temp;
    x_ws_old = x_ws;
    x_ds_old = x_ds;
    x_c_old  = x_c;
    x_a_old  = x_a;
    Y_O2_old = Y_O2;
    pres_old = pres;
    
    temp_surf_old = temp_surf;
    
    for i=1:nx_old         
        Y_O2s   = max(Y_O2_old(i),0);
        
        K_R1 = ...
           ( max(0,rho_ws*(1-psi_ws)*x_ws_old(i)*dV_old(i))^n_R1 ) ...
              *sum_R1(i)^(1-n_R1) ...
              *A_R1*exp(-Ta_R1/temp_old(i));
        K_R2 = ...
           ( max(0,rho_ds*(1-psi_ds)*x_ds_old(i)*dV_old(i))^n_R2 ) ...
              *sum_R2(i)^(1-n_R2) ...
              *A_R2*exp(-Ta_R2/temp_old(i));
        K_R3 = ...
           ( max(0,rho_ds*(1-psi_ds)*x_ds_old(i)*dV_old(i))^n_R3 ) ...
                 *sum_R3(i)^(1-n_R3) ...
                 *( (Y_O2s/0.226)^n_O2_R3 ) ...
                 *A_R3*exp(-Ta_R3/temp_old(i));
        %{
        K_R4 = ...
           ( max(0,rho_c*(1-psi_c)*x_c_old(i)*dV_old(i))^n_R4 ) ...
              *sum_R4(i)^(1-n_R4) ...
              *( (Y_O2s/0.226)^n_O2_R4 ) ...
              *A_R4*exp(-Ta_R4/temp_old(i));
              %%AT *A_R4*exp(-Ta_R4/min(temp_old(i),700));
        %}

        sum_R1(i) = sum_R1(i);
        sum_R2(i) = sum_R2(i) + dt * eta_ds_R1*K_R1;
        sum_R3(i) = sum_R3(i) + dt * eta_ds_R1*K_R1;
        sum_R4(i) = sum_R4(i) + dt * (eta_c_R2*K_R2+eta_c_R3*K_R3);
    end
    
    % Solution at new time step (t = time)
    
    time = time + dt;
    if(mod(n,i_output)==0)
        n_output            = n_output+1;
        time_save(n_output) = time;
        fprintf(' \n');
        fprintf(' n_output = %g, time = %g \n', ...
                  n_output,time_save(n_output));
    end

    %%AT End of flaming residence time
    t_flaming = 120;
    if(time >= t_flaming)
        G = sigma*T_g^4;
    end
    %%AT
    
    flag_iter = 0;
    
    %--- Begin iterative loop ---    
    % Initial guess
    temp_newiter = temp_old;
    x_ws_newiter = x_ws_old;
    x_ds_newiter = x_ds_old;
    x_c_newiter  = x_c_old;
    x_a_newiter  = x_a_old;
    Y_O2_newiter = Y_O2_old;
    pres_newiter = pres_old;
    
    temp_surf_newiter = temp_surf_old;
    %%AT
    temp_olditer = temp_newiter;
    x_ws_olditer = x_ws_newiter;
    x_ds_olditer = x_ds_newiter;
    x_c_olditer  = x_c_newiter;
    x_a_olditer  = x_a_newiter;
    Y_O2_olditer = Y_O2_newiter;
    pres_olditer = pres_newiter;
    
    temp_surf_olditer = temp_surf_newiter;
    %%AT
        
    for iter=1:1000
        
    % Solution at new time step, at previous iteration
    %%AT
    %{
    temp_olditer = temp_newiter;
    x_ws_olditer = x_ws_newiter;
    x_ds_olditer = x_ds_newiter;
    x_c_olditer  = x_c_newiter;
    x_a_olditer  = x_a_newiter;
    Y_O2_olditer = Y_O2_newiter;
    pres_olditer = pres_newiter;
    
    temp_surf_olditer = temp_surf_newiter;
    %}
    %%AT
    %%{
    temp_olditer = 0.5*temp_newiter + 0.5*temp_olditer;
    x_ws_olditer = 0.5*x_ws_newiter + 0.5*x_ws_olditer;
    x_ds_olditer = 0.5*x_ds_newiter + 0.5*x_ds_olditer;
    x_c_olditer  = 0.5*x_c_newiter  + 0.5*x_c_olditer;
    x_a_olditer  = 0.5*x_a_newiter  + 0.5*x_a_olditer;
    Y_O2_olditer = 0.5*Y_O2_newiter + 0.5*Y_O2_olditer;
    pres_olditer = 0.5*pres_newiter + 0.5*pres_olditer;
    
    temp_surf_olditer = 0.5*temp_surf_newiter + 0.5*temp_surf_olditer;
    %}
    %%AT
        
    % Calculate solution at new time step, at new iteration
    
    % - Calculate temperature
    %
    % Discretization (semi-implicit):
    %
    %   (Tp(n+1)-Tp(n))/dt = RHS1 + RHS2 where RHS1 is energy change
    %                            due to reactions (R1)-(R4) and RHS2 is
    %                            energy change due to convection/diffusion
    %
    %   Evaluation of RHS1:
    %      RHS1 = RHS1(n+1)
    %   Evaluation of RHS2:
    %      RHS2 = 0.5*RHS2(n+1) + 0.5*RHS2(n)
    %    
    
    [a, b, c, d, temp_surf_newiter] = ...
             energy_conservation(dt, temp_old, temp_olditer, ...
             x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
             x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
             Y_O2_olditer, temp_surf_olditer, h_conv, ...
             sum_R1, sum_R2, sum_R3, sum_R4, ...
             lambda, xRight, xCenter, dV_old, nx_old);         
    temp_newiter = tri(a,b,c,d);
    temp_newiter = temp_newiter';
        
    % Check for convergence of iterative loop
    DeltaTemp_max_iter = max( abs(temp_newiter-temp_olditer) );
    DeltaTemp_max_dt   = max( abs(temp_newiter-temp_old    ) );
    if (mod(n,i_output)==0)
        fprintf(' DeltaTemp_max_iter,iter = %g %g \n', ...
                  DeltaTemp_max_iter,iter);
    end
    %AT if (iter > 1)
    if (iter >= 50)
        fprintf(fid,' Time,DeltaTemp_max_iter,iter = %g %g %g \n', ...
                      time,DeltaTemp_max_iter,iter);
        flag_iter = 1;
    end
    
    % - Calculate solid species volume fractions and cell volumes
    [x_ws_newiter, x_ds_newiter, x_c_newiter, x_a_newiter, dV] = ...       
    mass_conservation_v2(dt, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, sum_R1, sum_R2, sum_R3, sum_R4, ...
    lambda_xk, dV_old, nx_old);
    x_ws_newiter = x_ws_newiter';
    x_ds_newiter = x_ds_newiter';
    x_c_newiter  = x_c_newiter';
    x_a_newiter  = x_a_newiter';
    
    % Check for convergence of iterative loop
    Deltax_ws_max_iter = max( abs(x_ws_newiter-x_ws_olditer) );
    Deltax_ds_max_iter = max( abs(x_ds_newiter-x_ds_olditer) );
    Deltax_c_max_iter  = max( abs(x_c_newiter - x_c_olditer) );
    Deltax_a_max_iter  = max( abs(x_a_newiter - x_a_olditer) );
    Deltax_k_max_iter  = max( [ Deltax_ws_max_iter,Deltax_ds_max_iter, ...
                                Deltax_c_max_iter ,Deltax_a_max_iter ] );
    Deltax_ws_max_dt   = max( abs(x_ws_newiter-x_ws_old) );
    Deltax_ds_max_dt   = max( abs(x_ds_newiter-x_ds_old) );
    Deltax_c_max_dt    = max( abs(x_c_newiter - x_c_old) );
    Deltax_a_max_dt    = max( abs(x_a_newiter - x_a_old) );
    Deltax_k_max_dt    = max( [ Deltax_ws_max_dt,Deltax_ds_max_dt, ...
                                Deltax_c_max_dt ,Deltax_a_max_dt ] );
    if (mod(n,i_output)==0)
        fprintf(' Deltax_k_max_iter,iter  = %g %g \n', ...
                  Deltax_k_max_iter,iter);
    end
    %AT if (iter > 1)
    if (iter >= 50)
        fprintf(fid,' Time,Deltax_k_max_iter,iter  = %g %g %g \n', ...
                      time,Deltax_k_max_iter,iter);
        flag_iter = 1;
    end
    
    if( ( A_R3 == 0 ) && ( A_R4 == 0 ) )   % (R1)-(R2) reaction model
    %%AT if( 1 == 0 ) % Uncomment this line to force calc. of Y_O2 & pres
        Y_O2_newiter       = Y_O2_old;
        pres_newiter       = pres_old;
        DeltaYO2_max_iter  = 0;
        DeltaYO2_max_dt    = 0;
        DeltaPres_max_iter = 0;
        DeltaPres_max_dt   = 0;
    else
        % *** Begin calculation of Y_O2 and pres ***
        
    % - Calculate oxygen mass fraction
    %
    % Discretization (semi-implicit):
    %
    %   (Y_O2(n+1)-Y_O2(n))/dt = RHS1 + RHS2 where RHS1 is O2 mass change
    %                            due to reactions (R1)-(R4) and RHS2 is
    %                            O2 mass change due to convection/diffusion
    %
    %   Evaluation of RHS1:
    %      RHS1 = RHS1(n+1)
    %   Evaluation of RHS2:
    %      RHS2 = 0.5*RHS2(n+1) + 0.5*RHS2(n)
    %
    [a, b, c, d] = oxygen_mass_conservation(dt, Y_O2_old, ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, lambda, h_conv, ...
    sum_R1, sum_R2, sum_R3, sum_R4, ...
    xRight, xCenter, dV_old, nx_old);
    Y_O2_newiter = tri(a,b,c,d);
    Y_O2_newiter = Y_O2_newiter';

    %AT
    Y_O2_newiter = max(0,min(1,Y_O2_newiter)); % Enforce 0 <= Y_O2 <= 1
    %AT
        
    % Check for convergence of iterative loop
    DeltaYO2_max_iter = max( abs(Y_O2_newiter-Y_O2_olditer) );
    DeltaYO2_max_dt   = max( abs(Y_O2_newiter-Y_O2_old    ) );
    if(mod(n,i_output)==0)
        fprintf(' DeltaYO2_max_iter,iter  = %g %g \n', ...
                  DeltaYO2_max_iter,iter);
    end
    %AT if (iter > 1)
    if (iter >= 50)
        fprintf(fid,' Time,DeltaYO2_max_iter,iter  = %g %g %g \n',...
                      time,DeltaYO2_max_iter,iter);
    end
        
    % - Calculate pressure
    %
    % Discretization (semi-implicit):
    %
    %   (pres(n+1)-pres(n))/dt = RHS1 + RHS2 where RHS1 is pressure
    %                              change due to reactions (R1)-(R4) and
    %                              RHS2 is pressure change due to diffusion
    %   Evaluation of RHS1:
    %      RHS1 = RHS1(n)
    %   Evaluation of RHS2:
    %      RHS2 = 0.5*RHS2(n+1) + 0.5*RHS2(n)
    %
    [a, b, c, d] = pressure_equation_QS( ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, pres_olditer, lambda_pres, ...
    sum_R1, sum_R2, sum_R3, sum_R4, ...
    xRight, xCenter, dV_old, nx_old);
    pres_newiter = tri(a,b,c,d);
    pres_newiter = pres_newiter';
        
    % Check for convergence of iterative loop
    DeltaPres_max_iter = max( abs(pres_newiter-pres_olditer) );
    DeltaPres_max_dt   = max( abs(pres_newiter-pres_old    ) );
    if(mod(n,i_output)==0)
        fprintf(' DeltaPres_max_iter,iter = %g %g \n', ...
                  DeltaPres_max_iter,iter);
    end
    %AT if (iter > 1)
    if (iter >= 50)
        fprintf(fid,' Time,DeltaPres_max_iter,iter = %g %g %g \n',...
                      time,DeltaPres_max_iter,iter);
    end
    
        % *** End calculation of Y_O2 and pres ***
    end
    
    % Convergence criteria
    if( ( DeltaTemp_max_iter <= (0.01*DeltaTemp_max_dt/lambda) ) & ...
        ( Deltax_k_max_iter  <= (0.01*Threshold_xk) )            & ...
        ( DeltaYO2_max_iter  <= (0.01*DeltaYO2_max_dt /lambda) ) & ...
        ( DeltaPres_max_iter <= (0.01*Threshold_pres) ) )
        break;
    end

    if(iter==1000)
        fprintf(' Problem in iterative loop: break at time %g \n',time);
        fprintf(' \n');
        return;
    end
    
    %--- End iterative loop ---
    end
    
    % Solution at new time step
    temp = temp_newiter;
    x_ws = x_ws_newiter;
    x_ds = x_ds_newiter;
    x_c  = x_c_newiter;
    x_a  = x_a_newiter;
    Y_O2 = Y_O2_newiter;
    pres = pres_newiter;
    
    % Time step restriction:
    % - Update dt such that abs( Tp(n+1)-Tp(n) ) <= Threshold_Temp   
    % Safety: (1) Limit the time step dt so that variations in temperature
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    %AT
    coefmin = 0.90;
    coefmax = 1.10;
    %AT
    if( DeltaTemp_max_dt > 0 )
        dt_Temp = dt_old*Threshold_Temp/DeltaTemp_max_dt;  % Constraint (1)
        dt_Temp = max(coefmin*dt_old, ...
                             min(coefmax*dt_old,dt_Temp)); % Constraint (2)
    else
        dt_Temp = coefmax*dt_old;
    end

    % Time step restriction:
    % - Update dt such that abs( xk(n+1)-xk(n) ) <= Threshold_xk
    % Safety: (1) Limit the time step dt so that variations in x_k
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( Deltax_k_max_dt > 0 )
        dt_xk = dt_old*Threshold_xk/Deltax_k_max_dt;     % Constraint (1)
        dt_xk = max(coefmin*dt_old,min(coefmax*dt_old,dt_xk));   % Constraint (2)
    else
        dt_xk = coefmax*dt_old;
    end
     
    % Time step restriction:
    % - Update dt such that abs( Y_O2(n+1)-Y_O2(n) ) <= Threshold_YO2         
    % Safety: (1) Limit the time step dt so that variations in Y_O2
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( DeltaYO2_max_dt > 0 )
        dt_YO2 = dt_old*Threshold_YO2/DeltaYO2_max_dt;   % Constraint (1)
        dt_YO2 = max(coefmin*dt_old, ...
                            min(coefmax*dt_old,dt_YO2)); % Constraint (2)
    else
        dt_YO2 = coefmax*dt_old;
    end

    % Time step restriction:
    % - Update dt such that abs( pres(n+1)-pres(n) ) <= Threshold_pres
    % Safety: (1) Limit the time step dt so that variations in pres
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( DeltaPres_max_dt > 0 )
        dt_pres = dt_old*Threshold_pres/DeltaPres_max_dt;  % Constraint (1)
        dt_pres = max(coefmin*dt_old, ...
                             min(coefmax*dt_old,dt_pres)); % Constraint (2)
    else
        dt_pres = coefmax*dt_old;
    end
    
    dt = min( [ dt_Temp,dt_xk,dt_YO2,dt_pres ] );
    dt = min(dt,dt_max); % Safety: do not let dt go to very large values
    %%AT dt = dt_old; % Uncomment this line for tests with fixed dt

    iter_max = max(iter,iter_max);
 
    if(mod(n,i_output)==0)
        fprintf(' time, dt             = %g %g \n',time,dt);
        fprintf(' iter                 = %g \n',iter);
        fprintf(' max(|temp-temp_old|) = %g \n',DeltaTemp_max_dt);
        fprintf(' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max_dt );
        fprintf(' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max_dt );
        fprintf(' max(|pres-pres_old|) = %g \n',DeltaPres_max_dt);
    end        
    if( ( flag_iter == 1 ) & ( iter >= 50 ) )
        fprintf(' *** WARNING *** \n');
        fprintf(' flag_iter = %g \n', flag_iter);
        fprintf(' iter                 = %g \n',iter);
        fprintf(' time, dt             = %g %g \n',time,dt);
        fprintf(' max(|temp-temp_old|) = %g \n',DeltaTemp_max_dt);
        fprintf(' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max_dt );
        fprintf(' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max_dt );
        fprintf(' max(|pres-pres_old|) = %g \n',DeltaPres_max_dt);
    end
    %AT if( ( flag_iter == 1 ) & ( iter > 1 ) )
    if( ( flag_iter == 1 ) & ( iter >= 50 ) )
        fprintf(fid,' *** WARNING *** \n');
        fprintf(fid,' flag_iter = %g \n', flag_iter);
        fprintf(fid,' iter                 = %g \n',iter);
        fprintf(fid,' time, dt             = %g %g \n',time,dt);
        fprintf(fid,' max(|temp-temp_old|) = %g \n',DeltaTemp_max_dt);
        fprintf(fid,' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max_dt );
        fprintf(fid,' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max_dt );
        fprintf(fid,' max(|pres-pres_old|) = %g \n',DeltaPres_max_dt);
    end
 
    % Update computational grid (in case of volume change)
    [xRight, xCenter, dx] = moveMesh(nx_old, dV);
    
    % Check
    if(mod(n,i_output)==0)
        sample_thickness = 0;
        for i=1:nx_old
            sample_thickness = sample_thickness + dx_old(i);
        end
        fprintf(' delta,     sample_thickness (radius) = %g %g \n', ...
                                          xRight(nx_old),sample_thickness);
    end
    
    % Begin remeshing    
    % Remeshing, step 1: update nx, dx
    dx_old = dx;
    dV_old = dV;
    [nx_new, dx_new] = remeshing(nx_old, xRight(nx_old));
    if(mod(n,i_output)==0)
        fprintf(' dx = %g, nx = %g \n',dx_new(1),nx_new);
    end
    
    % Remeshing, step 2: update computational grid
    xCenter_old = xCenter;
    [xRight, xCenter, xLeft1, dV] = mesh(nx_new, dx_new);
    dx     = dx_new;
    volume = sum(dV);
    
    % Check
    if(mod(n,i_output)==0)
        sample_thickness = 0;
        for i=1:nx_new
            sample_thickness = sample_thickness + dx_new(i);
        end
        fprintf(' delta_new, sample_thickness (radius) = %g %g \n', ...
                                          xRight(nx_new),sample_thickness);
    end
    
    % Remeshing, step 3: interpolate solution on new computational grid
    temp_old   = temp;
    x_ws_old   = x_ws;
    x_ds_old   = x_ds;
    x_c_old    = x_c;
    x_a_old    = x_a;
    Y_O2_old   = Y_O2;
    pres_old   = pres;
    sum_R1_old = sum_R1;
    sum_R2_old = sum_R2;
    sum_R3_old = sum_R3;
    sum_R4_old = sum_R4;
   
   [temp, x_ws, x_ds, x_c, x_a, Y_O2, pres, ...
          sum_R1, sum_R2, sum_R3, sum_R4]=interpolateOnNewMesh( ...
              nx_new, xCenter, volume, dV, nx_old, xCenter_old, dV_old, ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, Y_O2_old, pres_old, ...
                        sum_R1_old, sum_R2_old, sum_R3_old, sum_R4_old);
   
    temp = temp';
    x_ws = x_ws';
    x_ds = x_ds';
    x_c  = x_c';
    x_a  = x_a';
    Y_O2 = Y_O2';
    pres = pres';
    % End remeshing
    
    for i=1:nx_new
        psi_sg(i) = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
                  + psi_c *x_c(i)  + psi_a *x_a(i);
        if( (1-psi_sg(i)) > 0 )
            rho_s(i) = ( rho_ws*(1-psi_ws)*x_ws(i) ...
                       + rho_ds*(1-psi_ds)*x_ds(i) ...
                       + rho_c *(1-psi_c) *x_c(i)  ...
                       + rho_a *(1-psi_a) *x_a(i) )/(1-psi_sg(i));
        else
            rho_s(i) = rho_a;
        end
    end
     
    % Calculate the net surface heat flux and the surface temperature
    
    % Porosity
    psi       = psi_sg(nx_new);
    % Effective thermal conductivity
    %  - Porous medium treatment: 
    %    keff = k_s + psi x k_g
    k_ws      = k0_ws*(temp(nx_new)/300)^nk_ws;
    k_ds      = k0_ds*(temp(nx_new)/300)^nk_ds;
    k_c       = k0_c *(temp(nx_new)/300)^nk_c;
    k_a       = k0_a *(temp(nx_new)/300)^nk_a;    
    k_s       = k_ws*x_ws(nx_new) + k_ds*x_ds(nx_new) ...
              + k_c*x_c(nx_new) + k_a*x_a(nx_new); % Conduct. in sol. phase
    gamma     = gamma_ws*x_ws(nx_new) + gamma_ds*x_ds(nx_new) ...
              + gamma_c *x_c(nx_new)  + gamma_a *x_a(nx_new);
    k_rad     = gamma*sigma*temp(nx_new)^3; 
    k_s       = k_s + k_rad;
    k_g       = k_g0*(temp(nx_new)/300)^(0.76);    % Conduct. in gas phase
    keff_surf = k_s + psi*k_g;
    eps_surf  = eps_ws*x_ws(nx_new) + eps_ds*x_ds(nx_new) ...
              + eps_c *x_c(nx_new)  + eps_a *x_a(nx_new);
    dx_surf   = xRight(nx_new)-xCenter(nx_new);
    tempc_nx  = temp(nx_new);
    
    % Calculation of convective heat transfer coefficient
    T_film       = 0.5*(temp_surf_old+T_g);      % Film temperature [K]
    nu_g         = nu_g0*(T_film/300)^(1.76); % Kin. visc. at film temp.
    k_g          = k_g0*(T_film/300)^(0.76);  % Conductivity at film temp.
    if geometry=="rectangle"
        D_eff    = 2*xRight(nx_new); % Effective diameter of particle [m]
        Re_D     = u_g*D_eff/nu_g;
        Nu_D     = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
    elseif geometry=="cylinder"
        D_eff    = 2*xRight(nx_new);
        Re_D     = u_g*D_eff/nu_g;
        Nu_D     = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
    elseif geometry=="sphere"
        D_eff    = 2*xRight(nx_new);
        k_g      = k_g0*(T_g/300)^(0.76);        % Conductivity at gas temp.
        nu_g     = nu_g0*(T_g/300)^(1.76);       % Kin. visc. at gas temp.
        nu_s     = nu_g0*(temp_surf_old/300)^(1.76); % Visc. at surf. temp.
        rho_g    = pres_g*MW_g/R/T_g;
        rho_surf = pres_g*MW_g/R/temp_surf_old;
        mu_g     = rho_g*nu_g;
        mu_s     = rho_surf*nu_s;
        Re_D     = u_g*D_eff/nu_g; % Calculated at gas temperature
        Nu_D     = 2+(0.4*Re_D^0.5+0.06*Re_D^(2/3)) ...
                     *Pr^0.4*(mu_g/mu_s)^0.25;                                                     
    else
        fprintf(...
    ' Invalid geometry, valid geometries are: cylinder/sphere/rectangle ');
        return;
    end
    h_conv = Nu_D*k_g/D_eff;   % Convective heat transfer coef. [W/m2/K]
    %%AT
    %%{
    %AT h_conv0 = 10; % Model of Lautenberger & Fernandez-Pello (2009)
    h_conv0 = h_conv;
    %AT
    
    mdotPUA_surf = 0;
    % Permeability
    Kperm = Kperm_ws*x_ws(nx_new) + Kperm_ds*x_ds(nx_new) ...
          + Kperm_c *x_c(nx_new)  + Kperm_a *x_a(nx_new);
    % Gas kinematic viscosity
    nu_g  = nu_g0*(temp(nx_new)/300)^(1.76);

    if( pres(nx_new-1) > pres(nx_new) )
        % Case for which mdotPUA > 0
        mdotPUA_surf = (Kperm/nu_g) ...
        *(pres(nx_new-1)-pres(nx_new))/(xCenter(nx_new)-xCenter(nx_new-1));
    end
    %%AT if(mdotPUA_surf > 0)
    if(mdotPUA_surf > 1e-14)
        h_conv = mdotPUA_surf*cp_g0/( exp(mdotPUA_surf*cp_g0/h_conv0)-1 );
    else
        h_conv = h_conv0;
    end
    %}
    %%AT
    
    % Calculate the surface heat flux
    h_rad     = eps_surf*sigma*(temp_surf_newiter^3);
    Bi        = (h_conv+h_rad)*dx_surf/keff_surf;

    temp_surf = ( tempc_nx ...
                    + (h_conv*T_g+eps_surf*G)*dx_surf/keff_surf )/(1+Bi);
    q_surf    = eps_surf*G - eps_surf*sigma*temp_surf^4 ...
              + h_conv*(T_g-temp_surf);
    
    % Calculate the surface oxygen mass flux
    psi          = psi_sg(nx_new);
    rho_g        = pres_g*MW_g/R/temp(nx_new);
    D_g          = nu_g0*(temp(nx_new)/300)^(1.76);
    rhopsiD_surf = psi*rho_g*D_g;
    dx_surf      = xRight(nx_new)-xCenter(nx_new);
    h_mass       = h_conv/cp_g0;
    Bi           = h_mass*dx_surf/rhopsiD_surf;
    Y_O2c_nx     = Y_O2(nx_new);
    
    Y_O2_surf    = (Bi*Y_g_O2+Y_O2c_nx)/(1+Bi);
    mdot_surf    = h_mass*(Y_g_O2-Y_O2_surf);
    
    % Save data
    if(mod(n,i_output)==0)            
        % Calculate the mass loss rate 
        %  - Mass loss rate per unit volume [kg/s/m3]
        for i=1:nx_new
            Y_O2s = max(Y_O2(i),0);        
        
            K_R1 = ( max(0,rho_ws*(1-psi_ws)*x_ws(i)*dV(i))^n_R1 ) ...
                   *sum_R1(i)^(1-n_R1) ...
                   *A_R1*exp(-Ta_R1/temp(i));
            K_R2 = ( max(0,rho_ds*(1-psi_ds)*x_ds(i)*dV(i))^n_R2 ) ...
                   *sum_R2(i)^(1-n_R2) ...
                   *A_R2*exp(-Ta_R2/temp(i));
            K_R3 = ( max(0,rho_ds*(1-psi_ds)*x_ds(i)*dV(i))^n_R3 ) ...
                   *sum_R3(i)^(1-n_R3) ...
                   *( (Y_O2s/0.226)^n_O2_R3 ) ...
                   *A_R3*exp(-Ta_R3/temp(i));
            K_R4 = ( max(0,rho_c*(1-psi_c)*x_c(i)*dV(i))^n_R4 ) ...
                   *sum_R4(i)^(1-n_R4) ...
                   *( (Y_O2s/0.226)^n_O2_R4 ) ...
                   *A_R4*exp(-Ta_R4/min(temp(i),Tmax_R4));
                   %%AT *A_R4*exp(-Ta_R4/temp(i));
      
            MLRPUV(i) = (1-eta_ds_R1)*K_R1/dV(i) ...
                      + (1-eta_c_R2) *K_R2/dV(i) ...
                      + (1-eta_c_R3) *K_R3/dV(i) ...
                      + (1-eta_a_R4) *K_R4/dV(i);

            %AT
            MLR_R1PUV(i) = (1-eta_ds_R1)*K_R1/dV(i);
            MLR_R2PUV(i) = (1-eta_c_R2) *K_R2/dV(i);
            MLR_R3PUV(i) = (1-eta_c_R3) *K_R3/dV(i);
            MLR_R4PUV(i) = (1-eta_a_R4) *K_R4/dV(i);
            %AT

            HRRPUV(i) = (1-eta_ds_R1)*K_R1*DeltaH_R1/dV(i) ...
                      + (1-eta_c_R2) *K_R2*DeltaH_R2/dV(i) ...
                      + (1-eta_c_R3) *K_R3*DeltaH_R3/dV(i) ...
                      + (1-eta_a_R4) *K_R4*DeltaH_R4/dV(i);
                  
            RR1_save(i,n_output) = K_R1/dV(i);
            RR2_save(i,n_output) = K_R2/dV(i);
            RR3_save(i,n_output) = K_R3/dV(i);
            RR4_save(i,n_output) = K_R4/dV(i);
        end
        
        mdotPUA(1) = 0;
        % Permeability
        Kperm = Kperm_ws*x_ws(1) + Kperm_ds*x_ds(1) ...
              + Kperm_c *x_c(1)  + Kperm_a *x_a(1);
        % Gas kinematic viscosity
        nu_g  = nu_g0*(temp(1)/300)^(1.76);
        if (pres(1) < pres(2))
            % Case for which mdotPUA < 0
            mdotPUA(1) = (Kperm/nu_g) ...
                         *(pres(1)-pres(2))/(xCenter(2)-xCenter(1));
        end
        for i=2:(nx_new-1)
            % Permeability
            Kperm = Kperm_ws*x_ws(i) + Kperm_ds*x_ds(i) ...
                  + Kperm_c *x_c(i)  + Kperm_a* x_a(i);   
            % Gas kinematic viscosity
            nu_g  = nu_g0*(temp(i)/300)^(1.76);          
            if( (pres(i-1) > pres(i)) & (pres(i) > pres(i+1)) )
                % Case fhor which mdotPUA > 0
                mdotPUA(i) = (Kperm/nu_g) ...
                            *(pres(i-1)-pres(i))/(xCenter(i)-xCenter(i-1));
            elseif ( (pres(i-1) < pres(i)) & (pres(i) < pres(i+1)) )
                % Case for which mdotPUA < 0
                mdotPUA(i) = (Kperm/nu_g) ...
                            *(pres(i)-pres(i+1))/(xCenter(i+1)-xCenter(i));
            else
                % Case for which dp/dx ~ 0
                mdotPUA(i) = 0;
            end
        end
        mdotPUA(nx_new) = 0;
        % Permeability
        Kperm = Kperm_ws*x_ws(nx_new) + Kperm_ds*x_ds(nx_new) ...
              + Kperm_c *x_c(nx_new)  + Kperm_a *x_a(nx_new);
        % Gas kinematic viscosity
        nu_g  = nu_g0*(temp(nx_new)/300)^(1.76);
        if( pres(nx_new-1) > pres(nx_new) )
            % Case for which mdotPUA > 0
            mdotPUA(nx_new) = (Kperm/nu_g) ...
        *(pres(nx_new-1)-pres(nx_new))/(xCenter(nx_new)-xCenter(nx_new-1));
        end
        
        for i=1:nx_new
            rho_g  = (pres_g+pres(i))*MW_g/R/temp(i);   % Gas mass density
            vel(i) = mdotPUA(i)/rho_g;
        end
        
        %  - Total mass loss rate [kg/s]
        MLRtot_new = 0;
        for i=1:nx_new
            MLRtot_new = MLRtot_new + MLRPUV(i)*dV(i);
        end
        if geometry=="rectangle"
            MLRtot_new = MLRtot_new*2; % We solve for half of the particle
        end

        %AT
        MLR_R1tot_new = 0;
        MLR_R2tot_new = 0;
        MLR_R3tot_new = 0;
        MLR_R4tot_new = 0;
        for i=1:nx_new
            MLR_R1tot_new = MLR_R1tot_new + MLR_R1PUV(i)*dV(i);
            MLR_R2tot_new = MLR_R2tot_new + MLR_R2PUV(i)*dV(i);
            MLR_R3tot_new = MLR_R3tot_new + MLR_R3PUV(i)*dV(i);
            MLR_R4tot_new = MLR_R4tot_new + MLR_R4PUV(i)*dV(i);
        end
        if geometry=="rectangle"
            MLR_R1tot_new = MLR_R1tot_new*2;
            MLR_R2tot_new = MLR_R2tot_new*2;
            MLR_R3tot_new = MLR_R3tot_new*2;
            MLR_R4tot_new = MLR_R4tot_new*2;
        end
        %AT
        
        DeltaQ_max_save(1,n_output) = DeltaTemp_max_dt;
        DeltaQ_max_save(2,n_output) = Deltax_k_max_dt;
        DeltaQ_max_save(3,n_output) = DeltaYO2_max_dt;
        DeltaQ_max_save(4,n_output) = DeltaPres_max_dt;
         
        MLR_save(n_output)       = MLRtot_new;
        q_surf_save(n_output)    = q_surf;
        temp_surf_save(n_output) = temp_surf;
        temp_back_save(n_output) = temp(1);
        delta_save(n_output)     = xRight(nx_new);
        dt_save(n_output)        = dt;     
        nx_save(n_output)        = nx_new;       
        YO2_surf_save(n_output)  = Y_O2_surf;
        YO2_nx_save(n_output)    = Y_O2(nx_new);
        
        XCells_save(1:nx_new,n_output) = xCenter(1:nx_new);
        TEMP_save(1:nx_new,n_output)   = temp(1:nx_new);
        RHOS_save(1:nx_new,n_output)   = rho_s(1:nx_new);
        PSI_save(1:nx_new,n_output)    = psi_sg(1:nx_new);
        XWS_save(1:nx_new,n_output)    = x_ws(1:nx_new);
        XDS_save(1:nx_new,n_output)    = x_ds(1:nx_new);
        XC_save(1:nx_new,n_output)     = x_c(1:nx_new);
        XA_save(1:nx_new,n_output)     = x_a(1:nx_new);
        YO2G_save(1:nx_new,n_output)   = Y_O2(1:nx_new);
        PRESG_save(1:nx_new,n_output)  = pres(1:nx_new);
        MLRPUV_save(1:nx_new,n_output) = MLRPUV(1:nx_new);
        MFLUXG_save(1:nx_new,n_output) = mdotPUA(1:nx_new);
        VELG_save(1:nx_new,n_output)   = vel(1:nx_new);
        Qdot_save(1:nx_new,n_output)   = HRRPUV(1:nx_new);
        
        [RR1_peak_save(n_output),i1] = max(RR1_save(:,n_output));
        xRR1_peak_save(n_output)     = xCenter(i1);
        [RR2_peak_save(n_output),i2] = max(RR2_save(:,n_output));
        xRR2_peak_save(n_output)     = xCenter(i2);
        [RR3_peak_save(n_output),i3] = max(RR3_save(:,n_output));
        xRR3_peak_save(n_output)     = xCenter(i3);
        [RR4_peak_save(n_output),i4] = max(RR4_save(:,n_output));
        xRR4_peak_save(n_output)     = xCenter(i4);
        
        h_conv_save(n_output)        = h_conv;

        %AT
        MLR_R1_save(n_output)    = MLR_R1tot_new;
        MLR_R2_save(n_output)    = MLR_R2tot_new;
        MLR_R3_save(n_output)    = MLR_R3tot_new;
        MLR_R4_save(n_output)    = MLR_R4tot_new;
        temp_max_save(n_output)  = max(TEMP_save(:,n_output));
        YO2_back_save(n_output)  = YO2G_save(1,n_output);
        YO2_min_save(n_output)   = min(YO2G_save(:,n_output));
        pres_surf_save(n_output) = PRESG_save(nx_new,n_output);
        pres_back_save(n_output) = PRESG_save(1,n_output);
        pres_max_save(n_output)  = max(PRESG_save(:,n_output));
        MLR_surf_save(n_output)  = MFLUXG_save(nx_new,n_output);
        vel_surf_save(n_output)  = VELG_save(nx_new,n_output);
        vel_back_save(n_output)  = VELG_save(1,n_output);
        vel_max_save(n_output)   = max(VELG_save(:,n_output));
        temp_xRR1_peak_save(n_output) = TEMP_save(i1,n_output);
        temp_xRR2_peak_save(n_output) = TEMP_save(i2,n_output);
        temp_xRR3_peak_save(n_output) = TEMP_save(i3,n_output);
        temp_xRR4_peak_save(n_output) = TEMP_save(i4,n_output);
        RR1_surf_save(n_output)  = RR1_save(nx_new,n_output);
        RR1_back_save(n_output)  = RR1_save(1,n_output);
        RR2_surf_save(n_output)  = RR2_save(nx_new,n_output);
        RR2_back_save(n_output)  = RR2_save(1,n_output);
        RR3_surf_save(n_output)  = RR3_save(nx_new,n_output);
        RR3_back_save(n_output)  = RR3_save(1,n_output);
        RR4_surf_save(n_output)  = RR4_save(nx_new,n_output);
        RR4_back_save(n_output)  = RR4_save(1,n_output);
        mass_ws_save(n_output) = 0;
        for i=1:nx_new
	        mass_ws_save(n_output) = mass_ws_save(n_output) ...
                                   + rho_ws*(1-psi_ws)*x_ws(i)*dV(i);
        end
        mass_ds_save(n_output) = 0;
        for i=1:nx_new
	        mass_ds_save(n_output) = mass_ds_save(n_output) ...
                                   + rho_ds*(1-psi_ds)*x_ds(i)*dV(i);
        end
        mass_c_save(n_output) = 0;
        for i=1:nx_new
	        mass_c_save(n_output) = mass_c_save(n_output) ...
                                   + rho_c*(1-psi_c)*x_c(i)*dV(i);
        end
        mass_a_save(n_output) = 0;
        for i=1:nx_new
	        mass_a_save(n_output) = mass_a_save(n_output) ...
                                   + rho_a*(1-psi_a)*x_a(i)*dV(i);
        end
        %AT
    end
    
    % Burnout criterion
    %   (R1)-(R2) reaction model
    if( ( flag_burnout == 0 ) && ( A_R3 == 0 ) && ( A_R4 == 0 ) )
        if( ( eta_c_R2 == 0 ) && ( (xRight(nx_new)/delta_i) < 0.01 ) )
            % Particle has no solid residue and is shrinking to zero size
            flag_burnout  = 1;
            n_last_output = n_output;   % Stop simulation
            fprintf(' \n');
            fprintf(' Burnout time = %g \n',time);
            fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
            fprintf(' Thickness at burnout time = %g \n',delta_f);
        end
        if( ( eta_c_R2 ~= 0 ) && ( x_c(1) > 0.999 ) )
            % Particle has char residue and is not shrinking to zero size
            flag_burnout  = 1;
            n_last_output = n_output+1; % Sim. until time for last output
            % n_last_output = n_output+10;
            fprintf(' \n');
            fprintf(' Burnout time = %g \n',time);
            fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
            fprintf(' Thickness at burnout time = %g \n',delta_f);
        end
    end

    %   (R1)-(R4) reaction model
    if( ( flag_burnout == 0 ) && ( A_R4 ~= 0 ) )
        if( ( eta_a_R4 == 0 ) && ( (xRight(nx_new)/delta_i) < 0.01 ) )
            % Particle has no solid residue and is shrinking to zero size
            flag_burnout  = 1;
            n_last_output = n_output;   % Stop simulation
            fprintf(' \n');
            fprintf(' Burnout time = %g \n',time);
            fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
            fprintf(' Thickness at burnout time = %g \n',delta_f);
        end
        %AT if( ( eta_a_R4 ~= 0 ) && ( x_a(1) > 0.99 ) )
        if( ( eta_a_R4 ~= 0 ) && ( x_a(1) > 0.99 ) )
            % When close to burnout,  allow for larger time steps in order
            % to simulate the particle cooling with larger values of dt
            Threshold_Temp = 1;
            dt_max         = 10*dt_max_i;
        end
        %AT
        if( ( eta_a_R4 ~= 0 ) && ( x_a(1) > 0.999 ) )
            % Particle has ash residue and is not shrinking to zero size
            flag_burnout  = 1;
            n_last_output = n_output+1; % Sim. until time for last output
            % n_last_output = n_output+10;
            fprintf(' \n');
            fprintf(' Burnout time = %g \n',time);
            fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
            fprintf(' Thickness at burnout time = %g \n',delta_f);
        end
    end
    
    if( ( flag_burnout == 1 ) & ( n_output == n_last_output ) )
        flag_burnout = 2;
    end

end

%--- Enf Time loop ---
fclose(fid);

fprintf(' \n');
fprintf(' iter_max = %g \n',iter_max);

%% Output data
%{
csvwrite("time.csv",            time_save(1:n_output));
csvwrite("MLR.csv",             MLR_save(1:n_output));
csvwrite("number_cells.csv",    nx_save(1:n_output));
csvwrite("time.csv",            time_save(1:n_output));
csvwrite("MLR.csv",             MLR_save(1:n_output));
csvwrite("delta.csv",           delta_save(1:n_output));
csvwrite("Temp_surf.csv",       temp_surf_save(1:n_output));
csvwrite("Temp_core.csv",       temp_back_save(1:n_output)); 
csvwrite("q_surf.csv",          q_surf_save(1:n_output));
csvwrite("YO2_surf.csv",        YO2_surf_save(1:n_output));
csvwrite("YO2_nx.csv",          YO2_nx_save(1:n_output));
csvwrite("RR1_peak.csv",        RR1_peak_save(1:n_output));
csvwrite("RR2_peak.csv",        RR2_peak_save(1:n_output));
csvwrite("RR3_peak.csv",        RR3_peak_save(1:n_output));
csvwrite("RR4_peak.csv",        RR4_peak_save(1:n_output));
csvwrite("xRR1_peak.csv",       xRR1_peak_save(1:n_output));
csvwrite("xRR2_peak.csv",       xRR2_peak_save(1:n_output));
csvwrite("xRR3_peak.csv",       xRR3_peak_save(1:n_output));
csvwrite("xRR4_peak.csv",       xRR4_peak_save(1:n_output));
csvwrite("h_conv_surf.csv",     h_conv_save(1:n_output));

csvwrite("cell_center.csv",     XCells_save(:,1:n_output));
csvwrite("Temp.csv",            TEMP_save(:,1:n_output));
csvwrite("x_ws.csv",            XWS_save(:,1:n_output));
csvwrite("x_ds.csv",            XDS_save(:,1:n_output));
csvwrite("x_c.csv",             XC_save(:,1:n_output));
csvwrite("x_a.csv",             XA_save(:,1:n_output));
csvwrite("Y_O2.csv",            YO2G_save(:,1:n_output));
csvwrite("Pres.csv",            PRESG_save(:,1:n_output));
csvwrite("MLRPUV.csv",          MLRPUV_save(:,1:n_output));
csvwrite("MassFlux.csv",        MFLUXG_save(:,1:n_output));
csvwrite("Velocity.csv",        VELG_save(:,1:n_output));
csvwrite("Qdot.csv",            Qdot_save(:,1:n_output));
%}

% Check global conservation statements
fprintf(' \n');
fprintf(' Final thickness of the particle (theory)     = %g \n', ...
                                                                  delta_f);
fprintf(' Final thickness of the particle (simulation) = %g \n', ...
                                                           xRight(nx_new));

%{
% Initial mass of wet solid (t = 0) [kg]
if geometry=="rectangle"
    Mass_ws_i = rho_ws*(1-psi_ws)*2*delta_i*A_rectangle;
elseif geometry=="cylinder"
    Mass_ws_i = rho_ws*(1-psi_ws)*pi*delta_i^2*L_cylinder;
elseif geometry=="sphere"
    Mass_ws_i = rho_ws*(1-psi_ws)*(4/3)*pi*delta_i^3;
end
%}

% Total mass loss at final time
int_MLR = 0;
for n=2:n_output
    int_MLR = int_MLR + 0.5*(MLR_save(n-1)+MLR_save(n)) ...
                           *(time_save(n)-time_save(n-1));
end

% Total solid mass at final time
for i=1:nx_new
    rho_s_bulk(i) = rho_s(i)*(1-psi_sg(i));
end
Mass_s_f = 0;
for i=1:nx_new
	Mass_s_f = Mass_s_f + rho_s_bulk(i)*dV(i);
end
if geometry=="rectangle"
    Mass_ws_i = Mass_ws_i*2;
    Mass_s_f  = Mass_s_f*2;
end
fprintf(' \n');
fprintf(' Total mass of volatiles released (simulation) = %g \n', ...
                                                                  int_MLR);
fprintf(' Total mass of volatiles released (simulation) = %g \n', ...
                                                     (Mass_ws_i-Mass_s_f));

if( ( A_R3 == 0 ) && ( A_R4 == 0 ) )
    Mass_s_f = eta_c_R2*eta_ds_R1*Mass_ws_i;
    fprintf(' Total mass of volatiles released (theory) = %g \n', ...
                                                     (Mass_ws_i-Mass_s_f));
end
if( ( A_R3 ~= 0 ) && ( A_R4 ~= 0 ) )
    Mass_s_f = eta_a_R4*eta_c_R2*eta_ds_R1*Mass_ws_i;
    fprintf(' Total mass of volatiles released (theory) = %g \n', ...
                                                     (Mass_ws_i-Mass_s_f));
end

if geometry=="rectangle"
    mass_ws_save = mass_ws_save*2;
    mass_ds_save = mass_ds_save*2;
    mass_c_save  = mass_c_save*2;
    mass_a_save  = mass_a_save*2;
end


% Plot results
% - Time variations of key quantities

figure(1);   % Mass loss rate [kg/s]
hold on;
plot(time_save(1:n_output),MLR_save(1:n_output));
xlabel('Time (s)');
ylabel('MLR (kg/s)');

figure(2);   % Mass [kg]
hold on;
plot(time_save(1:n_output),mass_ws_save(1:n_output),'--k');
plot(time_save(1:n_output),mass_ds_save(1:n_output),'-k');
plot(time_save(1:n_output),mass_c_save(1:n_output),'-r');
plot(time_save(1:n_output),mass_a_save(1:n_output),'-b');
xlabel('Time (s)');
ylabel('Mass (kg)');

fprintf(' \n');
fprintf(' Total mass of ash at burnout time (theory) = %g \n', ...
                                                                 Mass_s_f);
fprintf(' Total mass of ash at final time (simulation) = %g \n', ...
                                                    mass_a_save(n_output));

figure(3);   % Net surface heat flux [kW/m2]
hold on;
plot(time_save(1:n_output),(q_surf_save(1:n_output)/1000));
xlabel('Time (s)');
ylabel('Net surface heat flux (kW/m2)');

figure(4);   % Temperature at exposed surface and at center of particle [K]
hold on;
plot(time_save(1:n_output),temp_surf_save(1:n_output));
hold on;
plot(time_save(1:n_output),temp_back_save(1:n_output));
xlabel('Time (s)');
ylabel('Temperatures at surface and center (K)');

figure(5);   % Half-thickness or radius of particle
hold on;
plot(time_save(1:n_output),delta_save(1:n_output));
xlabel('Time (s)');
ylabel('Half-thickness or radius (m)');

figure(6);   % Time increment [s]
hold on;
plot(time_save(1:n_output),dt_save(1:n_output));
xlabel('Time (s)');
ylabel('Time increment dt (s)');

figure(7);   % Number of cells [s]
hold on;
plot(time_save(1:n_output),nx_save(1:n_output));
xlabel('Time (s)');
ylabel('Number of cells nx');

figure(8);   % Max variation of Temp, x_k, Y_O2, pres
hold on;
plot(time_save(1:n_output),DeltaQ_max_save(1,1:n_output)','-r');
plot(time_save(1:n_output),DeltaQ_max_save(2,1:n_output)','-k');
plot(time_save(1:n_output),DeltaQ_max_save(3,1:n_output)','-b');
plot(time_save(1:n_output),DeltaQ_max_save(4,1:n_output)','--ok');
xlabel('Time (s)');
ylabel('Max variation during dt');

figure(9); % Y_O2 at x = Delta (surface) and i = nx (center of last cell)
hold on;
plot(time_save(1:n_output),YO2_surf_save(1:n_output)','-r');
plot(time_save(1:n_output),  YO2_nx_save(1:n_output)','-k');
xlabel('Time (s)');
ylabel('YO2');

figure(10); % Peak values of R1-R2-R3-R4 reaction rates
hold on;
plot(time_save(1:n_output),RR1_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),RR2_peak_save(1:n_output)','-b');
plot(time_save(1:n_output),RR3_peak_save(1:n_output)','-r');
plot(time_save(1:n_output),RR4_peak_save(1:n_output)','--or');
xlabel('Time (s)');
ylabel('Peak value of reaction rate (kg/s/m3)');

figure(101); % Peak, surface and back values of R1 reaction rate
hold on;
plot(time_save(1:n_output),RR1_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),RR1_surf_save(1:n_output)','--r');
plot(time_save(1:n_output),RR1_back_save(1:n_output)','--b');
xlabel('Time (s)');
ylabel('Peak, surf, back values of R1 reaction rate (kg/s/m3)');

figure(102); % Peak, surface and back values of R2 reaction rate
hold on;
plot(time_save(1:n_output),RR2_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),RR2_surf_save(1:n_output)','--r');
plot(time_save(1:n_output),RR2_back_save(1:n_output)','--b');
xlabel('Time (s)');
ylabel('Peak, surf, back values of R2 reaction rate (kg/s/m3)');

figure(103); % Peak, surface and back values of R3 reaction rate
hold on;
plot(time_save(1:n_output),RR3_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),RR3_surf_save(1:n_output)','--r');
plot(time_save(1:n_output),RR3_back_save(1:n_output)','--b');
xlabel('Time (s)');
ylabel('Peak, surf, back values of R3 reaction rate (kg/s/m3)');

figure(104); % Peak, surface and back values of R4 reaction rate
hold on;
plot(time_save(1:n_output),RR4_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),RR4_surf_save(1:n_output)','--r');
plot(time_save(1:n_output),RR4_back_save(1:n_output)','--b');
xlabel('Time (s)');
ylabel('Peak, surf, back values of R4 reaction rate (kg/s/m3)');

figure(11); % Location of peak values of R1-R2-R3-R4 reaction rates
hold on;
plot(time_save(1:n_output),xRR1_peak_save(1:n_output)','-k');
plot(time_save(1:n_output),xRR2_peak_save(1:n_output)','-b');
plot(time_save(1:n_output),xRR3_peak_save(1:n_output)','-r');
plot(time_save(1:n_output),xRR4_peak_save(1:n_output)','--or');
xlabel('Time (s)');
ylabel('Location of peak value of reaction rate (m)');

figure(12); % Convective heat transfer coefficient [W/m2/K]
hold on;
plot(time_save(1:n_output),h_conv_save(1:n_output));
xlabel('Time (s)');
ylabel('Convective heat transfer coefficient (W/m2/K)');

% - Spatial profiles

figure(21);   % Temperature
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),TEMP_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Temperature (K)');

figure(22);   % Solid mass density (ms/Vs)
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RHOS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Solid mass density (kg/m3)');

figure(23);   % Porosity
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),PSI_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Porosity');

figure(24);   % Volume fraction of wet solid
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XWS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of wet solid');

figure(25);   % Volume fraction of dry solid
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XDS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of dry solid');

figure(26);   % Volume fraction of char
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XC_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of char');

figure(27);   % Volume fraction of ash
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XA_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of ash');

figure(28);   % Check
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n), (1 ...
        -XWS_save(1:nx_save(n),n)-XDS_save(1:nx_save(n),n) ...
        -XC_save (1:nx_save(n),n)-XA_save (1:nx_save(n),n)) );
    hold on
end
xlabel('Spatial distance (m)');
ylabel('1-Sumxk (Check)');

figure(29);   % Mass reaction rate for moisture evaporation reaction R1
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR1_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR1 (kg/s/m3)');

figure(30);   % Mass reaction rate for thermal pyrolysis reaction R2
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR2_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR2 (kg/s/m3)');

figure(31);   % Mass reaction rate for oxidative pyrolysis reaction R3
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR3_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR3 (kg/s/m3)');

figure(32);   % Mass reaction rate for char oxidation reaction R4
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR4_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR4 (kg/s/m3)');

figure(33);   % Mass fraction of gaseous oxygen
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),YO2G_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Mass fraction of gaseous oxygen');

figure(34);   % Gauge pressure
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),PRESG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Gauge pressure (Pa)');

figure(35);   % Volumetric mass loss rate
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),MLRPUV_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('MLRPUV (kg/s/m3)');

figure(36);   % Volumetric HRR
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),Qdot_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('HRRPUV (kg/s/m3)');

figure(37);   % Gas Mass Flux
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),MFLUXG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Mass Flux (kg/s/m2)');

figure(38);   % Gas Velocity
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),VELG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Velocity (m/s)');

%{
n=n_output
plot(XCells_save(1:nx_save(n),n),TEMP_save(1:nx_save(n),n),'-k','LineWidth',3);
xlabel('Spatial distance (m)','FontSize',16);
ylabel('Temperature (K)','FontSize',16);

hold on;
plot(XCells_save(1:nx_save(n),n),XWS_save(1:nx_save(n),n),'--k','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),XDS_save(1:nx_save(n),n),'-k','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),XC_save(1:nx_save(n),n),'-r','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),XA_save(1:nx_save(n),n),'-b','LineWidth',3);
xlabel('Spatial distance (m)','FontSize',16);
ylabel('Volume fraction of solid species','FontSize',16);

hold on;
plot(XCells_save(1:nx_save(n),n),RR1_save(1:nx_save(n),n),'-k','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),RR2_save(1:nx_save(n),n),'-r','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),RR3_save(1:nx_save(n),n),'--r','LineWidth',3);
plot(XCells_save(1:nx_save(n),n),RR4_save(1:nx_save(n),n),'-','Color','#D95319','LineWidth',3);
xlabel('Spatial distance (m)','FontSize',16);
ylabel('Mass reaction rates kg/s/m3)','FontSize',16);

plot(time_save(1:n_output),(0.5*MLR_save(1:n_output)),'-k','LineWidth',3);
xlabel('Time (s)','FontSize',16);
ylabel('MLR (kg/s/m2)','FontSize',16);
%}

%AT
% Save
% ----
FMC    = (1-eta_ds_R1)/eta_ds_R1;        % Moisture content of the particle
x_g_O2 = (Y_g_O2/32)/( (Y_g_O2/32)+((1-Y_g_O2)/28) ); % x_O2 in ambient gas
yy     = [ delta_i; (G_i/1000); u_g; x_g_O2; T_g; FMC ];

fid = fopen('PBR_Case.txt','w');
fprintf(fid,['%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e \n'],yy);

for n=1:n_output
    zz = [ time_save(n); MLR_save(n); ...
                         MLR_R1_save(n); ...
                         MLR_R2_save(n); ...
                         MLR_R3_save(n); ...
                         MLR_R4_save(n); ...
                         RR1_peak_save(n); ...
                         RR2_peak_save(n); ...
                         RR3_peak_save(n); ...
                         RR4_peak_save(n); ...
                         temp_xRR1_peak_save(n); ...
                         temp_xRR2_peak_save(n); ...
                         temp_xRR3_peak_save(n); ...
                         temp_xRR4_peak_save(n); ...
                         xRR1_peak_save(n); ...
                         xRR2_peak_save(n); ...
                         xRR3_peak_save(n); ...
                         xRR4_peak_save(n); ...
                         RR1_surf_save(n); ...
                         RR1_back_save(n); ...
                         RR2_surf_save(n); ...
                         RR2_back_save(n); ...
                         RR3_surf_save(n); ...
                         RR3_back_save(n); ...
                         RR4_surf_save(n); ...
                         RR4_back_save(n); ...
                         mass_ws_save(n); ...
                         mass_ds_save(n); ...
                         mass_c_save(n); ...
                         mass_a_save(n); ...
                         (q_surf_save(n)/1000); ...
                         h_conv_save(n); ...
                         temp_surf_save(n); ...
                         temp_back_save(n); ...
                         temp_max_save(n); ...
                         YO2_surf_save(n); ...
                         YO2_back_save(n); ...
                         YO2_min_save(n); ...
                         pres_surf_save(n); ...
                         pres_back_save(n); ...
                         pres_max_save(n); ...
                         MLR_surf_save(n); ...
                         vel_surf_save(n); ...
                         vel_back_save(n); ...
                         vel_max_save(n) ];

    fprintf(fid,[ ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ' ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ' ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ' ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ' ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e ' ...
     '%12.8e %12.8e %12.8e %12.8e %12.8e \n'],zz);
end

fclose(fid);
%AT


return;
