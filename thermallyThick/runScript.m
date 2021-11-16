% Particle Burning Rate Model - PBR_model.m

clc;
clear;
close all;

global geometry A_rectangle L_cylinder ...
       T_g u_g G Y_g_O2 pres_g k_g0 cp_g0 nu_g0 MW_g R sigma ...
       rho_ws rho_ds rho_c rho_a ...
       k_ws k_ds k_c k_a ...
       c_ws c_ds c_c c_a ...
       eps_ws eps_ds eps_c eps_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 DeltaH_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 DeltaH_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 DeltaH_R3 eta_c_R3 eta_O2_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 DeltaH_R4 eta_a_R4 eta_O2_R4 ...
       dx_i 

% Read input parameters
[ T_end, geometry, delta_i, A_rectangle, L_cylinder, ...
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
         = input_parameters;

% NOTE: delta_i = Initial half-thickness (geometry = "rectangle") or
%                 initial radius (geometry = "cylinder" or "sphere")
%       Select spatial resolution ~ 100 microns
dx0  = 100e-6;               % Spatial resolution [m]
nx_i = round(delta_i/dx0);   % Initial number of grid cells
dx_i = delta_i/nx_i;              % Initial grid cell size
fprintf(' dx_i = %g \n',dx_i);
fprintf(' nx_i = %g \n',nx_i);
fprintf(' \n');

% Calculate initial time step
dt_i = timeStep(rho_ws, rho_ds, rho_c, rho_a, ...
                k_ws, k_ds, k_c, k_a, ...
                c_ws, c_ds, c_c, c_a, ...
                psi_ws, psi_ds, psi_c, psi_a, ...
                Kperm_ws, Kperm_ds, Kperm_c, Kperm_a, ...
                A_R1, Ta_R1, ...
                A_R2, Ta_R2, ...
                A_R3, Ta_R3, n_O2_R3, ...
                A_R4, Ta_R4, n_O2_R4, ...
                Y_g_O2, nu_g0, MW_g, R, dx_i);
%n_f = round(T_end/dt_i);   % Total number of time steps
n_f = 1e7;   % Total number of time steps
fprintf(' n_f = %g \n',n_f);
fprintf(' \n');
pause;

% Parameter that controls the frequency at which quantities are saved
% to output files
%%AT i_output  = 50;
i_output  = 500;
%%AT

% Parameters that control the time step and solution accuracy
%%AT Threshold_Temp = 0.01;    % Max. value of temp variation during dt [K]
%%AT Threshold_xk   = 0.001;   % Max. value of x_k variation during dt [-]
%%AT Threshold_YO2  = 0.01;    % Max. value of Y_O2 variation during dt [-]
%%AT Threshold_pres = 0.1;     % Max. value of pres variation during dt [-]
Threshold_Temp = 0.1;     % Max. value of temp variation during dt [K]
Threshold_xk   = 0.01;    % Max. value of x_k variation during dt [-]
Threshold_YO2  = 0.01;    % Max. value of Y_O2 variation during dt [-]
Threshold_pres = 0.1;     % Max. value of pres variation during dt [-]
lambda         = 5.0;     % Under-relaxation parameter fore temp & Y_O2 [-]
lambda_pres    = 5.0;     % Under-relaxation parameter for pressure [-]
absTol_Temp    = Threshold_Temp/10;    % Abs tolerance in iteration loop for temp [-]
absTol_YO2     = Threshold_YO2/10;     % Abs. tolerance in iteration loop for Y_O2 [-]
absTol_pres    = Threshold_pres/10;    % Abs. tolerance in iteration loop for pres [Pa]
%% AT
dt_max         = 0.1;     % Max. value of dt [s]

% Set Initial conditions (t = 0)
n         = 0;
time      = 0; 
n_output  = 0;

nx_new    = nx_i;

temp_surf = temp_i;    % Surface temperature [K] (t = 0)

% Calculation of convective heat transfer coefficient
T_film       = 0.5*(temp_surf+T_g);      % Film temperature [K]
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

q_surf_i = eps_ws*G - eps_ws*sigma*temp_surf^4 ...
         + h_conv*(T_g-temp_surf);  % Net surface heat flux [W/m2] (t = 0)
q_surf   = q_surf_i;

mdot_surf_i = 0;
mdot_surf   = mdot_surf_i;   % Surface oxygen mass flux [kg/s] (t = 0)

% Solid mass density (t = 0)
rho_s_i = rho_ws*x_ws_i + rho_ds*x_ds_i + rho_c*x_c_i + rho_a*x_a_i;

% Porosity (t = 0)
psi_i   = psi_ws*x_ws_i + psi_ds*x_ds_i + psi_c*x_c_i + psi_a*x_a_i;

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

XCells_save =         ones(nx_i,round(n_f/i_output));  % Coord. cell ctrs.
TEMP_save   = temp_i *ones(nx_i,round(n_f/i_output));  % Solid temperature
RHOS_save   = rho_s_i*ones(nx_i,round(n_f/i_output));  % Solid mass density
PSI_save    = psi_i  *ones(nx_i,round(n_f/i_output));  % Porosity
XWS_save    = x_ws_i *ones(nx_i,round(n_f/i_output));  % Vol. frac. wet s.
XDS_save    = x_ds_i *ones(nx_i,round(n_f/i_output));  % Vol. frac. dry s.
XC_save     = x_c_i  *ones(nx_i,round(n_f/i_output));  % Vol. frac. char
XA_save     = x_a_i  *ones(nx_i,round(n_f/i_output));  % Vol. frac. ash
RR1_save    =         zeros(nx_i,round(n_f/i_output)); % Reaction rate R1
RR2_save    =         zeros(nx_i,round(n_f/i_output)); % Reaction rate R2
RR3_save    =         zeros(nx_i,round(n_f/i_output)); % Reaction rate R3
RR4_save    =         zeros(nx_i,round(n_f/i_output)); % Reaction rate R4
Qdot_save   =         zeros(nx_i,round(n_f/i_output)); % Volumetric HRR
YO2G_save   = Y_O2_i *ones(nx_i,round(n_f/i_output));  % Mass frac. of O2
PRESG_save  = pres_i *ones(nx_i,round(n_f/i_output));  % Pressure
MLRPUV_save =         zeros(nx_i,round(n_f/i_output)); % Volumetric MLR
MFLUXG_save =         zeros(nx_i,round(n_f/i_output)); % Gas mass flux
VELG_save   =         zeros(nx_i,round(n_f/i_output)); % Gas velocity

% Set the computational grid
[xRight, xCenter, xLeft1, dV] = mesh(nx_i, dx);

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

fid = fopen('LogFile.txt','w');
%--- Begin Time loop ---
while ((n ~= n_f) && (flag_burnout ~= 2) && (time < T_end))
    
    n = n + 1;
    
    % Solution at previous time step (t = time-dt)
    nx_old    = nx_new;
    dV_old    = dV;
    dx_old    = dx;
    dt_old    = dt;    
    temp_old  = temp;
    x_ws_old  = x_ws;
    x_ds_old  = x_ds;
    x_c_old   = x_c;
    x_a_old   = x_a;
    Y_O2_old  = Y_O2;
    pres_old  = pres;
    
    q_surf_old    = q_surf;
    temp_surf_old = temp_surf;
    mdot_surf_old = mdot_surf;
    
    % Solution at new time step (t = time)
    
    time = time + dt;
    if(mod(n,i_output)==0)
        n_output            = n_output+1;
        time_save(n_output) = time;
        fprintf(' \n');
        fprintf(' n_output = %g, time = %g \n', ...
                  n_output,time_save(n_output));
    end
    
    flag_iter = 0;
    
    %--- Begin iterative loop ---    
    % Initial guess
    temp_newiter = temp_old;
    temp1        = temp_old;
    x_ws_newiter = x_ws_old;
    x_ds_newiter = x_ds_old;
    x_c_newiter  = x_c_old;
    x_a_newiter  = x_a_old;
    Y_O2_newiter = Y_O2_old;
    pres_newiter = pres_old;
        
    for iter=1:1000
        
    % Solution at new time step, at previous iteration
    temp_olditer = temp_newiter;
    x_ws_olditer = x_ws_newiter;
    x_ds_olditer = x_ds_newiter;
    x_c_olditer  = x_c_newiter;
    x_a_olditer  = x_a_newiter;
    Y_O2_olditer = Y_O2_newiter;
    pres_olditer = pres_newiter;
        
    % Calculate solution at new time step, at new iteration
        
    % - Calculate temperature
    %
    % Note: solve dTp/dt = RHS = (RHS1+RHS2) where RHS1 is energy change
    % due to reactions (R1)-(R4) and RHS2 is energy change due to
    % diffusive/convective transport of heat and to gas-to-solid heat
    % transfer at the exposed surface
    %
    % Discretization (implicit):
    %   (Tp(n+r1)-Tp(n+r0))/dt = RHS1 where r0 and r1 designate fractional
    %                            steps, 0 <= r0,r1 <= 1 and where RHS1 is
    %                            evaluated at time (n+r1)
    %   (Tp(n+r1)-Tp(n+r0))/dt = RHS2 where RHS2 is evaluated using
    %                            Crank-Nicolson
    %
    % Evaluation of RHS2 (energy change due to transport and heat transfer
    % at the exposed surface):
    %   RHS2 = 0.5*RHS2(n+r1) + 0.5*RHS2(n+r0)
    %
    % Operator splitting method (Strang's method)
    %   Step 1: 0.25*d(Tp)/dt = (0.5*RHS1)
    %           t(n) <= t <= t(n+0.25)
    %   Step 2: 0.5 *d(Tp)/dt = (1.0*RHS2)
    %           t(n+0.25) <= t <= t(n+0.75)
    %   Step 3: 0.25*d(Tp)/dt = (0.5*RHS1)
    %           t(n+0.75) <= t <= t(n+1)
    %
        
    % Step 1: calculate temp1 as the temperature at sub-step t(n+0.25)
    %   tempr0_newiter = temp_old;
    %   tempr1_olditer = temp1;
    [temp1] = ...
        energy_conservation_reactionstep(dt, temp_old, temp1, ...
    x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, Y_O2_olditer, ...
    lambda, temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, nx_old);
    temp1   = temp1';
    
    % Step 2: calculate temp2 as the temperature at sub-step t(n+0.75)
    %    tempr0_newiter = temp1
    tempc_nx = temp_old(nx_old);
    
    [a, b, c, d] = ...
              energy_conservation_diffusionstep(dt, temp1, ...
              temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
              q_surf_old, tempc_nx, temp_surf_old, h_conv, ...
              xRight, xCenter, dV_old, nx_old);
    temp2 = tri(a,b,c,d);
    temp2 = temp2';
    
    % Step 3: calculate temp3 as the temperature at sub-step t(n+1)
    %   tempr0_newiter = temp2;
    %   tempr1_olditer = temp_olditer;
    [temp_newiter] = ...
        energy_conservation_reactionstep(dt, temp2, temp_olditer, ...
    x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, Y_O2_olditer, ...
    lambda, temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, nx_old);
    temp_newiter   = temp_newiter';
        
    % Apply under-relaxation;
    %%AT temp_newiter = lambda*temp_newiter + (1-lambda)*temp_olditer;
    %%AT
        
    % Check for convergence of iteration loop
    DeltaTemp_max = max( abs(real(temp_newiter-temp_olditer)) );
    if (mod(n,i_output)==0)
        fprintf(' DeltaTemp_max,iter = %g %g \n', ...
                  DeltaTemp_max,iter);
    end
    if (iter > 1)
        fprintf(fid,' Time,DeltaTemp_max,iter = %g %g %g \n', ...
                      time,DeltaTemp_max,iter);
        flag_iter = 1;
    end
    
    % - Calculate solid species volume fractions and cell volumes
    [x_ws_newiter, x_ds_newiter, x_c_newiter, x_a_newiter, dV] = ...       
    mass_conservation(dt, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, dV_old, nx_old);
    x_ws_newiter = x_ws_newiter';
    x_ds_newiter = x_ds_newiter';
    x_c_newiter  = x_c_newiter';
    x_a_newiter  = x_a_newiter';
    
    %%AT if( ( A_R3 == 0 ) && ( A_R4 == 0 ) )   % (R1)-(R2) reaction model
    if( 1 == 0 ) % Uncomment this line to force calc. of Y_O2 & pres
        Y_O2_newiter  = Y_O2_old;
        pres_newiter  = pres_old;
        DeltaYO2_max  = 0;
        DeltaPres_max = 0;
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
    Y_O2_olditer, lambda, mdot_surf_old, h_conv, ...
    xRight, xCenter, dV_old, nx_old);
    Y_O2_newiter = tri(a,b,c,d);
    Y_O2_newiter = Y_O2_newiter';
        
    % Apply under-relaxation
    %%AT Y_O2_newiter = lambda*Y_O2_newiter + (1-lambda)*Y_O2_olditer;
    %%AT
        
    % Check for convergence of iteration loop
    DeltaYO2_max = max( abs(real(Y_O2_newiter-Y_O2_olditer)) );
    if(mod(n,i_output)==0)
        fprintf(' DeltaYO2_max,iter  = %g %g \n', ...
                  DeltaYO2_max,iter);
    end
    if (iter > 1)
        fprintf(fid,' Time,DeltaYO2_max,iter  = %g %g %g \n',...
                      time,DeltaYO2_max,iter);
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
    %%AT
    lambda_pres = 0;
    %%AT
    [a, b, c, d] = pressure_equation(dt, pres_old, ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, pres_olditer, lambda_pres, ...
    xRight, xCenter, dV_old, nx_old);
    pres_newiter = tri(a,b,c,d);
    pres_newiter = pres_newiter';
        
    % Apply under-relaxation
    %%AT pres_newiter = lambda_pres*pres_newiter + (1-lambda_pres)*pres_olditer;
    pres_newiter = 0.5*pres_newiter + (1-0.5)*pres_olditer;
    %%AT 
        
    % Check for convergence of iteration loop
    DeltaPres_max = max( abs(real(pres_newiter-pres_olditer)) );
    if(mod(n,i_output)==0)
        fprintf(' DeltaPres_max,iter = %g %g \n', ...
                  DeltaPres_max,iter);
    end
    if (iter > 1)
        fprintf(fid,' Time,DeltaPres_max,iter = %g %g %g \n',...
                      time,DeltaPres_max,iter);
    end
    
        % *** End calculation of Y_O2 and pres ***
    end
    
    if( ( DeltaTemp_max < (0.1*absTol_Temp/lambda) ) & ...
        ( DeltaYO2_max  < (0.1*absTol_YO2 /lambda) ) & ...
        ( DeltaPres_max < (0.1*absTol_pres/lambda) ) )
        break;
    end
    if(iter==1000)
        fprintf(' Problem in iteration loop: break at time %g \n',time);
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
    % - Update dt such that abs( Tp(n+1)-Tp(n) ) <= Threshold     
    DeltaTemp_max = max( abs(real(temp-temp_old)) );
    
    % Safety: (1) Limit the time step dt so that variations in temperature
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( DeltaTemp_max > 0 )
        dt_Temp = dt_old*Threshold_Temp/DeltaTemp_max;     % Constraint (1)
        dt_Temp = max(0.9*dt_old,min(1.1*dt_old,dt_Temp)); % Constraint (2)
    else
        dt_Temp = 1.1*dt_old;
    end

    % Time step restriction:
    % - Update dt such that abs( xk(n+1)-xk(n) ) <= Threshold
    Deltax_ws_max = max( abs(real(x_ws-x_ws_old)) );
    Deltax_ds_max = max( abs(real(x_ds-x_ds_old)) );
    Deltax_c_max  = max( abs(real(x_c -x_c_old) ) );
    Deltax_a_max  = max( abs(real(x_a -x_a_old) ) );
    Deltax_k_max  = max( [ Deltax_ws_max,Deltax_ds_max, ...
                           Deltax_c_max,Deltax_a_max ] );
    % Safety: (1) Limit the time step dt so that variations in x_k
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( Deltax_k_max > 0 )
        dt_xk = dt_old*Threshold_xk/Deltax_k_max;        % Constraint (1)
        dt_xk = max(0.9*dt_old,min(1.1*dt_old,dt_xk));   % Constraint (2)
    else
        dt_xk = 1.1*dt_old;
    end
     
    % Time step restriction:
    % - Update dt such that abs( Y_O2(n+1)-Y_O2(n) ) <= Threshold     
    DeltaYO2_max = max( abs(real(Y_O2-Y_O2_old)) );
    
    % Safety: (1) Limit the time step dt so that variations in Y_O2
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( DeltaYO2_max > 0 )
        dt_YO2 = dt_old*Threshold_YO2/DeltaYO2_max;      % Constraint (1)
        dt_YO2 = max(0.9*dt_old,min(1.1*dt_old,dt_YO2)); % Constraint (2)
    else
        dt_YO2 = 1.1*dt_old;
    end

    % Time step restriction:
    % - Update dt such that abs( pres(n+1)-pres(n) ) <= Threshold
    DeltaPres_max = max( abs(real(pres-pres_old)) );

    % Safety: (1) Limit the time step dt so that variations in pres
    %             are less than a user-defined Threshold
    %         (2) Limit the time step dt to a maximum of 10% variations
    if( DeltaPres_max > 0 )
        dt_pres = dt_old*Threshold_pres/DeltaPres_max;     % Constraint (1)
        dt_pres = max(0.9*dt_old,min(1.1*dt_old,dt_pres)); % Constraint (2)
    else
        dt_pres = 1.1*dt_old;
    end
    
    dt = min( [ dt_Temp,dt_xk,dt_YO2,dt_pres ] );
    %%AT dt = min( [ dt_Temp,dt_xk,dt_YO2 ] );
    dt = min(dt,dt_max); % Safety: do not let dt go to very large values
    %%AT
    %%AT dt = dt_old; % Uncomment this line for tests with fixed dt
    %%AT 
 
    if(mod(n,i_output)==0)
        fprintf(' dt                   = %g \n',dt);
        fprintf(' iter                 = %g \n',iter);
        fprintf(' max(|temp-temp_old|) = %g \n',DeltaTemp_max);
        fprintf(' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max);
        fprintf(' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max);
        fprintf(' max(|pres-pres_old|) = %g \n',DeltaPres_max);
    end        
    if( ( flag_iter == 1 ) & ( iter >= 50 ) )
        fprintf(' *** WARNING *** \n');
        fprintf(' flag_iter = %g \n', flag_iter);
        fprintf(' iter                 = %g \n',iter);
        fprintf(' time, dt             = %g %g \n',time,dt);
        fprintf(' max(|temp-temp_old|) = %g \n',DeltaTemp_max);
        fprintf(' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max);
        fprintf(' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max);
        fprintf(' max(|pres-pres_old|) = %g \n',DeltaPres_max);
    end
    if( ( flag_iter == 1 ) & ( iter > 1 ) )
        fprintf(fid,' *** WARNING *** \n');
        fprintf(fid,' flag_iter = %g \n', flag_iter);
        fprintf(fid,' iter                 = %g \n',iter);
        fprintf(fid,' time, dt             = %g %g \n',time,dt);
        fprintf(fid,' max(|temp-temp_old|) = %g \n',DeltaTemp_max);
        fprintf(fid,' max(|x_k-x_k_old|)   = %g \n',Deltax_k_max);
        fprintf(fid,' max(|YO2-YO2_old|)   = %g \n',DeltaYO2_max);
        fprintf(fid,' max(|pres-pres_old|) = %g \n',DeltaPres_max);
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
    temp_old = temp;
    x_ws_old = x_ws;
    x_ds_old = x_ds;
    x_c_old  = x_c;
    x_a_old  = x_a;
    Y_O2_old = Y_O2;
    pres_old = pres;
    
    [temp, x_ws, x_ds, x_c, x_a, Y_O2, pres] = interpolateOnNewMesh( ...
              nx_new, xCenter, volume, dV, nx_old, xCenter_old, dV_old, ...
       temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, Y_O2_old, pres_old);
    temp = temp';
    x_ws = x_ws';
    x_ds = x_ds';
    x_c  = x_c';
    x_a  = x_a';
    Y_O2 = Y_O2';
    pres = pres';
    % End remeshing
    
    for i=1:nx_new
        rho_s(i)  = rho_ws*x_ws(i) + rho_ds*x_ds(i) ...
                  + rho_c*x_c(i)+ rho_a*x_a(i);
        psi_sg(i) = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
                  + psi_c*x_c(i) + psi_a*x_a(i);
    end
     
    % Calculate the net surface heat flux and the surface temperature
    
    % Porosity
    psi       = psi_sg(nx_new);
    % Effective thermal conductivity
    %  - Porous medium treatment: 
    %    keff = (1-psi) x k_s + psi x k_g
    k_s       = k_ws*x_ws(nx_new) + k_ds*x_ds(nx_new) ...
              + k_c*x_c(nx_new) + k_a*x_a(nx_new); % Conduct. in sol. phase
    k_g       = k_g0*(temp(nx_new)/300)^(0.76);    % Conduct. in gas phase
    keff_surf = (1-psi)*k_s + psi*k_g;
    eps_surf  = eps_ws*x_ws(nx_new) + eps_ds*x_ds(nx_new) ...
              + eps_c*x_c(nx_new) + eps_a*x_a(nx_new);
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
        rho_surf = pres_g*MW_g/R/temp_surf;
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
    
    [q_surf, temp_surf] = particle_surface(temp_surf_old, ...
                           keff_surf, eps_surf, dx_surf, tempc_nx, h_conv);
    
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
            K_R1  = ( (rho_ws*x_ws(i)*(1-psi))^n_R1 ) ...
                    *A_R1*exp(-Ta_R1/temp(i));
            K_R2  = ( (rho_ds*x_ds(i)*(1-psi))^n_R2 ) ...
                    *A_R2*exp(-Ta_R2/temp(i));
            K_R3  = ( (rho_ds*x_ds(i)*(1-psi))^n_R3 ) ...
                    *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp(i));
            K_R4  = ( (rho_c *x_c(i) *(1-psi))^n_R4 ) ...
                    *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp(i));
      
            MLRPUV(i) = (1-eta_ds_R1)*K_R1  + (1-eta_c_R2) *K_R2 ...
                      + (1-eta_c_R3) *K_R3  + (1-eta_a_R4) *K_R4;

            HRRPUV(i) = (1-eta_ds_R1)*K_R1*DeltaH_R1 ...
                      + (1-eta_c_R2) *K_R2*DeltaH_R2 ...
                      + (1-eta_c_R3) *K_R3*DeltaH_R3 ...
                      + (1-eta_a_R4) *K_R4*DeltaH_R4;
        end
        
        mdotPUA(1) = 0;
        % Permeability
        Kperm = Kperm_ws*x_ws(1) + Kperm_ds*x_ds(1) ...
              + Kperm_c*x_c(1) + Kperm_a*x_a(1);
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
                  + Kperm_c*x_c(i) + Kperm_a*x_a(i);   
            % Gas kinematic viscosity
            nu_g  = nu_g0*(temp(i)/300)^(1.76);
            
            if( (pres(i-1) > pres(i)) & (pres(i) > pres(i+1)) )
                % Case for which mdotPUA > 0
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
              + Kperm_c*x_c(nx_new) + Kperm_a*x_a(nx_new);
        % Gas kinematic viscosity
        nu_g  = nu_g0*(temp(nx_new)/300)^(1.76);
        if( pres_old(nx_new-1) > pres_old(nx_new) )
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
            MLRtot_new = MLRtot_new*2 ; % We solve for half of the particle
        end
        
        DeltaQ_max_save(1,n_output) = DeltaTemp_max;
        DeltaQ_max_save(2,n_output) = Deltax_k_max;
        DeltaQ_max_save(3,n_output) = DeltaYO2_max;
        DeltaQ_max_save(4,n_output) = DeltaPres_max;
         
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
         
        for i=1:nx_new
            Y_O2s = max(Y_O2(i),0);
            K_R1  = ( (rho_ws*x_ws(i)*(1-psi))^n_R1 ) ...
                    *A_R1*exp(-Ta_R1/temp(i));
            K_R2  = ( (rho_ds*x_ds(i)*(1-psi))^n_R2 ) ...
                    *A_R2*exp(-Ta_R2/temp(i));
            K_R3  = ( (rho_ds*x_ds(i)*(1-psi))^n_R3 ) ...
                    *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp(i));
            K_R4  = ( (rho_c *x_c(i) *(1-psi))^n_R4 ) ...
                    *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp(i));
                    
            RR1_save(i,n_output) = K_R1;
            RR2_save(i,n_output) = K_R2;
            RR3_save(i,n_output) = K_R3;
            RR4_save(i,n_output) = K_R4;   
        end
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

%% Output data
%{
csvwrite("number_cells.csv",    nx_save(1:n_output));
csvwrite("time.csv",            time_save(1:n_output));
csvwrite("MLR.csv",             MLR_save(1:n_output));
csvwrite("delta.csv",           delta_save(1:n_output));
csvwrite("Temp_surf.csv",       temp_surf_save(1:n_output));
csvwrite("Temp_core.csv",       temp_back_save(1:n_output)); 
csvwrite("q_surf.csv",          q_surf_save(1:n_output));
csvwrite("YO2_surf.csv",        YO2_surf_save(1:n_output));
csvwrite("YO2_nx.csv",          YO2_nx_save(1:n_output));

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

% Initial mass of wet solid (t = 0) [kg]
if geometry=="rectangle"
    Mass_ws_i = rho_ws*(1-psi_ws)*2*delta_i*A_rectangle;
elseif geometry=="cylinder"
    Mass_ws_i = rho_ws*(1-psi_ws)*pi*delta_i^2*L_cylinder;
elseif geometry=="sphere"
    Mass_ws_i = rho_ws*(1-psi_ws)*(4/3)*pi*delta_i^3;
end

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
    Mass_s_f = Mass_s_f*2;
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
if( ( A_R3 == 0 ) && ( A_R4 ~= 0 ) )
    Mass_s_f = eta_a_R4*eta_c_R2*eta_ds_R1*Mass_ws_i;
    fprintf(' Total mass of volatiles released (theory) = %g \n', ...
                                                     (Mass_ws_i-Mass_s_f));
end


% Plot results
% - Time variations of key quantities

figure(1);   % Mass loss rate [kg/s/m2]
hold on;
plot(time_save(1:n_output),(MLR_save(1:n_output)));
xlabel('Time (s)');
ylabel('MLR (kg/s)');

figure(2);   % Net surface heat flux [kW/m2]
hold on;
plot(time_save(1:n_output),(q_surf_save(1:n_output)/1000));
xlabel('Time (s)');
ylabel('Net surface heat flux (kW/m2)');

figure(3);   % Temperature at exposed surface and at center of particle [K]
hold on;
plot(time_save(1:n_output),temp_surf_save(1:n_output));
hold on;
plot(time_save(1:n_output),temp_back_save(1:n_output));
xlabel('Time (s)');
ylabel('Temperatures at surface and center (K)');

figure(4);   % Half-thickness or radius of particle
hold on;
plot(time_save(1:n_output),delta_save(1:n_output));
xlabel('Time (s)');
ylabel('Half-thickness or radius (m)');

figure(5);   % Time increment [s]
hold on;
plot(time_save(1:n_output),dt_save(1:n_output));
xlabel('Time (s)');
ylabel('Time increment dt (s)');

figure(6);   % Number of cells [s]
hold on;
plot(time_save(1:n_output),nx_save(1:n_output));
xlabel('Time (s)');
ylabel('Number of cells nx');

figure(7);   % Max variation of Temp, x_k, Y_O2, pres
hold on;
plot(time_save(1:n_output),DeltaQ_max_save(1,1:n_output)','-r');
plot(time_save(1:n_output),DeltaQ_max_save(2,1:n_output)','-k');
plot(time_save(1:n_output),DeltaQ_max_save(3,1:n_output)','-b');
plot(time_save(1:n_output),DeltaQ_max_save(4,1:n_output)','--ok');
xlabel('Time (s)');
ylabel('Max variation during dt');

figure(8); % Y_O2 at x = Delta (surface) and i = nx (center of last cell)
hold on;
plot(time_save(1:n_output),YO2_surf_save(1:n_output)','-r');
plot(time_save(1:n_output),  YO2_nx_save(1:n_output)','-k');
xlabel('Time (s)');
ylabel('YO2');

% - Spatial profiles

figure(11);   % Temperature
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),TEMP_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Temperature (K)');

figure(12);   % Solid mass density (ms/Vs)
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RHOS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Solid mass density (kg/m3)');

figure(13);   % Porosity
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),PSI_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Porosity');

figure(14);   % Volume fraction of wet solid
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XWS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of wet solid');

figure(15);   % Volume fraction of dry solid
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XDS_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of dry solid');

figure(16);   % Volume fraction of char
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XC_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of char');

figure(17);   % Volume fraction of ash
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),XA_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Volume fraction of ash');

figure(18);   % Check
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n), (1 ...
        -XWS_save(1:nx_save(n),n)-XDS_save(1:nx_save(n),n) ...
        -XC_save (1:nx_save(n),n)-XA_save (1:nx_save(n),n)) );
    hold on
end
xlabel('Spatial distance (m)');
ylabel('1-Sumxk (Check)');

figure(19);   % Mass reaction rate for moisture evaporation reaction R1
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR1_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR1 (kg/s/m3)');

figure(20);   % Mass reaction rate for thermal pyrolysis reaction R2
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR2_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR2 (kg/s/m3)');

figure(21);   % Mass reaction rate for oxidative pyrolysis reaction R3
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR3_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR3 (kg/s/m3)');

figure(22);   % Mass reaction rate for char oxidation reaction R4
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),RR4_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('RR4 (kg/s/m3)');

figure(23);   % Mass fraction of gaseous oxygen
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),YO2G_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Mass fraction of gaseous oxygen');

figure(24);   % Gauge pressure
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),PRESG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Gauge pressure (Pa)');

figure(25);   % Volumetric mass loss rate
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),MLRPUV_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('MLRPUV (kg/s/m3)');

figure(26);   % Volumetric HRR
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),Qdot_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('HRRPUV (kg/s/m3)');

figure(27);   % Gas Mass Flux
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),MFLUXG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Mass Flux (kg/s/m2)');

figure(28);   % Gas Velocity
for n=1:n_output
    plot(XCells_save(1:nx_save(n),n),VELG_save(1:nx_save(n),n));
    hold on
end
xlabel('Spatial distance (m)');
ylabel('Velocity (m/s)');


return;
