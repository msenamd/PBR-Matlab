% Particle Burning Rate Model - PBR_model.m

clc;
clear;
close all;

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       k_m k_vs k_c c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g ...
       T_g G eps_m eps_vs eps_c dx_i geometry 

% Read input parameters
[T_g, u_g, G, x_O2_g, delta_i, FMC, ...
 rho_m, rho_vs, rho_c, k_m, k_vs, k_c, c_m, c_vs, c_c, ...
 eps_m, eps_vs, eps_c, A_R1, Ta_R1, DeltaH_R1, ...
 A_R2, Ta_R2, DeltaH_R2, eta_c, A_R3, Ta_R3, DeltaH_R3,...
 temp_surf_i,T_end,geometry,A_rectangle,L_cylinder] = input_parameters;

% NOTE: delta_i = Initial half-thickness for rectagular particles 
% or intial radius for cylindrical and spherical particles
% Select numerical parameters: spatial resolution ~ 100-micron
nx_i = round(delta_i/(100e-6));   % Initial number of grid cells
dx_i = delta_i/nx_i;              % Initial grid cell size
fprintf(' dx_i = %g \n',dx_i);
fprintf(' nx_i = %g \n',nx_i);
fprintf(' \n');

% calculate initial time step
dt_i=timeStep(FMC, A_R1, Ta_R1, A_R2, Ta_R2, rho_m, rho_vs, rho_c, ...
       k_m, k_vs, k_c, c_m, c_vs, c_c, dx_i );
n_f   = round(T_end/dt_i);   % Total number of time steps
fprintf(' n_f = %g \n',n_f);
fprintf(' \n');
pause;

% Parameter that controls the frequency at which output quantities are saved
i_output = 50;

% Set Initial conditions (t = 0)
n         = 0;
time      = 0; 
n_output  = 0;

nx_new    = nx_i;

temp_surf   = temp_surf_i;        % Surface temperature [K] (t = 0)

sigma       = 5.67e-8;   % Stefan-Boltzmann constant [W/m2/K4]

% Convection coefficient calculation
Pr          = 0.7;      % Air Prandtl number
nu_0        = 1.6d-5;   % Air kinematic viscosity at ambient temp. [m2/s]
k_0         = 0.026;    % Air thermal conductivity at ambient temp. [W/m/K]
T_film      = 0.5*(temp_surf_i+T_g);   % Film temperature [K]
nu_g        = nu_0*(T_film/300)^(1.76); % Kin. visc. at film temp.
k_g         = k_0*(T_film/300)^(0.76); % Conductivity at film temp.
if geometry=="rectangle"
   D_eff       = (2*delta_i)*2/sqrt(pi); % Effective diameter of particle [m]
   Re_D        = u_g*D_eff/nu_g;
   Nu_D        = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
elseif geometry=="cylinder"
   D_eff       = 2*delta_i;
   Re_D        = u_g*D_eff/nu_g;
   Nu_D        = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
            *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
elseif geometry=="sphere"
   D_eff       = 2*delta_i;
   k_g         = k_0*(T_g/300)^(0.76); % Conductivity at gas temp.
   nu_g        = nu_0*(T_g/300)^(1.76); % Kin. visc. at gas temp.
   nu_s        = nu_0*(temp_surf_i/300)^(1.76); % Kin. visc. at surface temp.
   rho_g       = 101325/(287*T_g);
   rho_s       = 101325/(287*temp_surf_i);
   mu_g        = rho_g*nu_g;
   mu_s        = rho_s*nu_s;
   Re_D        = u_g*D_eff/nu_g;       % Calculated at gas temperature
   Nu_D        = 2+(0.4*Re_D^0.5+0.06*Re_D^(2/3))*Pr^0.4*(mu_g/mu_s)^0.25;
else
   fprintf(' Invalid geometry, valid geometries are: cylinder/sphere/rectangle ');
   return;
end
h_conv      = Nu_D*k_g/D_eff;   % Convective heat transfer coef. [W/m2/K]

q_surf_i    = eps_vs*G - eps_vs*sigma*temp_surf^4 ...
              + h_conv*(T_g-temp_surf);  % Net surface heat flux [W/m2] (t = 0)
q_surf      = q_surf_i;

x_vs_i = 1/(1+(FMC*rho_vs/rho_m));  % Vol. frac. of virgin solid (t = 0)
x_m_i  = 1-x_vs_i;                  % Vol. frac. of moisture (t = 0)
x_m_i  = max(0,min(1,x_m_i));       % Safety: enforce 0 <= x_m <= 1
x_c_i  = 0;                         % Vol. frac. of solid char (t = 0)

rho_p_i = rho_m*x_m_i + rho_vs*x_vs_i;   % Solid mass density (t = 0)

dx    = dx_i   *ones(nx_i,1); % Array containing grid cell sizes
dt    = dt_i;                 % Time increment
temp  = temp_surf_i*ones(nx_i,1); % Array containing solid temperatures
rho_p = rho_p_i*ones(nx_i,1); % Array containing solid mass densities
x_m   = x_m_i  *ones(nx_i,1); % Array containing vol. frac. of moisture
x_vs  = x_vs_i *ones(nx_i,1); % Array containing vol. frac. of virgin solid
x_c   = x_c_i  *ones(nx_i,1); % Array containing vol. frac. of solid char

% Allocate additional arrays for saved output quantities
time_save      =          zeros(round(n_f/i_output),1);   % Time
MLR_save       =          zeros(round(n_f/i_output),1);   % Mass loss rate
q_surf_save    = q_surf_i *ones(round(n_f/i_output),1);   % Surf. heat flux
temp_surf_save = temp_surf*ones(round(n_f/i_output),1);   % Surf. temp.
temp_back_save = temp_surf*ones(round(n_f/i_output),1);   % Back temp.
delta_save     = delta_i  *ones(round(n_f/i_output),1);   % Thickness
dt_save        = dt_i     *ones(round(n_f/i_output),1);   % Time increment

nx_save    = nx_i  *ones(round(n_f/i_output),1);    % Number of cells
XC_save    =        ones(nx_i,round(n_f/i_output)); % Coord. cell centers
TEMP_save  = temp_surf*ones(nx_i,round(n_f/i_output)); % Solid temperature
RHOP_save  = rho_p_i*ones(nx_i,round(n_f/i_output)); % Solid mass density
XM_save    = x_m_i *ones(nx_i,round(n_f/i_output)); % Vol. frac. moisture
XVS_save   = x_vs_i*ones(nx_i,round(n_f/i_output)); % Vol. frac. virgin s.
XChar_save = x_c_i*ones(nx_i,round(n_f/i_output));  % Vol. frac. char
RR1_save   =       zeros(nx_i,round(n_f/i_output)); % Reaction rate R1
RR2_save   =       zeros(nx_i,round(n_f/i_output)); % Reaction rate R2
RR3_save   =       zeros(nx_i,round(n_f/i_output)); % Reaction rate R3

% Set the computational grid
[xRight, xCenter, xLeft1, dV] = mesh(nx_i,dx,geometry,A_rectangle,L_cylinder);
volume_i = sum(dV);

% Save output quantities at initial time
n_output            = n_output+1;
time_save(n_output) = time;
fprintf(' n_output = %g, time = %g \n',n_output,time_save(n_output));

% Update array containing coordinates of cell centers
XC_save(1:nx_save(n_output),n_output) = xCenter(1:nx_save(n_output));

flag_burnout = 0;
% Calculate final value of half-thickness of particle (at burn-out)
if(eta_c == 0)      % Non-charring material
    delta_f = 0;
elseif(A_R3 ~= 0)   % Charring material with char oxidation
    delta_f = 0;
else                % Charring material with no char oxidation
    delta_f = delta_i*(eta_c*rho_vs/rho_c)/(1+(FMC*rho_vs/rho_m));
end


%% Time loop
while ((n ~= n_f) && (flag_burnout ~= 2) && (time<T_end))
    
    n      = n + 1;
    
    % Solution at previous time step (t = time-dt)
    nx_old    = nx_new;
    dV_old    = dV;
    dx_old    = dx;
    dt_old    = dt;    
    temp_old  = temp;
    rho_p_old = rho_p;
    x_m_old   = x_m;
    x_vs_old  = x_vs;
    x_c_old   = x_c;
    
    q_surf_old    = q_surf;
    temp_surf_old = temp_surf;
    
    % Solution at new time step (t = time)
    
    % Time step restriction:
    % - Update dt such that abs(Tp(n+1)-Tp(n)) <= Threshold
    %   Estimate Tp(n+1) from explicit treatment    
    [temp] = energy_conservation_explicit(dt_old,q_surf_old, ...
                 temp_old,x_m_old,x_vs_old,x_c_old,...
                 xRight,xCenter,dV,nx_old,geometry,A_rectangle,L_cylinder);
    temp   = temp';
     
    DeltaTemp_max = max( abs(temp-temp_old) );
    
    Threshold = 10;
    if( DeltaTemp_max > 0 )
        dt = dt_old*Threshold/DeltaTemp_max;
        % Avoid strong variations in dt
        if( n > 1)
            dt = max(0.9*dt_old,min(1.1*dt_old,dt)); % Factor 10%
        end
    else
        dt = dt_old;
    end
    
    time   = time + dt;
    if(mod(n,i_output)==0)
        n_output            = n_output+1;
        time_save(n_output) = time;
        fprintf(' \n');
        fprintf(' n_output = %g, time = %g \n', ...
                  n_output,time_save(n_output));
    end    
    
    % - Calculate temperature
    %
    % Note: solve dTp/dt = RHS = (RHS1+RHS2) where RHS1 is energy release
    % due to pyrolysis/char oxidation and RHS2 is heat conduction plus
    % gas-to-solid heat transfer at the exposed surface
    %
    % Discretization (implicit):
    %   (Tp(n+1)-Tp(n))/dt = RHS1 where RHS1 is evaluated at Tp(n+1) using
    %                             linearization
    %   (Tp(n+1)-Tp(n))/dt = RHS2 where RHS2 is evaluated using
    %                             Crank-Nicolson
    %
    % Linearization of RHS1 (pyrolysis/char oxidation):
    %   RHS1 = RHS1_old + dRHS1dT_old*(Tp-temp_old)
    %
    % Operator splitting method (Strang's method)
    %   Step 1: 0.25*dTp/dt = (0.5*RHS1) + (0.5*dRHS1dT)*(Tp-temp_old)
    %           t(n) <= t <= t(n+0.25)
    %   Step 2: 0.5*dTp/dt = (1.0*RHS2)
    %           t(n+0.25) <= t <= t(n+0.75)
    %   Step 3: 0.25*dTp/dt = (0.5*RHS1) + (0.5*dRHS1dT)*(Tp-temp_old)
    %           t(n+0.75) <= t <= t(n+1)
    %    
    
    %   RHS1 = RHS1_old + dRHS1dT_old*(Tp-temp_old)
    %
    % Operator splitting method (Strang's method)
    %   Step 1: 0.25*dTp/dt = (0.5*RHS1) + (0.5*dRHS1dT)*(Tp-temp_old)
    %           t(n) <= t <= t(n+0.25)
    %   Step 2: 0.5*dTp/dt = (1.0*RHS2)
    %           t(n+0.25) <= t <= t(n+0.75)
    %   Step 3: 0.25*dTp/dt = (0.5*RHS1) + (0.5*dRHS1dT)*(Tp-temp_old)
    %           t(n+0.75) <= t <= t(n+1)
    %

    % Step 1: calculate temp1 as the temperature at sub-step t(n+0.25)
    temp0   = temp_old;
    [temp1] = energy_conservation_reactionstep(dt,temp0, ...
                                 temp_old,x_m_old,x_vs_old,x_c_old,nx_old);
    temp1   = temp1';    
    
    % Step 2: calculate temp2 as the temperature at sub-step t(n+0.75)
    temp0        = temp1;
    [a, b, c, d] = energy_conservation_diffusionstep(dt,q_surf_old, ...
                 temp_surf_old,temp0,x_m_old,x_vs_old,x_c_old, h_conv, ...
                 xRight,xCenter,dV,nx_old,geometry,A_rectangle,L_cylinder);                                                         
    temp2        = tri(a,b,c,d);
    temp2        = temp2';    
    
    % Step 3: calculate temp3 as the temperature at sub-step t(n+1)
    temp0   = temp2;
    [temp3] = energy_conservation_reactionstep(dt,temp0, ...
                                 temp_old,x_m_old,x_vs_old,x_c_old,nx_old);
    temp    = temp3';
    
    DeltaTemp_max = max( abs(temp-temp_old) );
    if(mod(n,i_output)==0)
        fprintf(' max(|temp-temp_old|), dt = %g %g \n',DeltaTemp_max,dt);
    end    
    
    % - Calculate species volume fractions and cell volumes
    [x_m, x_vs, x_c, dV] = ...
     mass_conservation(dt,temp_old,x_m_old,x_vs_old,x_c_old,dV_old,nx_old);
    volume = sum(dV);
 
    % Update computational grid (in case of volume change)
    [xRight, xCenter, dx] = moveMesh(nx_old,dV,geometry,A_rectangle,L_cylinder);
    
    % Check
    if(mod(n,i_output)==0)
        sample_thickness = 0;
        for i=1:nx_old
            sample_thickness = sample_thickness + dx_old(i);
        end
        fprintf(' delta, sample_thickness(radius) = %g %g \n', ...
                                          xRight(nx_old),sample_thickness);
    end

    
%Start remeshing    
    % Remeshing, step 1: update nx, dx
    dx_old = dx;
    dV_old = dV;
    [nx_new, dx_new] = remeshing(nx_old,xRight(nx_old));
    if(mod(n,i_output)==0)
        fprintf(' dx = %g, nx = %g \n',dx_new(1),nx_new);
    end
    
    % Remeshing, step 2: update computational grid
    xCenter_old = xCenter;
%    [xRight, xCenter, xLeft1] = mesh(nx_new,dx_new);
    [xRight, xCenter, xLeft1, dV] = mesh(nx_new,dx_new,geometry,A_rectangle,L_cylinder);
    dx = dx_new;
    volume = sum(dV);
    
    % Check
    if(mod(n,i_output)==0)
        sample_thickness = 0;
        for i=1:nx_new
            sample_thickness = sample_thickness + dx_new(i);
        end
        fprintf(' delta_new, sample_thickness = %g %g \n', ...
                                          xRight(nx_new),sample_thickness);
    end
    
    % Remeshing, step 3: interpolate solution on new computational grid
    temp_old = temp;
    x_m_old  = x_m;
    x_vs_old = x_vs;
    x_c_old  = x_c;
        
    [temp, x_m, x_vs, x_c] = interpolateOnNewMesh(nx_new,xCenter, ...
     volume,dV,nx_old,xCenter_old,dV_old,temp_old,x_m_old,x_vs_old,x_c_old);
    temp = temp';
    x_m  = x_m';
    x_vs = x_vs';
    x_c  = x_c';
    %fprintf(' size(temp) = %g %g \n',size(temp));
    %return;
%End remeshing
    
     for i=1:nx_new
         rho_p(i) = rho_m*x_m(i) + rho_vs*x_vs(i) + rho_c*x_c(i);
     end
    
     % Calculate the net surface heat flux and the surface temperature
     kp_surf  =   k_m*x_m(nx_new) +   k_vs*x_vs(nx_new)+   k_c*x_c(nx_new);
     eps_surf = eps_m*x_m(nx_new) + eps_vs*x_vs(nx_new)+ eps_c*x_c(nx_new);
     dx_surf  = xRight(nx_new)-xCenter(nx_new);
     tempc_nx = temp(nx_new);   
     T_film = 0.5*(temp_surf_old+T_g);   % Film temperature [K]
     nu_g   = nu_0*(T_film/300)^(1.76); % Kin. visc. at film temp.
     k_g    = k_0*(T_film/300)^(0.76); % Conductivity at film temp.    
     if geometry=="rectangle"
            D_eff  = (2*xRight(nx_new))*2/sqrt(pi); 
            Re_D   = u_g*D_eff/nu_g;
            Nu_D   = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
                *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
     elseif geometry=="cylinder"
            D_eff  = 2*xRight(nx_new); 
            Re_D   = u_g*D_eff/nu_g;
            Nu_D   = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
                *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
     elseif geometry=="sphere"
            D_eff  = 2*xRight(nx_new); 
            k_g    = k_0*(T_g/300)^(0.76); % Conductivity at gas temp.
            nu_g   = nu_0*(T_g/300)^(1.76); % Kin. visc. at gas temp.
            nu_s   = nu_0*(temp_surf_old/300)^(1.76); % Kin. visc. at surface temp.
            rho_g  = 101325/(287*T_g);
            rho_s  = 101325/(287*temp_surf_old);
            mu_g   = rho_g*nu_g;
            mu_s   = rho_s*nu_s;
            Re_D   = u_g*D_eff/nu_g;       % Calculated at gas temperature
            Nu_D   = 2+(0.4*Re_D^0.5+0.06*Re_D^(2/3))*Pr^0.4*(mu_g/mu_s)^0.25;
     else
           fprintf(' Invalid geometry, valid geometries are: cylinder/sphere/rectangle ');
           return;
     end
     h_conv = Nu_D*k_g/D_eff;   % Convective heat transfer coef. [W/m2/K]    
     [q_surf, temp_surf] = ...
         particle_surface(temp_surf_old,kp_surf,eps_surf,dx_surf,tempc_nx,h_conv);
 
     
     % Save data
     if(mod(n,i_output)==0) 
         % Calculate the mass loss rate 
         %  - Mass loss rate per unit volume [kg/s/m3]
         for i=1:nx_new
            MLRpuv(i) = rho_m *x_m(i) *A_R1*exp(-Ta_R1/temp(i)) ...
                      + rho_vs*x_vs(i)*A_R2*exp(-Ta_R2/temp(i)) ...
                                      *(1-eta_c) ...
                      + rho_c *x_c(i) *x_O2_g ...
                                      *A_R3*exp(-Ta_R3/temp(i));
         end
         %  - Total mass loss rate [kg/s]
         MLRtot_new = 0;
         for i=1:nx_new
             MLRtot_new = MLRtot_new + MLRpuv(i)*dV(i);
         end
         if geometry=="rectangle"
             MLRtot_new = MLRtot_new * 2 ; %because we solve half rectangles
         end
         
         MLR_save(n_output)       = MLRtot_new;
         q_surf_save(n_output)    = q_surf;
         temp_surf_save(n_output) = temp_surf;
         temp_back_save(n_output) = temp(1);
         delta_save(n_output)     = xRight(nx_new);        
         dt_save(n_output)        = dt;
        
         nx_save(n_output)            = nx_new;
         XC_save(1:nx_new,n_output)   = xCenter(1:nx_new);
         TEMP_save(1:nx_new,n_output) = temp(1:nx_new);
         RHOP_save(1:nx_new,n_output) = rho_p(1:nx_new);
         XM_save(1:nx_new,n_output)   = x_m(1:nx_new);
         XVS_save(1:nx_new,n_output)  = x_vs(1:nx_new);
         XChar_save(1:nx_new,n_output)  = x_c(1:nx_new);
         
         for i=1:nx_new
            RR1_save(i,n_output) = rho_m *x_m(i) ...
                                             *A_R1*exp(-Ta_R1/temp(i));
            RR2_save(i,n_output) = rho_vs*x_vs(i) ...
                                             *A_R2*exp(-Ta_R2/temp(i));
            RR3_save(i,n_output) = rho_c *x_c(i) *x_O2_g ...
                                             *A_R3*exp(-Ta_R3/temp(i));
         end
     end
     
     % Check for burnout time (Safety: do not allow shrinking to zero-size)
     if( ( flag_burnout == 0 ) && ...
         ( eta_c ~= 0 ) && ( x_vs(1) < 0.001 ) )
         flag_burnout  = 1;
         n_last_output = n_output+1; % Simulate until time for last output
         fprintf(' \n');
         fprintf(' Burnout time = %g \n',time);
         fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
         fprintf(' Thickness at burnout time = %g \n',delta_f);
     end
     if( ( flag_burnout == 0 ) && ...
        ( eta_c == 0 ) && ( (xRight(nx_new)/delta_i) < 0.02 ) )
         flag_burnout  = 1;
         n_last_output = n_output;   % Stop simulation
         fprintf(' \n');
         fprintf(' Burnout time = %g \n',time);
         fprintf(' Thickness at burnout time = %g \n',xRight(nx_new));
         fprintf(' Thickness at burnout time = %g \n',delta_f);
     end
     if( ( flag_burnout == 1 ) & ( n_output == n_last_output ) )
         flag_burnout = 2;
     end
     
end

%% Output data
csvwrite("num_cells.csv",nx_save(1:n_output));
csvwrite("time.csv",time_save(1:n_output));
csvwrite("cell_center.csv",XC_save(:,1:n_output));
csvwrite("MLR.csv",MLR_save(1:n_output));
csvwrite("T.csv",TEMP_save(:,1:n_output));
csvwrite("delta.csv",delta_save(1:n_output));
csvwrite("x_m.csv",XM_save(:,1:n_output));
csvwrite("x_vs.csv",XVS_save(:,1:n_output));
csvwrite("x_c.csv",XChar_save(:,1:n_output));
csvwrite("T_surf.csv",temp_surf_save(1:n_output));
csvwrite("T_core.csv",temp_back_save(1:n_output)); 
csvwrite("q_surf.csv", q_surf_save(1:n_output));

% Check global conservation statements
fprintf(' \n');
fprintf(' Final thickness of the particle (theory)     = %g \n', ...
                                                                  delta_f);
fprintf(' Final thickness of the particle (simulation) = %g \n', ...
                                                           xRight(nx_new));

% Initial mass of moisture and virgin solid (t = 0) [kg]
if geometry=="rectangle"
Mm_i  = FMC*rho_vs*2*delta_i*A_rectangle/(1+(FMC*rho_vs/rho_m)); % Mass of m
Mvs_i =     rho_vs*2*delta_i*A_rectangle/(1+(FMC*rho_vs/rho_m)); % Mass of vs
elseif geometry=="cylinder"
Mm_i  = FMC*rho_vs*pi*delta_i^2*L_cylinder/(1+(FMC*rho_vs/rho_m)); % Mass of m
Mvs_i =     rho_vs*pi*delta_i^2*L_cylinder/(1+(FMC*rho_vs/rho_m)); % Mass of vs
elseif geometry=="sphere"
Mm_i  = FMC*rho_vs*4/3*pi*delta_i^3/(1+(FMC*rho_vs/rho_m)); % Mass of m
Mvs_i =     rho_vs*4/3*pi*delta_i^3/(1+(FMC*rho_vs/rho_m)); % Mass of vs
end

% Total mass of water vapor and volatiles released (t = inf) [kg]
if(eta_c == 0)      % Non-charring material
    Mg_f = Mm_i + Mvs_i;
elseif(A_R3 ~= 0)   % Charring material with char oxidation
    Mg_f = Mm_i + Mvs_i;
else                % Charring material with no char oxidation
    Mg_f = Mm_i + (1-eta_c)*Mvs_i;
end

int_MLR = 0;
for n=2:n_output
    int_MLR = int_MLR + 0.5*(MLR_save(n-1)+MLR_save(n)) ...
                                 *(time_save(n)-time_save(n-1));
end  
    
fprintf(' \n');
fprintf(' Total mass of vapor/volatiles released (theory)     = %g \n', ...
                                                                  Mg_f);
fprintf(' Total mass of vapor/volatiles released (simulation) = %g \n', ...
                                                               int_MLR);
% Plot results
% - Time variations of key quantities

figure(1);   % Mass loss rate [kg/s/m2]
hold on;
plot(time_save(1:n_output),(MLR_save(1:n_output)));


figure(2);   % Net surface heat flux [kW/m2]
hold on;
plot(time_save(1:n_output),(q_surf_save(1:n_output)/1000));

figure(3);   % Temperature at exposed surface and at center of particle [K]
hold on;
plot(time_save(1:n_output),temp_surf_save(1:n_output));
hold on;
plot(time_save(1:n_output),temp_back_save(1:n_output));

figure(4);   % Half-thickness of particle
hold on;
plot(time_save(1:n_output),delta_save(1:n_output));

figure(5);   % Time increment [s]
hold on;
plot(time_save(1:n_output),dt_save(1:n_output));

% - Spatial profiles

figure(11);   % Temperature
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),TEMP_save(1:nx_save(n),n));
    hold on
end

figure(12);   % Mass density
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),RHOP_save(1:nx_save(n),n));
    hold on
end

figure(13);   % Volume fraction of moisture
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),XM_save(1:nx_save(n),n));
    hold on
end

figure(14);   % Volume fraction of virgin solid
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),XVS_save(1:nx_save(n),n));
    hold on
end

figure(15);   % Volume fraction of char
for n=1:n_output
    plot(XC_save(1:nx_save(n),n), ...
                     (1-XM_save(1:nx_save(n),n)-XVS_save(1:nx_save(n),n)));
    hold on
end

figure(16);   % Mass reaction rate for moisture evaporation reaction R1
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),RR1_save(1:nx_save(n),n));
    hold on
end

figure(17);   % Mass reaction rate for pyrolysis reaction R2
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),RR2_save(1:nx_save(n),n));
    hold on
end

figure(18);   % Mass reaction rate for char oxidation reaction R3
for n=1:n_output
    plot(XC_save(1:nx_save(n),n),RR3_save(1:nx_save(n),n));
    hold on
end

return; 