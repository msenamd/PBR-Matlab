% Particle Burning Rate Model - Thermally Thin Case - PBR_model_TTC.m

clc;
clear;
close all;

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       k_m k_vs k_c c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g ...
       T_g G eps_m eps_vs eps_c

% Read input parameters
[T_g, u_g, G, x_O2_g, delta_i, FMC, ...
 rho_m, rho_vs, rho_c, k_m, k_vs, k_c, c_m, c_vs, c_c, ...
 eps_m, eps_vs, eps_c, A_R1, Ta_R1, DeltaH_R1, ...
 A_R2, Ta_R2, DeltaH_R2, eta_c, A_R3, Ta_R3, DeltaH_R3,...
 temp_surf_i,T_end,geometry,A_rectangle,L_cylinder] = input_parameters;

sigma = 5.67e-8;   % Stefan-Boltzmann constant [W/m2/K4]

%
% Begin analysis of chemical time scales and select time step
tau_RR1_min = 1000;
for i=1:22
    tem_tmp = 300 + (i-1)*200/20; % Cover range 300 <= temp_tmp <= 500 K
    if(A_R1 ~= 0)
        tau_RR1 = 1/(A_R1*exp(-Ta_R1/tem_tmp));
        tau_RR1_min = min(tau_RR1_min,tau_RR1);
    end
end
fprintf(' \n');
fprintf(' tau_RR1_min = %g \n',tau_RR1_min);
tau_RR2_min = 1000;
for i=1:101
    tem_tmp = 300 + (i-1)*700/100; % Cover range 300 <= temp_tmp <= 1000 K
    if(A_R2 ~= 0)
        tau_RR2 = 1/(A_R2*exp(-Ta_R2/tem_tmp));
        tau_RR2_min = min(tau_RR2_min,tau_RR2);
    end
end
fprintf(' tau_RR2_min = %g \n',tau_RR2_min);

% Select time step (constant)
% - Enforce (tau_RR1_min/dt) and (tau_RR2_min/dt) >=10
dt = min( (tau_RR1_min/1), (tau_RR2_min/1) );
fprintf(' \n');
fprintf(' dt = %g \n',dt);
fprintf(' (tau_RR1_min/dt),(tau_RR2_min/dt) = %g %g \n', ...
                                        (tau_RR1_min/dt),(tau_RR2_min/dt));
% End analysis of chemical time scales
%

n_f   = round(T_end/dt);   % Total number of time steps
fprintf(' n_f = %g \n',n_f);
fprintf(' \n');

% Parameter that contols the frequency at which output quantities are saved
i_output = 10;

% Initial conditions (t = 0)
n         = 0;
time      = 0; 
n_output  = 0;

x_vs_i = 1/(1+(FMC*rho_vs/rho_m));  % Vol. frac. of virgin solid (t = 0)
x_m_i  = 1-x_vs_i;                  % Vol. frac. of moisture (t = 0)
x_m_i  = max(0,min(1,x_m_i));       % Safety: enforce 0 <= x_m <= 1
x_c_i  = 0;                         % Vol. frac. of solid char (t = 0)

rho_p_i = rho_m*x_m_i + rho_vs*x_vs_i;   % Solid mass density (t = 0)

dx    = delta_i;   % Thickness or radius of particle
% Volume of particle
if geometry=="rectangle"    
   dV=dx*A_rectangle;
elseif geometry=="cylinder"
    dV=pi*dx^2*L_cylinder;
elseif geometry=="sphere"
    dV=4/3*pi*dx^3;
else
    fprintf('Incorrect geometry type \n');
    fprintf('Available types: rectangle / sphere / cylinder \n');
    return;
end

temp  = temp_surf_i; % Solid temperature
rho_p = rho_p_i;     % Solid mass density
x_m   = x_m_i;       % Vol. frac. of moisture
x_vs  = x_vs_i;      % Vol. frac. of virgin solid
x_c   = x_c_i;       % Vol. frac. of solid char

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

q_surf_i    = eps_vs*G - eps_vs*sigma*temp^4 ...
              + h_conv*(T_g-temp);  % Net surface heat flux [W/m2] (t = 0)
q_surf      = q_surf_i;


% Allocate additional arrays for saved output quantities
time_save   =          zeros(round(n_f/i_output),1);   % Time
MLR_save    =          zeros(round(n_f/i_output),1);   % Mass loss rate
q_surf_save = q_surf_i *ones(round(n_f/i_output),1);   % Surf. heat flux
delta_save  = delta_i  *ones(round(n_f/i_output),1);   % Thickness

temp_save = temp_surf_i  *ones(round(n_f/i_output),1); % Solid temperature
rhop_save = rho_p_i*ones(round(n_f/i_output),1); % Solid mass density
xm_save   = x_m_i  *ones(round(n_f/i_output),1); % Vol. frac. moisture
xvs_save  = x_vs_i *ones(round(n_f/i_output),1); % Vol. frac. virgin solid
RR1_save  =        zeros(round(n_f/i_output),1); % Reaction rate R1
RR2_save  =        zeros(round(n_f/i_output),1); % Reaction rate R2
RR3_save  =        zeros(round(n_f/i_output),1); % Reaction rate R3

% Save output quantities at initial time
n_output            = n_output+1;
time_save(n_output) = time;
fprintf(' n_output = %g, time = %g \n',n_output,time_save(n_output));

flag_burnout = 0;
% Calculate final value of half-thickness of particle (at burn-out)
if(eta_c == 0)      % Non-charring material
    delta_f = 0;
elseif(A_R3 ~= 0)   % Charring material with char oxidation
    delta_f = 0;
else                % Charring material with no char oxidation
    delta_f = delta_i*(eta_c*rho_vs/rho_c)/(1+(FMC*rho_vs/rho_m));
end


% Time loop 
while ((n ~= n_f) & (flag_burnout ~= 2))
    
    n      = n + 1;
    time   = time + dt;
    
    if(mod(n,i_output)==0)
        n_output            = n_output+1;
        time_save(n_output) = time;
        fprintf(' \n');
        fprintf(' n_output = %g, time = %g \n', ...
                  n_output,time_save(n_output));
    end
    
    % Solution at previous time step (t = time-dt)
    dV_old    = dV;
    if geometry=="rectangle"
        dx_old=dV_old/A_rectangle;
    elseif geometry=="cylinder"
        dx_old=(dV_old/pi/L_cylinder)^0.5;
    elseif geometry=="sphere"
        dx_old=(dV_old/(4/3*pi))^(1/3);
    end
    temp_old  = temp;
    rho_p_old = rho_p;
    x_m_old   = x_m;
    x_vs_old  = x_vs;
    x_c_old   = x_c;
    
    q_surf_old = q_surf;
    
    % Solution at new time step (t = time)
    % - Calculate species volume fractions and cell volumes
    [x_m, x_vs, x_c, dV] = ...
        mass_conservation_TTC(dt,temp_old,x_m_old,x_vs_old,x_c_old,dV_old);
    
    % - Calculate particle thickness or radius
    if geometry=="rectangle"
        dx=dV/A_rectangle;
    elseif geometry=="cylinder"
        dx=(dV/pi/L_cylinder)^0.5;
    elseif geometry=="sphere"
        dx=(dV/(4/3*pi))^(1/3);
    end
    
    % - Calculate temperature
    [temp] = energy_conservation_TTC(dt,q_surf_old, ...
                            temp_old,x_m_old,x_vs_old,x_c_old,...
                            dx_old,dV_old,geometry,A_rectangle,L_cylinder);
    
    rho_p = rho_m*x_m + rho_vs*x_vs + rho_c*x_c;
    
    % Calculate the net surface heat flux and the surface temperature
    eps_surf = eps_m*x_m + eps_vs*x_vs + eps_c*x_c;
    T_film = 0.5*(temp+T_g);   % Film temperature [K]
    nu_g   = nu_0*(T_film/300)^(1.76); % Kin. visc. at film temp.
    k_g         = k_0*(T_film/300)^(0.76); % Conductivity at film temp.
    if geometry=="rectangle"
            D_eff  = (2*dx)*2/sqrt(pi); 
            Re_D   = u_g*D_eff/nu_g;
            Nu_D   = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
                *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
    elseif geometry=="cylinder"
            D_eff  = 2*dx; 
            Re_D   = u_g*D_eff/nu_g;
            Nu_D   = 0.3 + 0.62*(Re_D^0.5)*(Pr^0.333) ...
                *( (1+(Re_D/282000)^(5/8))^(4/5) )/( (1+(0.4/Pr)^(2/3))^0.25 );
    elseif geometry=="sphere"
            D_eff  = 2*dx; 
            k_g    = k_0*(T_g/300)^(0.76); % Conductivity at gas temp.
            nu_g   = nu_0*(T_g/300)^(1.76); % Kin. visc. at gas temp.
            nu_s   = nu_0*(temp/300)^(1.76); % Kin. visc. at surface temp.
            rho_g  = 101325/(287*T_g);
            rho_s  = 101325/(287*temp);
            mu_g   = rho_g*nu_g;
            mu_s   = rho_s*nu_s;
            Re_D   = u_g*D_eff/nu_g;       % Calculated at gas temperature
            Nu_D   = 2+(0.4*Re_D^0.5+0.06*Re_D^(2/3))*Pr^0.4*(mu_g/mu_s)^0.25;
    else
           fprintf(' Invalid geometry, valid geometries are: cylinder/sphere/rectangle ');
           return;
    end
    h_conv = Nu_D*k_g/D_eff;
    
    q_surf   = eps_surf*G - eps_surf*sigma*(temp^4) ...
                          + h_conv*(T_g-temp);
     
    if(mod(n,i_output)==0)
        % Calculate the mass loss rate
        %  - Mass loss rate per unit volume [kg/s/m3]
        MLRpuv = rho_m *x_m *A_R1*exp(-Ta_R1/temp) ...
               + rho_vs*x_vs*A_R2*exp(-Ta_R2/temp) *(1-eta_c) ...
               + rho_c *x_c *x_O2_g ...
                            *A_R3*exp(-Ta_R3/temp);
                           
        %  - Total mass loss rate per unit exposed surface area [kg/s]
        MLR_new = MLRpuv*dV;
        
         
        MLR_save(n_output)    = MLR_new;
        q_surf_save(n_output) = q_surf;
        delta_save(n_output)  = dx;
        
        temp_save(n_output)   = temp;
        rhop_save(n_output)   = rho_p;
        xm_save(n_output)     = x_m;
        xvs_save(n_output)    = x_vs;        
        RR1_save(n_output)    = rho_m *x_m *A_R1*exp(-Ta_R1/temp);
        RR2_save(n_output)    = rho_vs*x_vs*A_R2*exp(-Ta_R2/temp);
        RR3_save(n_output)    = rho_c *x_c *x_O2_g ...
                                           *A_R3*exp(-Ta_R3/temp);
    end
     
    % Check for burnout time (Safety: do not allow shrinking to zero-size)
    if( ( flag_burnout == 0 ) & ...
        ( eta_c ~= 0 ) & ( x_vs < 0.001 ) )
        flag_burnout  = 1;
        n_last_output = n_output+1; % Simulate until time for last output
        fprintf(' \n');
        fprintf(' Burnout time = %g \n',time);
        fprintf(' Thickness at burnout time = %g \n',dx);
        fprintf(' Thickness at burnout time = %g \n',delta_f);
    end
    if( ( flag_burnout == 0 ) & ...
        ( eta_c == 0 ) & ( (dx/delta_i) < 0.001 ) )
        flag_burnout  = 1;
        n_last_output = n_output;   % Stop simulation
        fprintf(' \n');
        fprintf(' Burnout time = %g \n',time);
        fprintf(' Thickness at burnout time = %g \n',dx);
        fprintf(' Thickness at burnout time = %g \n',delta_f);
    end
    if( ( flag_burnout == 1 ) & ( n_output == n_last_output ) )
        flag_burnout = 2;
    end
     
end

% Output data
csvwrite("n_output.csv",n_output);
csvwrite("time.csv",time_save(1:n_output));
csvwrite("MLR.csv",MLR_save(1:n_output));
csvwrite("T.csv",temp_save(1:n_output));
csvwrite("delta.csv",delta_save(1:n_output));
csvwrite("x_m.csv",xm_save(1:n_output));
csvwrite("x_vs.csv",xvs_save(1:n_output));
csvwrite("x_c.csv",1-xm_save(1:n_output)-xvs_save(1:n_output));
csvwrite("q_surf.csv", q_surf_save(1:n_output));

% Initial mass of moisture and virgin solid (t = 0) [kg]
if geometry=="rectangle"
Mm_i  = FMC*rho_vs*delta_i*A_rectangle/(1+(FMC*rho_vs/rho_m)); % Mass of m
Mvs_i =     rho_vs*delta_i*A_rectangle/(1+(FMC*rho_vs/rho_m)); % Mass of vs
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
return; 