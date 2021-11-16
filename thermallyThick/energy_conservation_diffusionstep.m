% Energy conservation
% - Sub-step associated with energy change due to diffusive/convective
%   transport of heat and to gas-to-solid heat transfer at the exposed
%   surface
% - Implicit time integration
%
%   (Tp(n+r1)-Tp(n+r0))/dt = RHS2 where RHS2 is evaluated using
%                            Crank-Nicolson
%
%   Evaluation of RHS2:
%      RHS2 = 0.5*RHS2(n+r1) + 0.5*RHS2(n+r0)
%      where r0 and r1 designate fractional steps, 0 <= r0,r1 <= 1
%
% - Formulation as a tri-diagonal matrix system
%
% - Definitions
%      tempr0_newiter = Tp(n+r0,newiter)

function [a, b, c, d] = ...
              energy_conservation_diffusionstep(dt, tempr0_newiter, ...
              temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
              q_surf_old, tempc_nx, temp_surf_old, h_conv, ...
              xRight, xCenter, dV_old, nx_old)
   
global geometry A_rectangle L_cylinder ...
       T_g G pres_g k_g0 cp_g0 nu_g0 MW_g R sigma ...  
       rho_ws rho_ds rho_c rho_a ...
       k_ws k_ds k_c k_a ...
       c_ws c_ds c_c c_a ...
       eps_ws eps_ds eps_c eps_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a
   
% Calculation of surface areas of cell boundary faces
if geometry=="rectangle"    
    S_pos(1)=A_rectangle;
    S_neg(1)=A_rectangle;
    for i=2:nx_old
    S_pos(i)=A_rectangle;
    S_neg(i)=A_rectangle;
    end
elseif geometry=="cylinder" 
    S_pos(1)=2*pi*xRight(1)*L_cylinder;
    S_neg(1)=0;
    for i=2:nx_old
    S_pos(i)=2*pi*xRight(i)*L_cylinder;
    S_neg(i)=2*pi*xRight(i-1)*L_cylinder;
    end
elseif geometry=="sphere"
    S_pos(1)=4*pi*xRight(1)^2;
    S_neg(1)=0;
    for i=2:nx_old
    S_pos(i)=4*pi*xRight(i)^2;
    S_neg(i)=4*pi*xRight(i-1)^2; 
    end   
end
   
for i = 1:nx_old
    % Porosity
    psi     = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
            + psi_c * x_c_old(i) + psi_a * x_a_old(i);
    
    % Effective thermal conductivity
    %  - Porous medium treatment: 
    %    keff = (1-psi) x k_s + psi x k_g
    k_s     = k_ws*x_ws_old(i) + k_ds*x_ds_old(i) ...
            + k_c * x_c_old(i) + k_a * x_a_old(i); % Conduct. in solid p.
    k_g     = k_g0*(temp_old(i)/300)^(0.76);    % Conductivity in gas phase
    keff(i) = (1-psi)*k_s + psi*k_g; % Estimated at time t(n)
    
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s      = rho_ws*c_ws*x_ws_old(i) ...
                 + rho_ds*c_ds*x_ds_old(i) ...
                 + rho_c *c_c * x_c_old(i) ...
                 + rho_a *c_a * x_a_old(i); % (rho*cp) in solid phase 
    rhocp_g      = (pres_g*MW_g/R/temp_old(i))*cp_g0; % (rho*cp) in gas p.
    rhocp_eff(i) = (1-psi)*rhocp_s + psi*rhocp_g; % Estimated at time t(n)
    
    % Permeability
    Kperm(i) = Kperm_ws*x_ws_old(i) + Kperm_ds*x_ds_old(i) ...
             + Kperm_c * x_c_old(i) + Kperm_a * x_a_old(i); % At time t(n)
end
               
% i = 1 (back surface)
i    = 1;
FO_L = 0;
FO_R = 0.5*(keff(i)+keff(i+1))/rhocp_eff(i) * S_pos(i)*dt / ...
                              (dV_old(i)*(xCenter(i+1)-xCenter(i)));

nu_g = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
CFL_L = 0;
if( pres_old(i) >= pres_old(i+1) )
    CFL_R = 0;
else
    CFL_R = dt*(Kperm(i)*cp_g0/nu_g/rhocp_eff(i)) ...
            *(pres_old(i)-pres_old(i+1))/(xCenter(i+1)-xCenter(i))^2;
    % Note: CFL_R < 0
end

a(i) = 0;
b(i) = 1 + 0.5*FO_R - 0.5*CFL_R;
c(i) =-0.5*FO_R + 0.5*CFL_R;
d(i) = (1 - 0.5*FO_R + 0.5*CFL_R)*tempr0_newiter(i) ...
     + (0.5*FO_R - 0.5*CFL_R)*tempr0_newiter(i+1);

% i = nx_old (exposed surface)
i    = nx_old;
FO_L = 0.5*(keff(i)+keff(i-1))/rhocp_eff(i) * S_neg(i)*dt / ...
                              (dV_old(i)*(xCenter(i)-xCenter(i-1)));
FO_R = 0;

nu_g  = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
if( pres_old(i-1) > pres_old(i) )
    CFL_L = dt*(Kperm(i)*cp_g0/nu_g/rhocp_eff(i)) ...
            *(pres_old(i-1)-pres_old(i))/(xCenter(i)-xCenter(i-1))^2;
    % Note: CFL_L > 0
else
    CFL_L = 0;
end
CFL_R = 0;

a(i) =-0.5*FO_L -0.5*CFL_L;
b(i) = 1 + 0.5*FO_L + 0.5*CFL_L;
c(i) = 0;
d(i) = (0.5*FO_L + 0.5*CFL_L)*tempr0_newiter(i-1) ...
     + (1 - 0.5*FO_L - 0.5*CFL_L)*tempr0_newiter(i);

% Contribution of gas-to-sold heat flux
keff_surf = keff(nx_old);
eps_surf  = eps_ws*x_ws_old(nx_old) + eps_ds*x_ds_old(nx_old) ...
          + eps_c * x_c_old(nx_old) + eps_a * x_a_old(nx_old);
dx_surf   = xRight(nx_old)-xCenter(nx_old);
h_rad     = eps_surf*sigma*(temp_surf_old^3);
Bi        = (h_conv+h_rad)*dx_surf/keff_surf;

temp_surf_newiter = ...
           ( tempc_nx + (h_conv*T_g+eps_surf*G)*dx_surf/keff_surf )/(1+Bi);

b(i)      = b(i) + ((h_conv+h_rad)/(1+Bi)) ...
                                 *0.5*dt/rhocp_eff(i) * S_pos(i)/dV_old(i);
d(i)      = d(i) + ( q_surf_old ...
                        + ((h_conv*T_g+eps_surf*G)/(1+Bi)) ) ...
                                 *0.5*dt/rhocp_eff(i) * S_pos(i)/dV_old(i);

% 2 <= i <= (nx_old-1) (interior nodes)
for i = 2:(nx_old-1)   
    FO_L = 0.5*(keff(i)+keff(i-1))/rhocp_eff(i) * S_neg(i)*dt / ...
                                  (dV_old(i)*(xCenter(i)-xCenter(i-1)));
    FO_R = 0.5*(keff(i)+keff(i+1))/rhocp_eff(i) * S_pos(i)*dt / ...
                                  (dV_old(i)*(xCenter(i+1)-xCenter(i)));
    
    nu_g = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
    if( (pres_old(i-1) > pres_old(i)) & (pres_old(i) > pres_old(i+1)) )
        % Case for which mdotPUA > 0. Note: CFL_L > 0
        CFL_L = dt*(Kperm(i)*cp_g0/nu_g/rhocp_eff(i)) ...
                *(pres_old(i-1)-pres_old(i))/(xCenter(i)-xCenter(i-1))^2;
        CFL_R = 0;
    elseif ( (pres_old(i-1) < pres_old(i)) ...
                                          & (pres_old(i) < pres_old(i+1)) )
        % Case for which mdotPUA < 0. Note: CFL_R < 0
        CFL_L = 0;
        CFL_R = dt*(Kperm(i)*cp_g0/nu_g/rhocp_eff(i)) ...
                *(pres_old(i)-pres_old(i+1))/(xCenter(i+1)-xCenter(i))^2;
    else
        % Case for which dp/dx ~ 0
        CFL_L = 0;
        CFL_R = 0;
    end
        
    a(i) =-0.5*FO_L -0.5*CFL_L;
    b(i) = 1 + 0.5*FO_L + 0.5*FO_R + 0.5*CFL_L - 0.5*CFL_R;
    c(i) =-0.5*FO_R + 0.5*CFL_R;
    d(i) = (0.5*FO_L + 0.5*CFL_L)*tempr0_newiter(i-1) ...
         + (1 - 0.5*FO_L - 0.5*FO_R - 0.5*CFL_L + 0.5*CFL_R) ...
           *tempr0_newiter(i) ...
         + (0.5*FO_R - 0.5*CFL_R)*tempr0_newiter(i+1);
end
            
end
