% Mass conservation
% - Semi-implicit time integration

function [x_m, x_vs, x_c, dV] = ...
         mass_conservation_TTC(dt,temp_old,x_m_old,x_vs_old,x_c_old,dV_old)

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c x_O2_g

% Species mass conservation statements
% - Notations:  xm_times_dV =  x_m x dV (at new time step, t = time)
%              xvs_times_dV = x_vs x dV (at new time step, t = time)
%               xc_times_dV =  x_c x dV (at new time step, t = time)

xm_times_dV  = (x_m_old*dV_old)*exp( -dt*A_R1*exp(-Ta_R1/temp_old) );

xvs_times_dV = (x_vs_old*dV_old)/(1+dt*A_R2*exp(-Ta_R2/temp_old));

if( eta_c ~=0 )
    store       = xvs_times_dV ...
                  *dt*A_R2*exp(-Ta_R2/temp_old)*(eta_c*rho_vs/rho_c);
    xc_times_dV = (x_c_old*dV_old + store) ...
                  /(1+dt*x_O2_g*A_R3*exp(-Ta_R3/temp_old));
else
    xc_times_dV = 0;
end

dV = xm_times_dV + xvs_times_dV + xc_times_dV;

x_m  = xm_times_dV /dV;
x_m  = max(0,min(1,x_m));    % Safety: enforce 0 <= x_m  <= 1
x_vs = xvs_times_dV/dV;
x_vs = max(0,min(1,x_vs));   % Safety: enforce 0 <= x_vs <= 1
x_c  = 1-x_m-x_vs;
x_c  = max(0,min(1,x_c));    % Safety: enforce 0 <= x_c  <= 1


