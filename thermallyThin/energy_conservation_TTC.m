% Energy conservation
% - Explicit time integration

function [temp] = energy_conservation_TTC(dt,q_surf_old, ...
                   temp_old,x_m_old,x_vs_old,x_c_old,...
                   dx_old,dV_old,geometry,A_rectangle,L_cylinder)

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       k_m k_vs k_c c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g

    
% Product of mass density times specific heat
% - Notations: rho_times_cp = rho_p x cp
rho_times_cp = rho_m*c_m*x_m_old + rho_vs*c_vs*x_vs_old ...
                                 + rho_c*c_c*x_c_old;

% Volumetric rate of heat production/consumption [W/m3]                    
Qdotp = rho_m * x_m_old*A_R1*exp(-Ta_R1/temp_old)*DeltaH_R1...
      + rho_vs*x_vs_old*A_R2*exp(-Ta_R2/temp_old)*DeltaH_R2*(1-eta_c) ...
      + rho_c * x_c_old*x_O2_g ...
                       *A_R3*exp(-Ta_R3/temp_old)*DeltaH_R3;

% Energy conservation statement
if geometry=="rectangle"
    temp = temp_old + (Qdotp*dt/rho_times_cp) ...
                    + (q_surf_old*A_rectangle*dt/rho_times_cp/dV_old);
elseif geometry=="cylinder"
    temp = temp_old + (Qdotp*dt/rho_times_cp) ...
                    + (q_surf_old*2*pi*dx_old*L_cylinder*dt/rho_times_cp/dV_old);                
elseif geometry=="sphere"
    temp = temp_old + (Qdotp*dt/rho_times_cp) ...
                    + (q_surf_old*4*pi*dx_old^2*dt/rho_times_cp/dV_old); 
end
            
end