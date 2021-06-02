% Energy conservation
% - Sub-step associated with energy release due to pyrolysis/char oxidation

function [temp] = energy_conservation_reactionstep(dt,temp0, ...
                                  temp_old,x_m_old,x_vs_old,x_c_old,nx_old)

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g

for i=1:nx_old
    
    % Product of mass density times specific heat
    % - Notations: rho_times_cp = rho_p x cp
    rho_times_cp = rho_m*c_m*x_m_old(i) + rho_vs*c_vs*x_vs_old(i) ...
                                           + rho_c*c_c*x_c_old(i);

    % Volumetric rate of heat production/consumption [W/m3]                    
    Qdotp = rho_m * x_m_old(i)*A_R1*exp(-Ta_R1/temp_old(i))*DeltaH_R1...
             + rho_vs*x_vs_old(i)*A_R2*exp(-Ta_R2/temp_old(i))*DeltaH_R2...
               *(1-eta_c) ...
             + rho_c * x_c_old(i)*x_O2_g ...
                                 *A_R3*exp(-Ta_R3/temp_old(i))*DeltaH_R3;
                             
    dQdotpdT = rho_m * x_m_old(i)*A_R1*exp(-Ta_R1/temp_old(i)) ...
                     *DeltaH_R1*(Ta_R1/(temp_old(i)^2)) ...
             + rho_vs*x_vs_old(i)*A_R2*exp(-Ta_R2/temp_old(i)) ...
                     *DeltaH_R2*(1-eta_c)*(Ta_R2/(temp_old(i)^2)) ...
             + rho_c * x_c_old(i)*x_O2_g*A_R3*exp(-Ta_R3/temp_old(i)) ...
                     *DeltaH_R3*(Ta_R3/(temp_old(i)^2));    
                 
    RHS1_old    =    Qdotp/rho_times_cp;
    dRHS1dT_old = dQdotpdT/rho_times_cp;
    
    if( abs(dRHS1dT_old) > 0 )
        ratio   = RHS1_old/dRHS1dT_old;
        temp(i) = temp0(i) ...
                + (ratio+temp0(i)-temp_old(i))*(exp(0.5*dRHS1dT_old*dt)-1);
    else % Case dRHS1dT_old = 0
        temp(i) = temp0(i) + 0.5*RHS1_old*dt;
    end
    
end
            
end