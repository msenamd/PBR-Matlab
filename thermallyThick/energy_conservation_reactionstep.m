% Energy conservation
% - Sub-step associated with energy change due to reactions (R1)-(R4)
% - Implicit time integration
%
%   (Tp(n+r1)-Tp(n+r0))/dt = RHS1(n+r1)
%                            where r0 and r1 designate fractional steps,
%                            0 <= r0,r1 <= 1, and where RHS1 is evaluated
%                            at time (n+r1)
%
%   Decomposition of RHS1 into positive and negative terms:
%      RHS1 = RHS1p + RHS1m, where RHS1p > 0 and RHS1m <= 0
%
%      Tp(n+r1) = ( Tp(n+r0) + dt*RHS1p(n+r1) ) 
%                /( 1 - dt*RHS1m(n+r1)/Tp(n+r1) )
% 
% - Definitions
%      tempr0_newiter = Tp(n+r0,newiter)
%      tempr1_olditer = Tp(n+r1,olditer)
%      tempr1_newiter = Tp(n+r1,newiter)

function [tempr1_newiter] = ...
   energy_conservation_reactionstep(dt, tempr0_newiter, tempr1_olditer, ...
    x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, Y_O2_olditer, ...
    lambda, temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, nx_old)


global pres_g cp_g0 MW_g R ...
       rho_ws rho_ds rho_c rho_a ...
       c_ws c_ds c_c c_a ...      
       psi_ws psi_ds psi_c psi_a ...
       A_R1 Ta_R1 n_R1 DeltaH_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 DeltaH_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 DeltaH_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 DeltaH_R4 eta_a_R4 ...
       IFilter nFilter


%%AT
if(IFilter == 1)
    % Filering of reaction rates for R3 and R4
    for i = 1:nx_old
        psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
                + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
        Y_O2s   = max(Y_O2_olditer(i),0);
    
        RR3(i)  = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
                 *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/tempr1_olditer(i));
        RR4(i)  = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
                 *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/tempr1_olditer(i));
    end
    
    for n = 1:nFilter
        RR3bar(1) = 0.5*RR3(1) + 0.5*RR3(2);
        for i=2:(nx_old-1)
            RR3bar(i) = 0.25*RR3(i-1) + 0.5*RR3(i) + 0.25*RR3(i+1);
        end
        RR3bar(nx_old) = 0.5*RR3(nx_old-1) + 0.5*RR3(nx_old);
    
        RR3 = RR3bar;
    
        RR4bar(1) = 0.5*RR4(1) + 0.5*RR4(2);
        for i=2:(nx_old-1)
            RR4bar(i) = 0.25*RR4(i-1) + 0.5*RR4(i) + 0.25*RR4(i+1);
        end
        RR4bar(nx_old) = 0.5*RR4(nx_old-1) + 0.5*RR4(nx_old);
    
        RR4 = RR4bar;
    end
end
%%AT

for i = 1:nx_old
    
    % Porosity
    psi = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
        + psi_c * x_c_old(i) + psi_a * x_a_old(i);
        
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s   = rho_ws*c_ws*x_ws_old(i) ...
              + rho_ds*c_ds*x_ds_old(i) ...
              + rho_c *c_c * x_c_old(i) ...
              + rho_a *c_a * x_a_old(i); % (rho*cp) in solid phase 
    rhocp_g   = (pres_g*MW_g/R/temp_old(i))*cp_g0; % (rho*cp) in gas phase
    rhocp_eff = (1-psi)*rhocp_s + psi*rhocp_g; % Estimated at time t(n)

    % Volumetric rate of heat production/consumption [W/m3]
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s   = max(Y_O2_olditer(i),0);
    K_R1    = ( max(0,rho_ws*x_ws_olditer(i)*(1-psi_new))^n_R1 ) ...
              *A_R1*exp(-Ta_R1/tempr1_olditer(i));
    K_R2    = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R2 ) ...
              *A_R2*exp(-Ta_R2/tempr1_olditer(i));
    %%AT
    if(IFilter == 1)
        K_R3  = RR3(i);
        K_R4  = RR4(i);
    else
        K_R3  = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
               *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/tempr1_olditer(i));
        K_R4  = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
               *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/tempr1_olditer(i));
    end
    %%AT
    
    Qdotp = 0;
    Qdotm = 0;
    if( (1-eta_ds_R1)*DeltaH_R1 > 0 )
        Qdotp = Qdotp + (1-eta_ds_R1)*K_R1*DeltaH_R1;
    else
        Qdotm = Qdotm + (1-eta_ds_R1)*K_R1*DeltaH_R1;
    end
    if( (1-eta_c_R2) *DeltaH_R2 > 0 )
        Qdotp = Qdotp + (1-eta_c_R2) *K_R2*DeltaH_R2;
    else
        Qdotm = Qdotm + (1-eta_c_R2) *K_R2*DeltaH_R2;
    end
    if( (1-eta_c_R3) *DeltaH_R3 > 0 )
        Qdotp = Qdotp + (1-eta_c_R3) *K_R3*DeltaH_R3;
    else
        Qdotm = Qdotm + (1-eta_c_R3) *K_R3*DeltaH_R3;
    end
    if( (1-eta_a_R4) *DeltaH_R4 > 0 )
        Qdotp = Qdotp + (1-eta_a_R4) *K_R4*DeltaH_R4;
    else
        Qdotm = Qdotm + (1-eta_a_R4) *K_R4*DeltaH_R4;
    end
    Qdotp = Qdotp/rhocp_eff;
    Qdotm = Qdotm/rhocp_eff;
    
    tempr1_newiter(i) = ...
    ( tempr0_newiter(i) + 0.5*dt*Qdotp + 0.5*lambda*tempr1_olditer(i) ) ...
    /( 1 - 0.5*dt*Qdotm/tempr1_olditer(i) + 0.5*lambda );
    
end
            
end
