% Energy conservation
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
% - Formulation as a tri-diagonal matrix system

function [a, b, c, d, temp_surf_newiter] = ...
             energy_conservation(dt, temp_old, temp_olditer, ...
             x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
             x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
             Y_O2_olditer, temp_surf_olditer, h_conv, ...
             lambda, xRight, xCenter, dV_old, nx_old)

   
global geometry A_rectangle L_cylinder ...
       T_g G pres_g k_g0 cp_g0 nu_g0 MW_g R sigma ...  
       rho_ws rho_ds rho_c rho_a ...
       k_ws k_ds k_c k_a ...
       c_ws c_ds c_c c_a ...
       eps_ws eps_ds eps_c eps_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 DeltaH_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 DeltaH_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 DeltaH_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 DeltaH_R4 eta_a_R4 ...
       IFilter nFilter
   
   
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

%%AT
if(IFilter == 1)
    % Filering of reaction rates for R3 and R4
    for i = 1:nx_old
        psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
                + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
        Y_O2s   = max(Y_O2_olditer(i),0);
    
        RR3(i)  = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
                 *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp_olditer(i));
        RR4(i)  = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
                 *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp_olditer(i));
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

    % Volumetric rate of heat production/consumption [W/m3]
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s   = max(Y_O2_olditer(i),0);
    K_R1    = ( max(0,rho_ws*x_ws_olditer(i)*(1-psi_new))^n_R1 ) ...
              *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2    = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R2 ) ...
              *A_R2*exp(-Ta_R2/temp_olditer(i));
    %%AT
    if(IFilter == 1)
        K_R3  = RR3(i);
        K_R4  = RR4(i);
    else
        K_R3  = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
               *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp_olditer(i));
        K_R4  = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
               *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp_olditer(i));
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
    dRRTemp(i) = Qdotp*dt/rhocp_eff(i);
    bRRTemp(i) =-Qdotm*dt/rhocp_eff(i)/temp_olditer(i);
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
b(i) = 1 + 0.5*FO_R - 0.5*CFL_R ...
     + bRRTemp(i);
c(i) =-0.5*FO_R + 0.5*CFL_R;
d(i) = (1 - 0.5*FO_R + 0.5*CFL_R)*temp_old(i) ...
     + (0.5*FO_R - 0.5*CFL_R)*temp_old(i+1) ...
     + dRRTemp(i);
 
%%AT
% Apply under-relaxation in iterative scheme
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*temp_olditer(i);
%%AT

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
b(i) = 1 + 0.5*FO_L + 0.5*CFL_L ...
     + bRRTemp(i);
c(i) = 0;
d(i) = (0.5*FO_L + 0.5*CFL_L)*temp_old(i-1) ...
     + (1 - 0.5*FO_L - 0.5*CFL_L)*temp_old(i) ...
     + dRRTemp(i);

%%AT
% Apply under-relaxation in iterative scheme
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*temp_olditer(i);
%%AT

% Contribution of gas-to-sold heat flux
keff_surf = keff(nx_old);
eps_surf  = eps_ws*x_ws_old(nx_old) + eps_ds*x_ds_old(nx_old) ...
          + eps_c * x_c_old(nx_old) + eps_a * x_a_old(nx_old);
dx_surf   = xRight(nx_old)-xCenter(nx_old);
%%AT h_rad     = eps_surf*sigma*(temp_surf_old^3);
%%AT Bi        = (h_conv+h_rad)*dx_surf/keff_surf;
h_rad     = eps_surf*sigma*(temp_surf_olditer^3);
Bi        = (h_conv+h_rad)*dx_surf/keff_surf;
tempc_nx  = temp_olditer(nx_old);
temp_surf_newiter = ( tempc_nx ...
                    + (h_conv*T_g+eps_surf*G)*dx_surf/keff_surf )/(1+Bi);
%%AT

%%AT
%{
b(i)      = b(i) + ((h_conv+h_rad)/(1+Bi)) ...
                                 *0.5*dt/rhocp_eff(i) * S_pos(i)/dV_old(i);
d(i)      = d(i) + ( q_surf_old ...
                        + ((h_conv*T_g+eps_surf*G)/(1+Bi)) ) ...
                                 *0.5*dt/rhocp_eff(i) * S_pos(i)/dV_old(i);
%}
b(i) = b(i) + ((h_conv+h_rad)/(1+Bi)) ...
                                     *dt/rhocp_eff(i) * S_pos(i)/dV_old(i);
d(i) = d(i) + ((h_conv*T_g+eps_surf*G)/(1+Bi)) ...
                                     *dt/rhocp_eff(i) * S_pos(i)/dV_old(i);
%%AT

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
    b(i) = 1 + 0.5*FO_L + 0.5*FO_R + 0.5*CFL_L - 0.5*CFL_R ...
         + bRRTemp(i);
    c(i) =-0.5*FO_R + 0.5*CFL_R;
    d(i) = (0.5*FO_L + 0.5*CFL_L)*temp_old(i-1) ...
         + (1 - 0.5*FO_L - 0.5*FO_R - 0.5*CFL_L + 0.5*CFL_R) ...
           *temp_old(i) ...
         + (0.5*FO_R - 0.5*CFL_R)*temp_old(i+1) ...
         + dRRTemp(i);
     %%AT
    % Apply under-relaxation in iterative scheme
    b(i) = b(i) + lambda;
    d(i) = d(i) + lambda*temp_olditer(i);
    %%AT
end
            
end
