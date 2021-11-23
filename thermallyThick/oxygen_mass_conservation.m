% Oxygen mass conservation
% - Implicit time integration
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
% - Formulation as a tri-diagonal matrix system

function [a, b, c, d] = oxygen_mass_conservation(dt, Y_O2_old, ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, pres_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, lambda, mdot_surf_old, h_conv, ...
    xRight, xCenter, dV_old, nx_old)
           

global geometry A_rectangle L_cylinder ...
       Y_g_O2 pres_g cp_g0 nu_g0 MW_g R ...
       rho_ws rho_ds rho_c rho_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 eta_c_R3 eta_O2_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 eta_a_R4 eta_O2_R4 ...
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
    % Porosity (estimated at time t(n))
    psi          = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
                 + psi_c * x_c_old(i) + psi_a * x_a_old(i);
    
    % Gas mass density (estimated at time t(n))
    rho_g        = (pres_g*MW_g/R/temp_old(i));
    
    % Gas mass diffusity (D_g = nu_g) (estimated at time t(n))
    D_g          = nu_g0*(temp_old(i)/300)^(1.76);
    
    % Product of mass density times porosity times diffusivity
    % (estimated at time t(n))
    rhopsiD_g(i) = psi*rho_g*D_g;
    
    % Product of mass density times porosity (estimated at time t(n))
    rhopsi_g(i)  = psi*rho_g;
    
    % Permeability (estimated at time t(n))
    Kperm(i)     = Kperm_ws*x_ws_old(i) + Kperm_ds*x_ds_old(i) ...
                 + Kperm_c * x_c_old(i) + Kperm_a * x_a_old(i);
    
    % Volumetric reaction rate [kg/s/m3]
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s   = max(Y_O2_olditer(i),0);
    K_R1    = ( max(0,rho_ws*x_ws_olditer(i)*(1-psi_new))^n_R1 ) ...
              *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2    = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R2 ) ...
              *A_R2*exp(-Ta_R2/temp_olditer(i));
    %%AT
    if(IFilter == 1)
        K_R3 = RR3(i);
        K_R4 = RR4(i);
    else
        K_R3 = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
              *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp_olditer(i));
        K_R4 = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
              *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp_olditer(i));
    end
    %%AT

    % Assume n_O2_R3 >= 1 and n_O2_R4 >=1
    if( abs(Y_O2s) > 0 )
        K_R3b = K_R3/Y_O2s;
        K_R4b = K_R4/Y_O2s;
    else
        K_R3b = 0;
        K_R4b = 0;
    end
    %%AT K_R3b  = ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
    %%AT          *( Y_O2s^(n_O2_R3-1) )*A_R3*exp(-Ta_R3/temp_olditer(i));
    %%AT K_R4b  = ( max(0,rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
    %%AT          *( Y_O2s^(n_O2_R4-1) )*A_R4*exp(-Ta_R4/temp_olditer(i));
    
    xstore = (1-eta_ds_R1)*K_R1  + (1-eta_c_R2) *K_R2 ...
           + (1-eta_c_R3) *K_R3  + (1-eta_a_R4) *K_R4 ...    
           + eta_O2_R3    *K_R3b + eta_O2_R4    *K_R4b;
          
    bRRO2(i) = xstore*dt/rho_g/psi;
end
               
% i = 1 (back surface)
i    = 1;
FO_L = 0;
FO_R = 0.5*(rhopsiD_g(i)+rhopsiD_g(i+1))/rhopsi_g(i) * S_pos(i)*dt / ...
                                     (dV_old(i)*(xCenter(i+1)-xCenter(i)));

nu_g = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
CFL_L = 0;
if( pres_old(i) >= pres_old(i+1) )
    CFL_R = 0;
else
    CFL_R = dt*(Kperm(i)/nu_g/rhopsi_g(i)) ...
            *(pres_old(i)-pres_old(i+1))/(xCenter(i+1)-xCenter(i))^2;
    % Note: CFL_R < 0
end

a(i) = 0;
b(i) = 1 + 0.5*FO_R - 0.5*CFL_R ...
     + bRRO2(i);
c(i) =-0.5*FO_R + 0.5*CFL_R;
d(i) = (1 - 0.5*FO_R + 0.5*CFL_R)*Y_O2_old(i) ...
     + (0.5*FO_R - 0.5*CFL_R)*Y_O2_old(i+1);
%%AT
% Apply under-relaxation in iterative scheme
%%AT d(i) = d(i) + (1-lambda)*b(i)*Y_O2_olditer(i)/lambda;
%%AT b(i) = b(i)/lambda;
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*Y_O2_olditer(i);
%%AT

% i = nx_old (exposed surface)
i    = nx_old;
FO_L = 0.5*(rhopsiD_g(i)+rhopsiD_g(i-1))/rhopsi_g(i) * S_neg(i)*dt / ...
                                     (dV_old(i)*(xCenter(i)-xCenter(i-1)));
FO_R = 0;

nu_g = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
if( pres_old(i-1) > pres_old(i) )
    CFL_L = dt*(Kperm(i)/nu_g/rhopsi_g(i)) ...
            *(pres_old(i-1)-pres_old(i))/(xCenter(i)-xCenter(i-1))^2;
    % Note: CFL_L > 0
else
    CFL_L = 0;
end
CFL_R = 0;

a(i) =-0.5*FO_L -0.5*CFL_L;
b(i) = 1 + 0.5*FO_L + 0.5*CFL_L ...
     + bRRO2(i);
c(i) = 0;
d(i) = (0.5*FO_L + 0.5*CFL_L)*Y_O2_old(i-1) ...
     + (1 - 0.5*FO_L - 0.5*CFL_L)*Y_O2_old(i);
%%AT
% Apply under-relaxation in iterative scheme
%%AT d(i) = d(i) + (1-lambda)*b(i)*Y_O2_olditer(i)/lambda;
%%AT b(i) = b(i)/lambda;
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*Y_O2_olditer(i);
%%AT

% Contribution of gas-to-sold oxygen mass flux
rhopsi_surf  = rhopsi_g(nx_old);
rhopsiD_surf = rhopsiD_g(nx_old);
dx_surf      = xRight(nx_old)-xCenter(nx_old);
h_mass       = h_conv/cp_g0;
Bi           = h_mass*dx_surf/rhopsiD_surf;

b(i) = b(i) + ( h_mass/(1+Bi) ) ...
                                  *0.5*dt/rhopsi_surf * S_pos(i)/dV_old(i);
d(i) = d(i) + ( mdot_surf_old + (h_mass*Y_g_O2/(1+Bi)) ) ...
                                  *0.5*dt/rhopsi_surf * S_pos(i)/dV_old(i);

% 2 <= i <= (nx_old-1) (interior nodes)
for i = 2:(nx_old-1)   
    FO_L = 0.5*(rhopsiD_g(i)+rhopsiD_g(i-1))/rhopsi_g(i) * S_neg(i)*dt /...
                                  (dV_old(i)*(xCenter(i)-xCenter(i-1)));
    FO_R = 0.5*(rhopsiD_g(i)+rhopsiD_g(i+1))/rhopsi_g(i) * S_pos(i)*dt /...
                                  (dV_old(i)*(xCenter(i+1)-xCenter(i)));
    
    nu_g = nu_g0*(temp_old(i)/300)^(1.76); % Gas kinematic viscosity
    if( (pres_old(i-1) > pres_old(i)) & (pres_old(i) > pres_old(i+1)) )
        % Case for which mdotPUA > 0. Note: CFL_L > 0
        CFL_L = dt*(Kperm(i)/nu_g/rhopsi_g(i)) ...
                *(pres_old(i-1)-pres_old(i))/(xCenter(i)-xCenter(i-1))^2;
        CFL_R = 0;
    elseif ( (pres_old(i-1) < pres_old(i)) ...
                                          & (pres_old(i) < pres_old(i+1)) )
        % Case for which mdotPUA < 0. Note: CFL_R < 0
        CFL_L = 0;
        CFL_R = dt*(Kperm(i)/nu_g/rhopsi_g(i)) ...
                *(pres_old(i)-pres_old(i+1))/(xCenter(i+1)-xCenter(i))^2;
    else
        % Case for which dp/dx ~ 0
        CFL_L = 0;
        CFL_R = 0;
    end
        
    a(i) =-0.5*FO_L -0.5*CFL_L;
    b(i) = 1 + 0.5*FO_L + 0.5*FO_R + 0.5*CFL_L - 0.5*CFL_R ...
         + bRRO2(i);
    c(i) =-0.5*FO_R + 0.5*CFL_R;
    d(i) = (0.5*FO_L + 0.5*CFL_L)*Y_O2_old(i-1) ...
         + (1 - 0.5*FO_L - 0.5*FO_R - 0.5*CFL_L + 0.5*CFL_R)*Y_O2_old(i)...
         + (0.5*FO_R - 0.5*CFL_R)*Y_O2_old(i+1);
    %%AT
    % Apply under-relaxation in iterative scheme
    %%AT d(i) = d(i) + (1-lambda)*b(i)*Y_O2_olditer(i)/lambda;
    %%AT b(i) = b(i)/lambda;
    b(i) = b(i) + lambda;
    d(i) = d(i) + lambda*Y_O2_olditer(i);
    %%AT
end
            
end