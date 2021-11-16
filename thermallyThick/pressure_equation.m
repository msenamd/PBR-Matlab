% Pressure equation
% - Semi-implicit time integration
%
%   (pres(n+1)-pres(n))/dt = RHS1 + RHS2 where RHS1 is pressure change
%                            due to reactions (R1)-(R4) and RHS2 is
%                            pressure change due to diffusion
%
%   Evaluation of RHS1:
%      RHS1 = RHS1(n)
%   Evaluation of RHS2:
%      RHS2 = 0.5*RHS2(n+1) + 0.5*RHS2(n)
%
% - Formulation as a tri-diagonal matrix system
%
% - Definitions
%      pres_old     = pres(n)

function [a, b, c, d] = pressure_equation(dt, pres_old, ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, pres_olditer, lambda_pres, ...
    xRight, xCenter, dV_old, nx_old)

   

global geometry A_rectangle L_cylinder ...
       nu_g0 MW_g R ...
       rho_ws rho_ds rho_c rho_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 eta_a_R4

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
    % Porosity (estimated at time t(n))
    psi        = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
               + psi_c * x_c_old(i) + psi_a * x_a_old(i);
    
    % Permeability (estimated at time t(n))
    Kperm      = Kperm_ws*x_ws_old(i) + Kperm_ds*x_ds_old(i) ...
               + Kperm_c * x_c_old(i) + Kperm_a * x_a_old(i);
    
    % Gas kinematic viscosity (estimated at time t(n))
    nu_g       = nu_g0*(temp_old(i)/300)^(1.76);
    
    xstore1(i) = Kperm/nu_g;
    xstore2(i) = MW_g*psi/R/temp_old(i);
    
    % Volumetric reaction rate [kg/s/m3]
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s  = max(Y_O2_olditer(i),0);
    K_R1   = ( (rho_ws*x_ws_olditer(i)*(1-psi_new))^n_R1 ) ...
             *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2   = ( (rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R2 ) ...
             *A_R2*exp(-Ta_R2/temp_olditer(i));
    K_R3   = ( (rho_ds*x_ds_olditer(i)*(1-psi_new))^n_R3 ) ...
             *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp_olditer(i));
    K_R4   = ( (rho_c *x_c_olditer(i) *(1-psi_new))^n_R4 ) ...
             *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp_olditer(i));
         
    Mdotg  = (1-eta_ds_R1)*K_R1  + (1-eta_c_R2) *K_R2 ...
           + (1-eta_c_R3) *K_R3  + (1-eta_a_R4) *K_R4;  
                 
    RRp(i) = Mdotg*dt*temp_old(i)*R/MW_g/psi;  
end
               
% i = 1 (back surface)
i    = 1;
FO_L = 0;
FO_R = 0.5*(xstore1(i)+xstore1(i+1))/xstore2(i) * S_pos(i)*dt / ...
                                     (dV_old(i)*(xCenter(i+1)-xCenter(i)));

a(i) = 0;
b(i) = 1 + 0.5*FO_R;
c(i) =-0.5*FO_R;
d(i) = (1 - 0.5*FO_R)*pres_old(i) ...
     + (0.5*FO_R)*pres_old(i+1) ...
     + RRp(i);
%%AT
% Apply under-relaxation in iterative scheme
%%AT d(i) = d(i) + (1-lambda_pres)*b(i)*pres_olditer(i)/lambda_pres;
%%AT b(i) = b(i)/lambda_pres;
b(i) = b(i) + lambda_pres;
d(i) = d(i) + lambda_pres*pres_olditer(i);
%%AT
     
% i = nx_old (exposed surface)
i       = nx_old;
FO_L    = 0.5*(xstore1(i)+xstore1(i-1))/xstore2(i) * S_neg(i)*dt / ...
                                     (dV_old(i)*(xCenter(i)-xCenter(i-1)));
dx_surf = xRight(nx_old)-xCenter(nx_old);
FO_R    = xstore1(i)/xstore2(i) * S_pos(i)*dt / (dV_old(i)*dx_surf);

a(i) =-0.5*FO_L;
b(i) = 1 + 0.5*FO_L + 0.5*FO_R;
c(i) = 0;
d(i) = (0.5*FO_L)*pres_old(i-1) ...
     + (1 - 0.5*FO_L - 0.5*FO_R)*pres_old(i) ...
     + (FO_R)*0 ...
     + RRp(i);

%%AT
% Apply under-relaxation in iterative scheme
%%AT d(i) = d(i) + (1-lambda_pres)*b(i)*pres_olditer(i)/lambda_pres;
%%AT b(i) = b(i)/lambda_pres;
b(i) = b(i) + lambda_pres;
d(i) = d(i) + lambda_pres*pres_olditer(i);
%%AT

% 2 <= i <= (nx_old-1) (interior nodes)
for i = 2:(nx_old-1)
    FO_L = 0.5*(xstore1(i)+xstore1(i-1))/xstore2(i) * S_neg(i)*dt / ...
                                     (dV_old(i)*(xCenter(i)-xCenter(i-1)));
    FO_R = 0.5*(xstore1(i)+xstore1(i+1))/xstore2(i) * S_pos(i)*dt / ...
                                     (dV_old(i)*(xCenter(i+1)-xCenter(i)));
        
    a(i) =-0.5*FO_L;
    b(i) = 1 + 0.5*FO_L + 0.5*FO_R;
    c(i) =-0.5*FO_R;
    d(i) = (0.5*FO_L)*pres_old(i-1) ...
         + (1 - 0.5*FO_L - 0.5*FO_R)*pres_old(i) ...
         + (0.5*FO_R)*pres_old(i+1) ...
         + RRp(i);
    %%AT
    % Apply under-relaxation in iterative scheme
    %%AT d(i) = d(i) + (1-lambda_pres)*b(i)*pres_olditer(i)/lambda_pres;
    %%AT b(i) = b(i)/lambda_pres;
    b(i) = b(i) + lambda_pres;
    d(i) = d(i) + lambda_pres*pres_olditer(i);
    %%AT
end
            
end
