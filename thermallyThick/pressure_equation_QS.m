% Pressure equation with Quasi-Steady model
% - Semi-implicit time integration
%
%   0 = RHS1 + RHS2 where RHS1 is pressure change due to reactions
%                   (R1)-(R4) and RHS2 is pressure change due to convection
%
% - Formulation as a tri-diagonal matrix system

function [a, b, c, d] = pressure_equation_QS( ...
    temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, pres_olditer, lambda, ...
    sum_R1, sum_R2, sum_R3, sum_R4, ...
    xRight, xCenter, dV_old, nx_old)


global geometry A_rectangle L_cylinder ...
       nu_g0 ...
       rho_ws rho_ds rho_c rho_a ...
       psi_ws psi_ds psi_c psi_a ...
       Kperm_ws Kperm_ds Kperm_c Kperm_a ...
       A_R1 Ta_R1 n_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 eta_a_R4 ...
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
        Y_O2s   = max(Y_O2_olditer(i),0);

        RR3(i)  = ...
           ( max(0,rho_ds*(1-psi_ds)*x_ds_olditer(i)*dV_old(i))^n_R3 ) ...
              *sum_R3(i)^(1-n_R3) ...
              *( (Y_O2s/0.226)^n_O2_R3 ) ...
              *A_R3*exp(-Ta_R3/temp_olditer(i));
        RR4(i)  = ...
           ( max(0,rho_c *(1-psi_c) *x_c_olditer(i) *dV_old(i))^n_R4 ) ...
              *sum_R4(i)^(1-n_R4) ...
              *( (Y_O2s/0.226)^n_O2_R4 ) ...
              *A_R4*exp(-Ta_R4/temp_olditer(i));
              %%AT *A_R4*exp(-Ta_R4/min(temp_olditer(i),700));
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
    % Permeability (estimated at time t(n))
    Kperm      = Kperm_ws*x_ws_old(i) + Kperm_ds*x_ds_old(i) ...
               + Kperm_c * x_c_old(i) + Kperm_a * x_a_old(i);
    
    % Gas kinematic viscosity (estimated at time t(n))
    nu_g       = nu_g0*(temp_old(i)/300)^(1.76);
    
    xstore1(i) = Kperm/nu_g;
    
    % Volumetric reaction rate [kg/s/m3]
    Y_O2s   = max(Y_O2_olditer(i),0);
    K_R1  = ...
         ( max(0,rho_ws*(1-psi_ws)*x_ws_olditer(i)*dV_old(i))^n_R1 ) ...
            *sum_R1(i)^(1-n_R1) ...
            *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2  = ...
         ( max(0,rho_ds*(1-psi_ds)*x_ds_olditer(i)*dV_old(i))^n_R2 ) ...
            *sum_R2(i)^(1-n_R2) ...
            *A_R2*exp(-Ta_R2/temp_olditer(i));
    %%AT
    if(IFilter == 1)
        K_R3 = RR3(i);     
        K_R4 = RR4(i);
    else
        K_R3 = ...
           ( max(0,rho_ds*(1-psi_ds)*x_ds_olditer(i)*dV_old(i))^n_R3 ) ...
              *sum_R3(i)^(1-n_R3) ...
              *( (Y_O2s/0.226)^n_O2_R3 ) ...
              *A_R3*exp(-Ta_R3/temp_olditer(i));
        K_R4 = ...
           ( max(0,rho_c *(1-psi_c) *x_c_olditer(i) *dV_old(i))^n_R4 ) ...
              *sum_R4(i)^(1-n_R4) ...
              *( (Y_O2s/0.226)^n_O2_R4 ) ...
              *A_R4*exp(-Ta_R4/temp_olditer(i));
              %%AT *A_R4*exp(-Ta_R4/min(temp_olditer(i),700));
    end
    %%AT
         
    Mdotg = (1-eta_ds_R1)*K_R1 + (1-eta_c_R2) *K_R2 ...
          + (1-eta_c_R3) *K_R3 + (1-eta_a_R4) *K_R4;
                 
    dRR(i) = Mdotg;
end
               
% i = 1 (back surface)
i    = 1;
a(i) = 0;
c(i) = -0.5*(xstore1(i+1)+xstore1(i)) ...
           *S_pos(i)/(xCenter(i+1)-xCenter(i));
b(i) = -c(i);
d(i) = dRR(i);

% Apply under-relaxation in iterative scheme
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*pres_olditer(i);
     
% i = nx_old (exposed surface)
i       = nx_old;
a(i)    = -0.5*(xstore1(i)+xstore1(i-1)) ...
              *S_neg(i)/(xCenter(i)-xCenter(i-1));   
dx_surf = xRight(nx_old)-xCenter(nx_old);
b(i)    = -a(i)+xstore1(i) *S_pos(i)/dx_surf;
c(i)    = 0;
d(i)    = (xstore1(i) *S_pos(i)/dx_surf)*0 + dRR(i);

% Apply under-relaxation in iterative scheme
b(i) = b(i) + lambda;
d(i) = d(i) + lambda*pres_olditer(i);

% 2 <= i <= (nx_old-1) (interior nodes)
for i = 2:(nx_old-1)    
    a(i) = -0.5*(xstore1(i)+xstore1(i-1)) ...
               *S_neg(i)/(xCenter(i)-xCenter(i-1));   
    c(i) = -0.5*(xstore1(i+1)+xstore1(i)) ...
               *S_pos(i)/(xCenter(i+1)-xCenter(i));
    b(i) = -a(i)-c(i);
    d(i) = dRR(i);

    % Apply under-relaxation in iterative scheme
    b(i) = b(i) + lambda;
    d(i) = d(i) + lambda*pres_olditer(i);
end
            
end
