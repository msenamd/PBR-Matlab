function dt = timeStep(rho_ws, rho_ds, rho_c, rho_a, ...
                    k0_ws, k0_ds, k0_c, k0_a, nk_ws, nk_ds, nk_c, nk_a, ...
                    gamma_ws, gamma_ds, gamma_c, gamma_a, ...
                    c0_ws, c0_ds, c0_c, c0_a, nc_ws, nc_ds, nc_c, nc_a, ...
                       psi_ws, psi_ds, psi_c, psi_a, ...
                       Kperm_ws, Kperm_ds, Kperm_c, Kperm_a, ...
                       A_R1, Ta_R1, ...
                       A_R2, Ta_R2, ...
                       A_R3, Ta_R3, n_O2_R3, ...
                       A_R4, Ta_R4, n_O2_R4, ...
                       Y_g_O2, nu_g0, MW_g, R, dx_i)
         

% Estimate dt based on heat diffusion inside solid phase
%  - Enforce FO = (k_s/rho_s c_s)*dt/dx^2 ~ 0.5
alpha_max = max( [(k0_ws/rho_ws/c0_ws),(k0_ds/rho_ds/c0_ds), ...
                  (k0_c /rho_c /c0_c ),(k0_a /rho_a /c0_a )] );
dt        = 0.5*dx_i^2/alpha_max;
dt_SPDiff = dt;
FO        = alpha_max*dt/dx_i^2;
fprintf(' dt                  = %g \n',dt);
fprintf(' FO (solid species)  = %g \n',FO);

% Estimate dt based on oxygen diffusion inside gas phase
%  - Enforce FO = D_g*dt/dx^2 ~ 0.5
D_g_max = 0;
for i = 1:101
    tem_tmp = 300 + (i-1)*900/100; % Cover range 300 <= temp <= 1200 K
    D_g     = nu_g0*(tem_tmp/300)^(1.76);
    D_g_max = max(D_g_max, D_g);
end
%dt = 0.5*dx_i^2/D_g_max;
FO = D_g_max*dt/dx_i^2;
fprintf(' \n');
fprintf(' dt                  = %g \n',dt);
fprintf(' FO (gaseous oxygen) = %g \n',FO);

% Estimate dt based on pressure diffusion inside gas phase
%  - Enforce FO = D_p*dt/dx^2 ~ 0.5
D_p_max = 0;
for i = 1:101
    tem_tmp = 300 + (i-1)*900/100; % Cover range 300 <= temp <= 1200 K
    nu_g    = nu_g0*(tem_tmp/300)^(1.76);
    D_p     = max( [(Kperm_ws/nu_g)/(MW_g*psi_ws/R/tem_tmp), ...
                    (Kperm_ds/nu_g)/(MW_g*psi_ds/R/tem_tmp), ...
                    (Kperm_c /nu_g)/(MW_g*psi_c /R/tem_tmp), ...
                    (Kperm_a /nu_g)/(MW_g*psi_a /R/tem_tmp)] );
    D_p_max = max(D_p_max, D_p);
end
%dt = 0.5*dx_i^2/D_p_max;
FO = D_p_max*dt/dx_i^2;
fprintf(' \n');
fprintf(' dt                  = %g \n',dt);
fprintf(' FO (gas pressure)   = %g \n',FO);

% Estimate dt based on heterogeneous chemistry
tau_RR1_min = 1000;
if(A_R1 ~= 0)
    for i = 1:21
        tem_tmp = 300 + (i-1)*200/20; % Cover range 300 <= temp <= 500 K
        tau_RR1 = 1/(A_R1*exp(-Ta_R1/tem_tmp));
        tau_RR1_min = min(tau_RR1_min,tau_RR1);
    end
end
fprintf(' \n');
fprintf(' tau_RR1_min,(tau_RR1_min/dt) = %g %g \n', ...
          tau_RR1_min,(tau_RR1_min/dt));

tau_RR2_min = 1000;
if(A_R2 ~= 0)
    for i = 1:101
        tem_tmp = 300 + (i-1)*900/100; % Cover range 300 <= temp <= 1200 K
        tau_RR2 = 1/(A_R2*exp(-Ta_R2/tem_tmp));
        tau_RR2_min = min(tau_RR2_min,tau_RR2);
    end
end
fprintf(' tau_RR2_min,(tau_RR2_min/dt) = %g %g \n', ...
          tau_RR2_min,(tau_RR2_min/dt));

tau_RR3_min = 1000;
if( (Y_g_O2 ~= 0) & (A_R3 ~= 0) )    
    for i = 1:101
        tem_tmp = 300 + (i-1)*900/100; % Cover range 300 <= temp <= 1200 K
        tau_RR3 = 1/(A_R3*exp(-Ta_R3/tem_tmp));
        tau_RR3_min = min(tau_RR3_min,tau_RR3);
    end
end
fprintf(' tau_RR3_min,(tau_RR3_min/dt) = %g %g \n', ...
          tau_RR3_min,(tau_RR3_min/dt));      

tau_RR4_min = 1000;
if( (Y_g_O2 ~= 0) & (A_R4 ~= 0) )
    for i = 1:101
        tem_tmp = 300 + (i-1)*1200/100; % Cover range 300 <= temp <= 1500 K
        tau_RR4 = 1/(A_R4*exp(-Ta_R4/min(tem_tmp,700)));
        tau_RR4_min = min(tau_RR4_min,tau_RR4);
    end
end
fprintf(' tau_RR4_min,(tau_RR4_min/dt) = %g %g \n', ...
          tau_RR4_min,(tau_RR4_min/dt));
      
% Enforce (tau_RRj_min/dt) >=10
scale = min( [(tau_RR1_min/dt),(tau_RR2_min/dt), ...
              (tau_RR3_min/dt),(tau_RR4_min/dt)] );
if( scale < 10 )
    dt = dt*0.1*scale;   % Modified time step
    FO = alpha_max*dt/dx_i^2;
    fprintf(' dt_RR = %g \n',dt);
    fprintf(' FO    = %g \n',FO);
end

% Choose initial value of dt
dt = dt_SPDiff;
fprintf(' dt                  = %g \n',dt);
fprintf('\n');
    

end