function dt=timeStep(FMC, A_R1, Ta_R1, A_R2, Ta_R2, rho_m, rho_vs, rho_c, ...
       k_m, k_vs, k_c, c_m, c_vs, c_c, dx_i )


         

% Select numerical parameters: time resolution such that FO ~ 0.5
% Calculate maximum value of heat diffusivity
alpha_max = 0;
for i=1:101
    x_vs_tmp = (i-1)/100;   % Cover range 0 <= x_vs_tmp <= 1
    for j=1:21
        FMC_tmp = (j-1)*FMC/20;   % Cover range 0 <= FMC_tmp <= FMC
        x_m_tmp = (rho_vs/rho_m)*FMC_tmp*x_vs_tmp;
        x_c_tmp = 1-x_vs_tmp-x_m_tmp;
    
        rhoc_p_tmp = rho_m*c_m*x_m_tmp + rho_vs*c_vs*x_vs_tmp ...
                                       + rho_c*c_c*x_c_tmp;
        k_p_tmp    = k_m*x_m_tmp + k_vs*x_vs_tmp +k_c*x_c_tmp;
    
        alpha_max = max(alpha_max,(k_p_tmp/rhoc_p_tmp));
    end
end
dt = 0.5*dx_i^2/alpha_max;   % Select time step (constant)
FO = alpha_max*dt/dx_i^2;
fprintf(' dt = %g \n',dt);
fprintf(' FO = %g \n',FO);

%
% Begin analysis of chemical time scales
tau_RR1_min = 1000;
for i=1:21
    tem_tmp = 300 + (i-1)*200/20; % Cover range 300 <= temp_tmp <= 500 K
    if(A_R1 ~= 0)
        tau_RR1 = 1/(A_R1*exp(-Ta_R1/tem_tmp));
        tau_RR1_min = min(tau_RR1_min,tau_RR1);
    end
end
fprintf(' \n');
fprintf(' tau_RR1_min, (tau_RR1_min/dt) = %g %g \n', ...
          tau_RR1_min,(tau_RR1_min/dt));
tau_RR2_min = 1000;
for i=1:101
    tem_tmp = 300 + (i-1)*700/100; % Cover range 300 <= temp_tmp <= 1000 K
    if(A_R2 ~= 0)
        tau_RR2 = 1/(A_R2*exp(-Ta_R2/tem_tmp));
        tau_RR2_min = min(tau_RR2_min,tau_RR2);
    end
end
fprintf(' tau_RR2_min, (tau_RR2_min/dt) = %g %g \n', ...
          tau_RR2_min,(tau_RR2_min/dt));
% Enforce (tau_RR1_min/dt) and (tau_RR2_min/dt) >=10
% scale = min((tau_RR1_min/dt),(tau_RR2_min/dt));
% if( scale < 10 )
%     dt = dt*0.1*scale;   % Modified time step
%     FO = alpha_max*dt/dx_i^2;
%     fprintf(' dt = %g \n',dt);
%     fprintf(' FO = %g \n',FO);
% end
% End analysis of chemical time scales
%

end