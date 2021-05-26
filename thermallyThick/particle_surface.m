% Calculate the net surface heat flux and the surface temperature
% - at exposed surface (x = delta)

function [q_surf, temp_surf] = ...
          particle_surface(temp_surf_old,kp_surf,eps_surf,dx_surf,tempc_nx,h_conv)
                      
global T_g G

sigma = 5.67e-8;   % Stefan-Boltzmann constant [W/m2/K4]

iter  = 0;
small = 1;
maxIter = 100;

T_old = temp_surf_old;   % Initial guess for particle surface temperature
while (small > 1e-3 && iter <= maxIter)
    iter  = iter+1;
    h_rad = eps_surf*sigma*(T_old^3);
    Bi    = (h_conv+h_rad)*dx_surf/kp_surf;
    T_new = (tempc_nx + h_conv*T_g*dx_surf/kp_surf + (eps_surf*G*dx_surf/kp_surf))/(1+Bi);
    
    small = abs(T_new-T_old);
    T_old = T_new;
    %fprintf(' *** iter,small = %g %g \n',iter,small);
end
if(iter >= maxIter)
    fprintf(' *** Error in particle_surface, iterations exceed the maximum limit, temp_surf = %g \n',T_new);
    return
end

temp_surf = T_new;
q_surf    = eps_surf*G - eps_surf*sigma*temp_surf^4 ...
                       + h_conv*(T_g-temp_surf);

end