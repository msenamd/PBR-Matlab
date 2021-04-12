% Energy conservation
% - Semi-implicit time integration
%   (implicit treatment of diffusion, explicit treatment of chemisty)
% - Formulation as a tri-diagonal matrix system

function [a, b, c, d] = energy_conservation(dt,q_surf_old, ...
           temp_old,x_m_old,x_vs_old,x_c_old,...
           xRight,xCenter,dV,nx_old,geometry,A_rectangle,L_cylinder)

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       k_m k_vs k_c c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g

%constructung face area , NON SIGNED (Note: r=xRight for cylinders and spheres)   
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
   
  
   
for i=1:nx_old
    % Thermal conductivity
    kp(i) = k_m*x_m_old(i) + k_vs*x_vs_old(i) + k_c*x_c_old(i);
    
    % Product of mass density times specific heat
    % - Notations: rho_times_cp = rho_p x cp
    rho_times_cp(i) = rho_m*c_m*x_m_old(i) + rho_vs*c_vs*x_vs_old(i) ...
                                           + rho_c*c_c*x_c_old(i);

    % Volumetric rate of heat production/consumption [W/m3]                    
    Qdotp(i) = rho_m * x_m_old(i)*A_R1*exp(-Ta_R1/temp_old(i))*DeltaH_R1...
             + rho_vs*x_vs_old(i)*A_R2*exp(-Ta_R2/temp_old(i))*DeltaH_R2...
               *(1-eta_c) ...
             + rho_c * x_c_old(i)*x_O2_g ...
                                 *A_R3*exp(-Ta_R3/temp_old(i))*DeltaH_R3;                       
end
               
% i = 1 (back surface)
i = 1;
FO_L = 0;
FO_R = 0.5*(kp(i)+kp(i+1))/rho_times_cp(i) * S_pos(i)*dt / ...
                    (dV(i)*(xCenter(i+1)-xCenter(i)));
a(1) =-0.5*FO_L;
b(1) = 1 + 0.5*FO_L + 0.5*FO_R;
c(1) =-0.5*FO_R;
d(1) = (1-0.5*FO_R)*temp_old(1) + 0.5*FO_R*temp_old(2) ...
                                + dt/rho_times_cp(i)*Qdotp(1);

% i = nx_old (exposed surface)
i = nx_old;
FO_L = 0.5*(kp(i)+kp(i-1))/rho_times_cp(i) * S_neg(i)*dt / ...
                                     (dV(i)*(xCenter(i)-xCenter(i-1)));
FO_R = 0;
a(nx_old) =-0.5*FO_L;
b(nx_old) = 1 + 0.5*FO_L + 0.5*FO_R;
c(nx_old) =-0.5*FO_R;

d(nx_old) = 0.5*FO_L*temp_old(nx_old-1) + (1-0.5*FO_L)*temp_old(nx_old) ...
 + dt/rho_times_cp(i)*Qdotp(nx_old) + dt/rho_times_cp(i)*S_pos(nx_old)/dV(nx_old)*q_surf_old;


% 2 <= i <= (nx_old-1) (interior nodes)
for i=2:(nx_old-1)
    FO_L = 0.5*(kp(i)+kp(i-1))/rho_times_cp(i) * S_neg(i)*dt / ...
                                         (dV(i)*(xCenter(i)-xCenter(i-1)));
    FO_R = 0.5*(kp(i)+kp(i+1))/rho_times_cp(i) * S_pos(i)*dt / ...
                        (dV(i)*(xCenter(i+1)-xCenter(i)));
    a(i) =-0.5*FO_L;
    b(i) = 1 + 0.5*FO_L + 0.5*FO_R;
    c(i) =-0.5*FO_R;
    d(i) = 0.5*FO_L*temp_old(i-1)+(1 - 0.5*FO_L - 0.5*FO_R)*temp_old(i) ...
                                 + 0.5*FO_R*temp_old(i+1) + dt/rho_times_cp(i)*Qdotp(i);
end
            
end