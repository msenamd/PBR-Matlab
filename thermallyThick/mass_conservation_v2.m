% Mass conservation v2
% - Semi-implicit time integration

function [x_ws, x_ds, x_c, x_a, dV] = ...
    mass_conservation_v2(dt, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, sum_R1, sum_R2, sum_R3, sum_R4, ...
    lambda, dV_old, nx_old)


global rho_ws rho_ds rho_c rho_a ...
       psi_ws psi_ds psi_c psi_a ...
       A_R1 Ta_R1 n_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 eta_a_R4 ...
       IFilter nFilter


%%AT
if(IFilter == 1)
    % Filering of reaction rates for R3 and R4
    for i = 1:nx_old
        psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
                + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);           
        Y_O2s   = max(Y_O2_olditer(i),0);
        
        RR3(i)  = ...
           ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new)*dV_old(i))^n_R3 ) ...
                 *sum_R3(i)^(1-n_R3) ...
                 *( (1+Y_O2s)^n_O2_R3 - 1 ) ...
                 *A_R3*exp(-Ta_R3/temp_olditer(i));
        RR4(i)  = ...
           ( max(0,rho_c *x_c_olditer(i) *(1-psi_new)*dV_old(i))^n_R4 ) ...
                 *sum_R4(i)^(1-n_R4) ...
                 *( (1+Y_O2s)^n_O2_R4 - 1 ) ...
                 *A_R4*exp(-Ta_R4/min(temp_olditer(i),700));
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
    psi_old = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
            + psi_c * x_c_old(i) + psi_a * x_a_old(i);
        
    % Species mass conservation statements
    % - Notations:  xws_1mpsi_dV = x_ws x (1-psi) x dV
    %               xds_1mpsi_dV = x_ds x (1-psi) x dV
    %               xc_1mpsi_dV  = x_c  x (1-psi) x dV
    %               xa_1mpsi_dV  = x_a  x (1-psi) x dV
    %               onempsi_dV   = (1-psi) x dV
          
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s   = max(Y_O2_olditer(i),0);
    K_R1    = ...
           ( max(0,rho_ws*x_ws_olditer(i)*(1-psi_new)*dV_old(i))^n_R1 ) ...
              *sum_R1(i)^(1-n_R1) ...
              *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2    = ...
           ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new)*dV_old(i))^n_R2 ) ...
              *sum_R2(i)^(1-n_R2) ...
              *A_R2*exp(-Ta_R2/temp_olditer(i));          
    %%AT
    if(IFilter == 1)
        K_R3 = RR3(i);     
        K_R4 = RR4(i);
    else
        K_R3 = ...
           ( max(0,rho_ds*x_ds_olditer(i)*(1-psi_new)*dV_old(i))^n_R3 ) ...
              *sum_R3(i)^(1-n_R3) ...
              *( (1+Y_O2s)^n_O2_R3 - 1 ) ...
              *A_R3*exp(-Ta_R3/temp_olditer(i));
        K_R4 = ...
           ( max(0,rho_c *x_c_olditer(i) *(1-psi_new)*dV_old(i))^n_R4 ) ...
              *sum_R4(i)^(1-n_R4) ...
              *( (1+Y_O2s)^n_O2_R4 - 1 ) ...
              *A_R4*exp(-Ta_R4/min(temp_olditer(i),700));
    end
    %%AT
    
    xws_1mpsi_dV_old     = ( x_ws_old(i)    *(1-psi_old)*dV_old(i) );
    xws_1mpsi_dV_olditer = ( x_ws_olditer(i)*(1-psi_new)*dV_old(i) );
    RHSp                 =  0;
    RHSm                 = -K_R1/rho_ws;
    
    xws_1mpsi_dV_newiter = ...
        ( xws_1mpsi_dV_old + dt*RHSp + lambda*xws_1mpsi_dV_olditer ) ...
        /( 1 - dt*RHSm/max(xws_1mpsi_dV_olditer,1e-14) + lambda );
    
    xds_1mpsi_dV_old     = ( x_ds_old(i)    *(1-psi_old)*dV_old(i) );
    xds_1mpsi_dV_olditer = ( x_ds_olditer(i)*(1-psi_new)*dV_old(i) );
    RHSp                 =  K_R1*eta_ds_R1/rho_ds;
    RHSm                 = -(K_R2+K_R3)/rho_ds;
    
    xds_1mpsi_dV_newiter = ...
        ( xds_1mpsi_dV_old + dt*RHSp + lambda*xds_1mpsi_dV_olditer ) ...
        /( 1 - dt*RHSm/max(xds_1mpsi_dV_olditer,1e-14) + lambda );
    
    xc_1mpsi_dV_old      = ( x_c_old(i)     *(1-psi_old)*dV_old(i) );
    xc_1mpsi_dV_olditer  = ( x_c_olditer(i) *(1-psi_new)*dV_old(i) );
    RHSp                 =  K_R2*eta_c_R2/rho_c ...
                          + K_R3*eta_c_R3/rho_c;
    RHSm                 = -K_R4/rho_c;
    
    xc_1mpsi_dV_newiter = ...
        ( xc_1mpsi_dV_old + dt*RHSp + lambda*xc_1mpsi_dV_olditer ) ...
        /( 1 - dt*RHSm/max(xc_1mpsi_dV_olditer,1e-14) + lambda );
    
    
    xa_1mpsi_dV_old      = ( x_a_old(i)     *(1-psi_old)*dV_old(i) );
    xa_1mpsi_dV_olditer  = ( x_a_olditer(i) *(1-psi_new)*dV_old(i) );
    RHSp                 =  K_R4*eta_a_R4/rho_a;
    RHSm                 =  0;
    
    xa_1mpsi_dV_newiter = ...
        ( xa_1mpsi_dV_old + dt*RHSp + lambda*xa_1mpsi_dV_olditer ) ...
        /( 1 - dt*RHSm/max(xa_1mpsi_dV_olditer,1e-14) + lambda );
    
    
    onempsi_dV = xws_1mpsi_dV_newiter + xds_1mpsi_dV_newiter ...
               +  xc_1mpsi_dV_newiter +  xa_1mpsi_dV_newiter;
    
    % Update mole fractions
    x_ws(i) = xws_1mpsi_dV_newiter/onempsi_dV;
    x_ds(i) = xds_1mpsi_dV_newiter/onempsi_dV;
    x_c(i)  =  xc_1mpsi_dV_newiter/onempsi_dV;
    x_a(i)  =  xa_1mpsi_dV_newiter/onempsi_dV;
    % NB: Sum(x_k) = 1 by construction
    
    x_ws(i) = max(0,min(1,x_ws(i)));   % Safety: enforce 0 <= x_ws <= 1
    x_ds(i) = max(0,min(1,x_ds(i)));   % Safety: enforce 0 <= x_ds <= 1
    x_c(i)  = max(0,min(1,x_c(i)));    % Safety: enforce 0 <= x_c <= 1
    x_a(i)  = max(0,min(1,x_a(i)));    % Safety: enforce 0 <= x_a <= 1

    % Update porosity
    psi     = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
            + psi_c * x_c(i) + psi_a * x_a(i);
    if( (psi <= 0) | (1 <= psi) )
        fprintf(...
         ' Error in psi (mass_conservation); cell i = %g \n',i);
        return;
    end
    
    % Update cell volume
    dV(i)   = onempsi_dV/(1-psi);

end
     
end

