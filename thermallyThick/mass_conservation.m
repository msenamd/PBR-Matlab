% Mass conservation
% - Semi-implicit time integration

function [x_ws, x_ds, x_c, x_a, dV] = ...
    mass_conservation(dt, x_ws_old, x_ds_old, x_c_old, x_a_old, ...
    temp_olditer, x_ws_olditer, x_ds_olditer, x_c_olditer, x_a_olditer, ...
    Y_O2_olditer, dV_old, nx_old)

global rho_ws rho_ds rho_c rho_a ...
       psi_ws psi_ds psi_c psi_a ...
       A_R1 Ta_R1 n_R1 eta_ds_R1 ...
       A_R2 Ta_R2 n_R2 eta_c_R2 ...
       A_R3 Ta_R3 n_R3 n_O2_R3 eta_c_R3 ...
       A_R4 Ta_R4 n_R4 n_O2_R4 eta_a_R4


for i = 1:nx_old
    
    % Porosity
    psi_old = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
            + psi_c * x_c_old(i) + psi_a *x_a_old(i);
        
    % Species mass conservation statements
    % - Notations:  xws_1mpsi_dV = x_ws x (1-psi) x dV
    %               xds_1mpsi_dV = x_ds x (1-psi) x dV
    %               xc_1mpsi_dV  = x_c  x (1-psi) x dV
    %               xa_1mpsi_dV  = x_a  x (1-psi) x dV
    %               onempsi_dV   = (1-psi) x dV
    
    psi_new = psi_ws*x_ws_olditer(i) + psi_ds*x_ds_olditer(i) ...
            + psi_c * x_c_olditer(i) + psi_a * x_a_olditer(i);
    Y_O2s  = max(Y_O2_olditer(i),0);
    K_R1   = ( (rho_ws*x_ws_olditer(i)*(1-psi_new))^(n_R1-1) ) ...
             *A_R1*exp(-Ta_R1/temp_olditer(i));
    K_R2   = ( (rho_ds*x_ds_olditer(i)*(1-psi_new))^(n_R2-1) ) ...
             *A_R2*exp(-Ta_R2/temp_olditer(i));
    K_R3   = ( (rho_ds*x_ds_olditer(i)*(1-psi_new))^(n_R3-1) ) ...
             *( Y_O2s^n_O2_R3 )*A_R3*exp(-Ta_R3/temp_olditer(i));
    K_R4   = ( (rho_c *x_c_olditer(i) *(1-psi_new))^(n_R4-1) ) ...
             *( Y_O2s^n_O2_R4 )*A_R4*exp(-Ta_R4/temp_olditer(i));
    
    xws_1mpsi_dV = ( x_ws_old(i)*(1-psi_old)*dV_old(i) )*exp( -dt*K_R1 );

    xtemp2 = K_R2+K_R3-K_R1;
    if( abs( xtemp2 ) > 0 )
        A_x2         = (eta_ds_R1*rho_ws/rho_ds) ...
                       *( x_ws_old(i)*(1-psi_old)*dV_old(i) ) ...
                       *K_R1/xtemp2;
        B_x2         = ( x_ds_old(i)*(1-psi_old)*dV_old(i) ) - A_x2;
        xds_1mpsi_dV = A_x2*exp( -dt*K_R1 ) ...
                     + B_x2*exp( -dt*(K_R2+K_R3) );
        
        xtemp31 = K_R4-K_R1;
        xtemp32 = K_R4-K_R2-K_R3;
        xtemp33 = (eta_c_R2*rho_ds/rho_c)*K_R2 ...
                + (eta_c_R3*rho_ds/rho_c)*K_R3;
        xtemp4  = (eta_a_R4*rho_c/rho_a)*K_R4;
        if( ( abs( xtemp31 ) > 0 ) & ( abs( xtemp32 ) > 0 ) )
            A_x3        = A_x2*xtemp33/xtemp31;
            B_x3        = B_x2*xtemp33/xtemp32;
            C_x3        = ( x_c_old(i)*(1-psi_old)*dV_old(i) ) ...
                        - A_x3 - B_x3;
            xc_1mpsi_dV = A_x3*exp( -dt*K_R1 ) ...
                        + B_x3*exp( -dt*(K_R2+K_R3) ) ...
                        + C_x3*exp( -dt*K_R4 );
            
            xa_1mpsi_dV = 0;
            D_x4        = ( x_a_old(i)*(1-psi_old)*dV_old(i) );
            if( abs( K_R1 ) > 0 )
                A_x4        = -A_x3*xtemp4/K_R1;
                D_x4        = D_x4 - A_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + A_x4*exp( -dt*K_R1 );
            else   % Case K_R1 = 0
                A_x4        = A_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + A_x4*dt;
            end
            if( abs( (K_R2+K_R3) ) > 0 )
                B_x4        = -B_x3*xtemp4/(K_R2+K_R3);
                D_x4        = D_x4 - B_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + B_x4*exp( -dt*(K_R2+K_R3) );
            else   % Case (K_R2+K_R3) = 0
                B_x4        = B_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + B_x4*dt;
            end
            if( abs( K_R4 ) > 0 )
                C_x4        = -C_x3*xtemp4/(K_R4);
                D_x4        = D_x4 - C_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + C_x4*exp( -dt*K_R4 );
            else   % Case K_R4 = 0
                C_x4        = C_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + C_x4*dt;
            end
            xa_1mpsi_dV = xa_1mpsi_dV + D_x4;
        end
        if( ( abs( xtemp31 ) == 0 ) & ( abs( xtemp32 ) > 0 ) )
            % Case xtemp31 = K_R4-K_R1 = 0
            A_x3        = A_x2*xtemp33;
            B_x3        = B_x2*xtemp33/xtemp32;
            C_x3        = ( x_c_old(i)*(1-psi_old)*dV_old(i) ) ...
                        - B_x3;
            xc_1mpsi_dV = A_x3*exp( -dt*K_R1 ) *dt ...
                        + B_x3*exp( -dt*(K_R2+K_R3) ) ...
                        + C_x3*exp( -dt*K_R4 );
                        
            xa_1mpsi_dV = 0;
            D_x4        = ( x_a_old(i)*(1-psi_old)*dV_old(i) );
            if( abs( K_R1 ) > 0 )
                A_x4        = -A_x3*xtemp4/K_R1;
                D_x4        = D_x4 - (A_x4/K_R1);
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + A_x4*exp( -dt*K_R1 ) *(dt+(1/K_R1));
            else   % Case K_R1 = 0
                A_x4        = A_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + A_x4*0.5*(dt^2);
            end
            if( abs( (K_R2+K_R3) ) > 0 )
                B_x4        = -B_x3*xtemp4/(K_R2+K_R3);
                D_x4        = D_x4 - B_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + B_x4*exp( -dt*(K_R2+K_R3) );
            else   % Case (K_R2+K_R3) = 0
                B_x4        = B_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + B_x4*dt;
            end
            if( abs( K_R4 ) > 0 )
                C_x4        = -C_x3*xtemp4/(K_R4);
                D_x4        = D_x4 - C_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + C_x4*exp( -dt*K_R4 );
            else   % Case K_R4 = 0
                C_x4        = C_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + C_x4*dt;
            end
            xa_1mpsi_dV = xa_1mpsi_dV + D_x4;
        end
        if( ( abs( xtemp31 ) > 0 ) & ( abs( xtemp32 ) == 0 ) )
            % Case xtemp32 = K_R4-K_R2-K_R3 = 0
            A_x3        = A_x2*xtemp33/xtemp31;
            B_x3        = B_x2*xtemp33;
            C_x3        = ( x_c_old(i)*(1-psi_old)*dV_old(i) ) ...
                        - A_x3;
            xc_1mpsi_dV = A_x3*exp( -dt*K_R1 ) ...
                        + B_x3*exp( -dt*(K_R2+K_R3) ) *dt ...
                        + C_x3*exp( -dt*K_R4 );
            
            xa_1mpsi_dV = 0;
            D_x4        = ( x_a_old(i)*(1-psi_old)*dV_old(i) );
            if( abs( K_R1 ) > 0 )
                A_x4        = -A_x3*xtemp4/K_R1;
                D_x4        = D_x4 - A_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + A_x4*exp( -dt*K_R1 );
            else   % Case K_R1 = 0
                A_x4        = A_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + A_x4*dt;
            end
            if( abs( (K_R2+K_R3) ) > 0 )
                B_x4        = -B_x3*xtemp4/(K_R2+K_R3);
                D_x4        = D_x4 - (B_x4/(K_R2+K_R3));
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + B_x4*exp( -dt*(K_R2+K_R3) ) ...
                                  *(dt+(1/(K_R2+K_R3)));
            else   % Case (K_R2+K_R3) = 0
                B_x4        = B_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + B_x4*0.5*(dt^2);
            end            
            if( abs( K_R4 ) > 0 )
                C_x4        = -C_x3*xtemp4/(K_R4);
                D_x4        = D_x4 - C_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + C_x4*exp( -dt*K_R4 );
            else   % Case K_R4 = 0
                C_x4        = C_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + C_x4*dt;
            end
            xa_1mpsi_dV = xa_1mpsi_dV + D_x4;    
        end
    else   % Case xtemp2 = K_R2+K_R3-K_R1 = 0
        A_x2         = (eta_ds_R1*rho_ws/rho_ds) ...
                       *( x_ws_old(i)*(1-psi_old)*dV_old(i) ) ...
                       *K_R1;
        B_x2         = ( x_ds_old(i)*(1-psi_old)*dV_old(i) );
        xds_1mpsi_dV = A_x2*exp( -dt*K_R1 ) *dt ...
                     + B_x2*exp( -dt*K_R1 );
        
        xtemp31 = K_R4-K_R1;
        xtemp32 = K_R4-K_R2-K_R3;
        xtemp33 = (eta_c_R2*rho_ds/rho_c)*K_R2 ...
                + (eta_c_R3*rho_ds/rho_c)*K_R3;
        xtemp4  = (eta_a_R4*rho_c/rho_a)*K_R4;      
        if( abs( xtemp31 ) > 0 )
            A_x3        = A_x2*xtemp33/xtemp31;
            B_x3        = B_x2*xtemp33/xtemp31;
            C_x3        = ( x_c_old(i)*(1-psi_old)*dV_old(i) ) ...
                        + (A_x3/xtemp31) - B_x3;
            xc_1mpsi_dV = A_x3*exp( -dt*K_R1 ) *(dt-(1/xtemp31)) ...
                        + B_x3*exp( -dt*K_R1 ) ...
                        + C_x3*exp( -dt*K_R4 );
            
            xa_1mpsi_dV = 0;
            D_x4        = ( x_a_old(i)*(1-psi_old)*dV_old(i) );
            if( abs( K_R1 ) > 0 )
                A_x4        = -A_x3*xtemp4/K_R1;
                D_x4        = D_x4 - A_x4((1/K_R1)-(1/xtemp31));
                xa_1mpsi_dV = xa_1mpsi_dV ...
                         + A_x4*exp( -dt*K_R1 ) *(dt+(1/K_R1)-(1/xtemp31));
                
                B_x4        = -B_x3*xtemp4/K_R1;
                D_x4        = D_x4 - B_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + B_x4*exp( -dt*K_R1 );
            else   % Case K_R1 = 0
                A_x4        = A_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + A_x4*0.5*(dt^2);
                
                B_x4        = B_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + B_x4*dt;
            end
            if( abs( K_R4 ) > 0 )
                C_x4        = -C_x3*xtemp4/(K_R4);
                D_x4        = D_x4 - C_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + C_x4*exp( -dt*K_R4 );
            else   % Case K_R4 = 0
                C_x4        = C_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV + C_x4*dt;
            end
            xa_1mpsi_dV = xa_1mpsi_dV + D_x4;
        else   % Case xtemp31 = K_R4-K_R1 = 0
            A_x3        = 0.5*A_x2*xtemp33;
            B_x3        = B_x2*xtemp33;
            C_x3        = ( x_c_old(i)*(1-psi_old)*dV_old(i) );
            xc_1mpsi_dV = A_x3*exp( -dt*K_R1 ) *(dt^2) ...
                        + B_x3*exp( -dt*K_R1 ) *dt ...
                        + C_x3*exp( -dt*K_R1 );
            
            xa_1mpsi_dV = 0;
            D_x4        = ( x_a_old(i)*(1-psi_old)*dV_old(i) );
            if( abs( K_R1 ) > 0 )                        
                A_x4        = -A_x3*xtemp4/K_R1;
                B_x4        = -B_x3*xtemp4/K_R1;
                C_x4        = -C_x3*xtemp4/K_R1;
                D_x4        = D_x4 ...
                            - (A_x4*2/(K_R1^2)) - (B_x4/K_R1) - C_x4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + A_x4*exp( -dt*K_R1 ) ...
                                  *((dt^2)+(dt*2/K_R1)+(2/(K_R1^2))) ...
                            + B_x4*exp( -dt*K_R1 ) *(dt+(1/K_R1)) ...
                            + C_x4*exp( -dt*K_R1 );
            else   % Case K_R1 = 0
                A_x4        = A_x3*xtemp4;
                B_x4        = B_x3*xtemp4;
                C_x4        = C_x3*xtemp4;
                xa_1mpsi_dV = xa_1mpsi_dV ...
                            + A_x4*(1/3)*(dt^3) ...
                            + B_x4*0.5*(dt^2) ...
                            + C_x4*dt;
            end
            xa_1mpsi_dV = xa_1mpsi_dV + D_x4;
        end
    end
    
    onempsi_dV   = xws_1mpsi_dV + xds_1mpsi_dV ...
                 +  xc_1mpsi_dV +  xa_1mpsi_dV;
    if( (onempsi_dV <= 0) | (1 <= onempsi_dV) )
        fprintf(...
         ' Error in onempsi_dV (mass_conservation); cell i = %g \n',i);
        return;
    end
    
    % Update mole fractions
    x_ws(i) = xws_1mpsi_dV/onempsi_dV;
    x_ds(i) = xds_1mpsi_dV/onempsi_dV;
    x_c(i)  =  xc_1mpsi_dV/onempsi_dV;
    x_a(i)  =  xa_1mpsi_dV/onempsi_dV;
    % or x_a(i) = 1 - x_ws(i) - x_ds(i) - x_c(i)
    
    x_ws(i) = max(0,min(1,x_ws(i)));   % Safety: enforce 0 <= x_ws <= 1
    x_ds(i) = max(0,min(1,x_ds(i)));   % Safety: enforce 0 <= x_ds <= 1
    x_c(i)  = max(0,min(1,x_c(i)));    % Safety: enforce 0 <= x_c <= 1
    x_a(i)  = max(0,min(1,x_a(i)));    % Safety: enforce 0 <= x_a <= 1

    % Update porosity
    psi     = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
            + psi_c * x_c(i) + psi_a * x_a(i);
    
    % Update cell volume
    dV(i)   = onempsi_dV/(1-psi);

end
     
end

