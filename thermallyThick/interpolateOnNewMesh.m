% Interpolate solution on new mesh
function [temp, x_ws, x_ds, x_c, x_a, Y_O2, pres]=interpolateOnNewMesh( ...
              nx_new, xCenter, volume, dV, nx_old, xCenter_old, dV_old, ...
        temp_old, x_ws_old, x_ds_old, x_c_old, x_a_old, Y_O2_old, pres_old)

global pres_g cp_g0 MW_g R ...
       rho_ws rho_ds rho_c rho_a ...
       c_ws c_ds c_c c_a ...
       psi_ws psi_ds psi_c psi_a


% Look for q that conserves int_{x=0}^{x=volume}(q dV)
% q = x_ws, x_ds, x_c, x_a, (rhocpTemp), Y_O2, pres

% First estimate of q based on linear interpolation: q0
% q0 = x_ws0, x_ds0, x_c0, x_a0, (rhocpTemp0), Y_O200, pres00

for i=1:nx_old
    % Porosity
    psi     = psi_ws*x_ws_old(i) + psi_ds*x_ds_old(i) ...
            + psi_c*x_c_old(i) + psi_a*x_a_old(i);
    
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s   = rho_ws*c_ws*x_ws_old(i) + rho_ds*c_ds*x_ds_old(i) ...
              + rho_c*c_c*x_c_old(i) ...
              + rho_a*c_a*x_a_old(i);   % (rho*cp) in solid phase
    rhocp_g   = (101325/287/temp_old(i))*cp_g0; % (rho*cp) in gas phase
    rhocp_eff = (1-psi)*rhocp_s + psi*rhocp_g;

    rhocpTemp_old(i) = rhocp_eff*temp_old(i);
end

for i=1:nx_new
    
    ii = 0;
    for j=1:(nx_old-1)
        if( ( (xCenter_old(j)-xCenter(i))   <= 0 ) ...
          & ( (xCenter(i)-xCenter_old(j+1)) < 0  ) )
            ii = j;
        end
    end
    if( (ii == 0) & (xCenter(i) < xCenter_old(1)) )
        ii = 1;
    end
    
    if( (ii == 0) & (xCenter_old(nx_old) <= xCenter(i)) )
        ii = nx_old-1;
    end
    
    if(ii > 0) 
        x_ws0(i) = x_ws_old(ii) ...
       +(x_ws_old(ii+1)-x_ws_old(ii))*(xCenter(i)     -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));

        x_ds0(i) = x_ds_old(ii) ...
     +(x_ds_old(ii+1)-x_ds_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));

        x_c0(i) = x_c_old(ii) ...
       +(x_c_old(ii+1)-x_c_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));

        x_a0(i) = x_a_old(ii) ...
       +(x_a_old(ii+1)-x_a_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));

        rhocpTemp0(i) = rhocpTemp_old(ii) ...
                      + (rhocpTemp_old(ii+1)-rhocpTemp_old(ii)) ...
                                   *(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        
        Y_O200(i) = Y_O2_old(ii) ...
     +(Y_O2_old(ii+1)-Y_O2_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        
        pres00(i) = pres_old(ii) ...
     +(pres_old(ii+1)-pres_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
    end
    
    if(ii == 0)
        fprintf(' *** Error in interpolation routine, i = %g \n',i);
        return;
    end

end

%AT
%
x_ws = x_ws0;
x_ds = x_ds0;
x_c  = x_c0;
x_a  = x_a0;
Y_O2 = Y_O200;
pres = pres00;
for i=1:nx_new
    % Porosity
    psi = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
        + psi_c*x_c(i) + psi_a*x_a(i);
    
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s   = rho_ws*c_ws*x_ws(i) + rho_ds*c_ds*x_ds(i) ...
              + rho_c*c_c*x_c(i) ...
              + rho_a*c_a*x_a(i);   % (rho*cp) in solid phase
    rhocp_g   = pres_g*MW_g/R/temp_old(i)*cp_g0; % (rho*cp) in gas phase
    rhocp_eff = (1-psi)*rhocp_s + psi*rhocp_g;

    temp(i) = rhocpTemp0(i)/rhocp_eff;
end
return;
%}
%AT


% Find q that minimizes cost function:
%   CF = W_I*(int_q - int_q_old)^2 ...
%      + W_P*sum_{i=1}^{i=nx_new}(q(i)-qref(i))^2
%   where int_q = sum_{i=1}^{i=nx_new}(q(i)*coef(i)*dV(i))
%   and where W_I and W_P are weight coefficients
%  [Find solution of sum{j=1}^{j=nx_new} A_ij*q(j) = b_i]
%
%   d(CF)/d(q(i)) = 2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
%                 + 2*W_P*(q(i)-qref(i))
%  [d(CF)/d(q(i)) = sum{j=1}^{j=nx_new} A_ij*q(j) - b_i]
%
%   d2(CF)/d(q(i))d(q(j)) = 2*W_I*coef(i)*dV(i)*coef(j)*dV(j) ...
%                         + 2*W_P*delta_ij
%  [d2(CF)/d(q(i))d(q(j)) = A_ij]
%
% Use Method of Conjugate Gradients
%
iter_max = 10;
eps      = 1e-6;

%------ Begin q = x_ws ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_ws_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = x_ws0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_ws0;
for i=1:nx_new
    q(i) = x_ws0(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= x_ws <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' x_ws, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
    
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_ws, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_ws dV)         = %g \n',int_q);
fprintf(' integral(x_ws_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_ws, iter = %g \n'...
            ,iter);
    return;
end

q   = max(0,min(1,q));   % Safety: enforce 0 <= x_ws <= 1
x_ws = q;

%-------- End q = x_ws ------

%------ Begin q = x_ds -----
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_ds_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = x_ds0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_ds0;
for i=1:nx_new
    q(i) = x_ds0(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= x_ds <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' x_ds, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
        
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_ds, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_ds dV)         = %g \n',int_q);
fprintf(' integral(x_ds_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_ds, iter = %g \n'...
            ,iter);
    return;
end

q    = max(0,min(1,q));   % Safety: enforce 0 <= x_ds <= 1
x_ds = q;

%-------- End q = x_ds -----

%------ Begin q = x_c ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_c_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = x_c0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_c0;
for i=1:nx_new
    q(i) = x_c0(i) + coef(i)*1e-5*Q*(2*rand(1)-1);  % Add noise
    q(i) = max(0,min(1,q(i)));                      % Safety: 0 <= x_c <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' x_c, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
        
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_c, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_c dV)         = %g \n',int_q);
fprintf(' integral(x_c_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_c, iter = %g \n' ...
            ,iter);
    return;
end

q   = max(0,min(1,q));   % Safety: enforce 0 <= x_c <= 1
x_c = q;

%-------- End q = x_c ------

%------ Begin q = x_a ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_a_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = x_a0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_a0;
for i=1:nx_new
    q(i) = x_a0(i) + coef(i)*1e-5*Q*(2*rand(1)-1);  % Add noise
    q(i) = max(0,min(1,q(i)));                      % Safety: 0 <= x_a <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' x_a, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
        
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_a, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_a dV)         = %g \n',int_q);
fprintf(' integral(x_a_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_a, iter = %g \n' ...
            ,iter);
    return;
end

q   = max(0,min(1,q));   % Safety: enforce 0 <= x_a <= 1
x_a = q;

%-------- End q = x_a ------

%------ Begin q = temp -----
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = rhocpTemp_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = rhocpTemp0;

% Initial guess (taken as equal to reference solution)
%q    = rhocpTemp0;
for i=1:nx_new
    q(i) = rhocpTemp0(i) + 1e-5*Q*(2*rand(1)-1); % Add noise

    %{
    % Porosity
    psi = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
        + psi_c*x_c(i) + psi_a*x_a(i);
    
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s   = rho_ws*c_ws*x_ws(i) + rho_ds*c_ds*x_ds(i) ...
              + rho_c*c_c*x_c(i) ...
              + rho_a*c_a*x_a(i);   % (rho*cp) in solid phase
    rhocp_g   = (101325/287/temp_old(i))*cp_g0; % (rho*cp) in gas phase
    rhocp_eff = (1-psi)*rhocp_s + psi*rhocp_g;

    q(i) = max((rhocp_eff*293),q(i));   % Safety: 293 <= T
    %}

end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*dV(i) -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' rhocpTemp, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*dV(i)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
    
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*dV(i) -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' rhocpTemp, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(rhocpTemp dV)         = %g \n',int_q);
fprintf(' integral(rhocpTemp_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for temp, iter = %g \n'...
            ,iter);
    return;
end

for i=1:nx_new
    % Porosity
    psi = psi_ws*x_ws(i) + psi_ds*x_ds(i) ...
        + psi_c*x_c(i) + psi_a*x_a(i);
    
    % Product of mass density times specific heat
    %  - Porous medium treatment:
    %    rhocp_eff = (1-psi) x rhocp_s + psi x rho_g x cp_g
    rhocp_s   = rho_ws*c_ws*x_ws(i) + rho_ds*c_ds*x_ds(i) ...
              + rho_c*c_c*x_c(i) ...
              + rho_a*c_a*x_a(i);   % (rho*cp) in solid phase
    rhocp_g   = (101325/287/temp_old(i))*cp_g0; % (rho*cp) in gas phase
    rhocp_eff = (1-psi)*rhocp_s + psi*rhocp_g;

    temp(i) = q(i)/rhocp_eff;
    %temp(i) = max(293,temp(i));   % Safety: 293 <= T
end

%-------- End q = temp -----

%------ Begin q = Y_O2 ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = Y_O2_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = Y_O200;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = Y_O200;
for i=1:nx_new
    q(i) = Y_O200(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= Y_O2 <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' Y_O2, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
    
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' Y_O2, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(Y_O2 dV)         = %g \n',int_q);
fprintf(' integral(Y_O2_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for Y_O2, iter = %g \n'...
            ,iter);
    return;
end

q   = max(0,min(1,q));   % Safety: enforce 0 <= Y_O2 <= 1
Y_O2 = q;

%-------- End q = Y_O2 ------

%------ Begin q = pres ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = pres_old;

int_q_old = 0;
for i=1:nx_old
    int_q_old = int_q_old + q_old(i)*dV_old(i);
end

% Scaling of CG coefficents
W_I      = 1e+6;
W_P      = 1e-3;
Q        = int_q_old/volume;
Q        = max(Q,1e-10);
W_I      = W_I*(1*0.01*0.0005/Q/volume/dV(1));
W_P      = W_P*(1/Q);

% Reference solution
qref = pres00;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = pres00;
for i=1:nx_new
    q(i) = pres00(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= pres <= 1
end

% Calculate initial residual
iter = 0;

res = zeros(nx_new,1);
d   = zeros(nx_new,1);

int_q0 = 0;
for i=1:nx_new
    int_q0 = int_q0 + q(i)*coef(i)*dV(i);
end
for i=1:nx_new
    res(i) =-2*W_I*(int_q0 - int_q_old)*coef(i)*dV(i) ...
            -2*W_P*(q(i)-qref(i));
    d(i)   = res(i);
end

magnitude_res = 0;
for i=1:nx_new
    magnitude_res = magnitude_res + res(i)*res(i);
end
magnitude_res0 = magnitude_res;
fprintf('\n');
fprintf(' pres, magnitude_res0 = %g \n',magnitude_res0);
if(magnitude_res0 <= 1e-14)
    % Initial guess is the solution
    int_q = int_q0;
else
    
    Ad = zeros(nx_new,1);
    for iter = 1:iter_max
    
        for i=1:nx_new
            Ad(i) = (2*W_P)*d(i);
            for j=1:nx_new
                Ad(i) = Ad(i) + (2*W_I*coef(i)*dV(i)*coef(j)*dV(j))*d(j);
            end
        end
    
        alpha = 0;
        for i=1:nx_new
            alpha = alpha + d(i)*Ad(i);
        end
        alpha = magnitude_res/alpha;
    
        for i=1:nx_new
            q(i) = q(i) + alpha*d(i);
        end
    
        int_q = 0;
        for i=1:nx_new
            int_q = int_q + q(i)*coef(i)*dV(i);
        end
        for i=1:nx_new
            res(i) =-2*W_I*(int_q - int_q_old)*coef(i)*dV(i) ...
                    -2*W_P*(q(i)-qref(i));
            %res(i) = res(i) -alpha*Ad(i);
        end
    
        magnitude_res_old = magnitude_res;
        magnitude_res = 0;
        for i=1:nx_new
            magnitude_res = magnitude_res + res(i)*res(i);
        end
        beta = magnitude_res/magnitude_res_old;
    
        if( (magnitude_res/magnitude_res0) <= eps^2 )
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n',...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' pres, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(pres dV)         = %g \n',int_q);
fprintf(' integral(pres_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for pres, iter = %g \n'...
            ,iter);
    return;
end

pres = q;

%-------- End q = pres ------


 