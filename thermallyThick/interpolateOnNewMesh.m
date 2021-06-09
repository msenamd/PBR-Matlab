% Interpolate solution on new mesh
function [temp, x_m, x_vs, x_c] = interpolateOnNewMesh(nx_new,xCenter, ...
      volume,dV,nx_old,xCenter_old,dV_old,temp_old,x_m_old,x_vs_old,x_c_old)

global rho_m rho_vs rho_c c_m c_vs c_c T_g


% Look for q that conserves int_{x=0}^{x=volume}(q dV)
% q = x_m, x_vs, x_c, (rho_times_cp_times_temp)

% First estimate of q based on linear interpolation: q0
% q0 = x_m0, x_vs0, x_c0, (rho_times_cp_times_temp0)

for i=1:nx_old
    rho_times_cp = rho_m*c_m*x_m_old(i) + rho_vs*c_vs*x_vs_old(i) ...
                                        + rho_c*c_c*x_c_old(i);
    rho_times_cp_times_temp_old(i) = rho_times_cp*temp_old(i);
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
        x_m0(i) = x_m_old(ii) ...
       +(x_m_old(ii+1)-x_m_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        x_vs0(i) = x_vs_old(ii) ...
     +(x_vs_old(ii+1)-x_vs_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        x_c0(i) = x_c_old(ii) ...
       +(x_c_old(ii+1)-x_c_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        rho_times_cp_times_temp0(i) = rho_times_cp_times_temp_old(ii) ...
  +(rho_times_cp_times_temp_old(ii+1)-rho_times_cp_times_temp_old(ii)) ...
                                   *(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
    end
    
    if(ii == 0)
        fprintf(' *** Error in interpolation routine, i = %g \n',i);
        return;
    end

end

%AT
%{
x_m  = x_m0;
x_vs = x_vs0;
x_c  = x_c0;
for i=1:nx_new
    rho_times_cp = rho_m*c_m*x_m(i) + rho_vs*c_vs*x_vs(i) ...
                                    + rho_c*c_c*x_c(i);
    temp(i) = rho_times_cp_times_temp0(i)/rho_times_cp;
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

%------ Begin q = x_m ------
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_m_old;

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
qref = x_m0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_m0;
for i=1:nx_new
    q(i) = x_m0(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= x_m <= 1
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
fprintf(' x_m, magnitude_res0 = %g \n',magnitude_res0);
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
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_m, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_m dV)         = %g \n',int_q);
fprintf(' integral(x_m_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_m, iter = %g \n' ...
            ,iter);
    return;
end

q   = max(0,min(1,q));   % Safety: enforce 0 <= x_m <= 1
x_m = q;

%-------- End q = x_m ------

%------ Begin q = x_vs -----
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = x_vs_old;

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
qref = x_vs0;

coef = zeros(nx_new,1);
for i=1:nx_new
    if(qref(i) < 1e-10)
        coef(i) = 0;
    else
        coef(i) = 1;
    end
end

% Initial guess (taken as equal to reference solution)
%q    = x_vs0;
for i=1:nx_new
    q(i) = x_vs0(i) + coef(i)*1e-5*Q*(2*rand(1)-1); % Add noise
    q(i) = max(0,min(1,q(i)));                     % Safety: 0 <= x_vs <= 1
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
fprintf(' x_vs, magnitude_res0 = %g \n',magnitude_res0);
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
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' x_vs, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(x_vs dV)         = %g \n',int_q);
fprintf(' integral(x_vs_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for x_vs, iter = %g \n'...
            ,iter);
    return;
end

q    = max(0,min(1,q));   % Safety: enforce 0 <= x_vs <= 1
x_vs = q;

%-------- End q = x_vs -----

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
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
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

%------ Begin q = temp -----
% Calculate integral int_{x=0}^{x=volume}(q_old dV_old)
q_old = rho_times_cp_times_temp_old;

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
qref = rho_times_cp_times_temp0;

% Initial guess (taken as equal to reference solution)
%q    = rho_times_cp_times_temp0;
for i=1:nx_new
    q(i) = rho_times_cp_times_temp0(i) + 1e-5*Q*(2*rand(1)-1); % Add noise
    rho_times_cp = rho_m*c_m*x_m(i) + rho_vs*c_vs*x_vs(i) ...
                                    + rho_c*c_c*x_c(i);
    q(i) = max((rho_times_cp*T_g),q(i));   % Safety: T_g <= T
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
fprintf(' rho_times_cp_times_temp, magnitude_res0 = %g \n',magnitude_res0);
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
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
            break;
        else
            %fprintf(' iter, (magnitude_res/magnitude_res0) = %g %g \n', ...
            %          iter,(magnitude_res/magnitude_res0));
        end
    
        for i=1:nx_new
            d(i) = res(i) + beta*d(i);
        end
    
    end
    max_res = max(abs(res));
    fprintf(' rhocptemp, iter, max(abs(res)) = %g %g \n',iter,max_res);
    
end

fprintf(' integral(rho_times_cp_times_temp dV)         = %g \n',int_q);
fprintf(' integral(rho_times_cp_times_temp_old dV_old) = %g \n',int_q_old);

if( iter == iter_max)
    fprintf(' *** Error in interpolation routine for temp, iter = %g \n'...
            ,iter);
    return;
end

for i=1:nx_new
    rho_times_cp = rho_m*c_m*x_m(i) + rho_vs*c_vs*x_vs(i) ...
                                    + rho_c*c_c*x_c(i);
    temp(i) = q(i)/rho_times_cp;
    temp(i) = max(T_g,temp(i));   % Safety: T_g <= T
end

%-------- End q = temp -----


 