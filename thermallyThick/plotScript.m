clc;
clear;
close all;

global A_R1 Ta_R1 A_R2 Ta_R2 eta_c A_R3 Ta_R3 rho_m rho_vs rho_c ...
       k_m k_vs k_c c_m c_vs c_c DeltaH_R1 DeltaH_R2 DeltaH_R3 x_O2_g ...
       T_g G eps_m eps_vs eps_c dx_i geometry 

% Read input parameters
[T_g, u_g, G, x_O2_g, delta_i, FMC, ...
 rho_m, rho_vs, rho_c, k_m, k_vs, k_c, c_m, c_vs, c_c, ...
 eps_m, eps_vs, eps_c, A_R1, Ta_R1, DeltaH_R1, ...
 A_R2, Ta_R2, DeltaH_R2, eta_c, A_R3, Ta_R3, DeltaH_R3,...
 temp_surf_i,T_end,geometry] = input_parameters;

% Read model results
num_cells=load("num_cells.csv");
time_model=load("time.csv");
cell_center=load("cell_center.csv");
T_model=load("T.csv");
x_m_model=load("x_m.csv");
x_vs_model=load("x_vs.csv");
x_c_model=load("x_c.csv");
delta_model=load("delta.csv");   %radius or length
MLR=load("MLR.csv");   %radius or length

n_output=length(time_model);

% Comparisons with analytical solutions
% comment or un-comment required secitions to plot required case

%% >>>>>>>>>>>>> START Genral Plots <<<<<<<<<<<<< 
figure(1);   % Temperature vs x at different times
m=1;
for n=round([0.3 0.6 0.9 1]*n_output)
    time_plot(m)=time_model(n);
    plot(cell_center(1:num_cells(n),n),T_model(1:num_cells(n),n),'linewidth',1.2);
    hold on
    legend_all(m)={[sprintf('%.0f',time_plot(m)),' (s)']};
    m=m+1;
end
title('Spacial Temperature Distribution','fontsize',12,'Interpreter','latex');
xlabel('r (m)','fontsize',12,'Interpreter','latex');
%xlim([0.0 0.01]);
ylabel('T (K)','fontsize',12,'Interpreter','latex');
%ylim([0 90]);
box on
set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.012 0.012],'ticklabelinterpreter','latex');
legend(legend_all,'fontsize',12,'Interpreter','latex');  

figure(2);   % particle length scale vs time
plot(time_model(1:n_output),delta_model(1:n_output),'-b','linewidth',1.2);
title('Particle Length(Radius) Variation with Time','fontsize',12,'Interpreter','latex');
hold on
xlabel('time (s)','fontsize',12,'Interpreter','latex');
%xlim([0.0 0.025]);
ylabel('sample radius or length (m)','fontsize',12,'Interpreter','latex');
%ylim([0.0 0.014]);
%legend('Model','Anlytical','fontsize',10,'Interpreter','latex');
box on
set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.01 0.01],'ticklabelinterpreter','latex');

figure(3)  %volume frations vs time
AA=plot(time_model(1:n_output),x_m_model(1,(1:n_output)),'-r','linewidth',1.2);
hold on;
BB=plot(time_model(1:n_output),x_vs_model(1,(1:n_output)),'-k','linewidth',1.2);
CC=plot(time_model(1:n_output),x_c_model(1,(1:n_output)),'-b','linewidth',1.2);
title('Centerline Volume Fraction Variation with Time','fontsize',12,'Interpreter','latex');
xlabel('time (s)','fontsize',12,'Interpreter','latex');
%xlim([0 3]);
ylabel('volume fraction','fontsize',12,'Interpreter','latex');
ylim([0 1.1]);
legend([AA BB CC],'moisture','virgin solid','char'...
    ,'fontsize',12,'Interpreter','latex');
box on
set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.01 0.01],'ticklabelinterpreter','latex');
hold off

figure(4);   % particle length scale vs time
plot(time_model(1:n_output),MLR(1:n_output),'-b','linewidth',1.2);
title('Mass Loss Rate','fontsize',12,'Interpreter','latex');
hold on
xlabel('time (s)','fontsize',12,'Interpreter','latex');
%xlim([0.0 0.025]);
ylabel('Mass Loss Rate','fontsize',12,'Interpreter','latex');
%ylim([0.0 0.014]);
%legend('Model','Anlytical','fontsize',10,'Interpreter','latex');
box on
set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.01 0.01],'ticklabelinterpreter','latex');


 
%% >>>>>>>>>>>>> START comparison with analytical solution of <<<<<<<<<<<<< 
% %% >>>>>>>>>>>>> heat transfer in rectangle, sphere, cylinder <<<<<<<<<<<<<
% figure(1);   % Temperature
% x_th=linspace(0,delta_i,20);
% m=1;
% for n=round([0.25 0.5 0.75 1]*n_output)
%     curve1= plot(cell_center(1:num_cells(n),n),T_model(1:num_cells(n),n),'-r','linewidth',1.2);
%     hold on   
%    
% time_th(m)=time_model(n);
% tau(m)=(k_c/rho_vs/c_vs) * time_th(m)/delta_i^2;
% if geometry=="rectangle"
%     lambda1=1.0769;
%     A1=1.1785;
%     for i=1:length(x_th)
%     T_th(i,m) = 300 +(900 - 300)*A1*exp(-lambda1^2*tau(m))*cos(lambda1*x_th(i)/delta_i);
%     end
% elseif geometry=="cylinder"
%     lambda1=1.5995;
%     A1=1.3384;
%     for i=1:length(x_th)
%     T_th(i,m) = 300 +(900 - 300)*A1*exp(-lambda1^2*tau(m))*besselj(0,lambda1*x_th(i)/delta_i);
%     end
% elseif geometry=="sphere"
%     lambda1=2.0288;
%     A1=1.4793;
%     for i=1:length(x_th)
%     T_th(i,m) = 300 +(900 - 300)*A1*exp(-lambda1^2*tau(m))*sin(lambda1*x_th(i)/delta_i)/(lambda1*x_th(i)/delta_i);
%     end
% end
% 
% figure(1)
% curve2=plot(x_th,T_th(:,m),'ok','linewidth',1,'MarkerSize',4);
% m=m+1;
% end
% 
% figure(1)
% xlabel('r (m)','fontsize',12,'Interpreter','latex');
% xlim([0.0 0.01]);
% ylabel('T (K)','fontsize',12,'Interpreter','latex');
% %ylim([0 90]);
% legend([curve1 curve2],'Model','Analytical','fontsize',12,'Interpreter','latex');
% box on
% set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.012 0.012],'ticklabelinterpreter','latex');
% %% >>>>>>>>>>>>>  END comparison with analytical solution of  <<<<<<<<<<<<< 
% %% >>>>>>>>>>>>> heat transfer in rectangle, sphere, cylinder <<<<<<<<<<<<<






%% >>>>>>>>>>>>> START comparison with analytical solution of <<<<<<<<<<<<< 
% %% >>>>>>>>>>>>>    heat transfer in semi-infinte solid       <<<<<<<<<<<<<
% figure(1);   % Temperature
% n_output=length(time_model);
% x_th=linspace(0,0.03,100);
% for n=round([0.1 0.3 0.6 0.9]*n_output)
%     AA= plot(flipud(cell_center(:,n)),T_model(:,n),'-r','linewidth',1.2);
%     hold on   
%     
% t(n)=time_model(n);
% X(n)= k_vs/(rho_vs*c_vs) * t(n);
% for i=1:length(x_th)
%     kk(i,n) = 1- erf(x_th(i)/2/sqrt(X(n))) - exp(h*x_th(i)/k_vs+ h^2*X(n)/k_vs^2)  * (1-erf(x_th(i)/2/sqrt(X(n))+h*sqrt(X(n))/k_vs));    
%     T_th(i,n) = 300 +(T_g - 300)*kk(i,n);
% end
% figure(1)
% BB=plot(x_th,T_th(:,n),'ok','linewidth',1,'MarkerSize',3);
% end
% 
% figure(1)
% xlabel('x (m)','fontsize',12,'Interpreter','latex');
% xlim([0.0 0.025]);
% ylabel('T (K)','fontsize',12,'Interpreter','latex');
% %ylim([0 90]);
% legend([AA BB],'Model','Analytical','fontsize',12,'Interpreter','latex');
% box on
% set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.012 0.012],'ticklabelinterpreter','latex');
% % >>>>>>>>>>>>> END comparison with analytical solution of <<<<<<<<<<<<< 
% % >>>>>>>>>>>>>    heat transfer in semi-infinte solid       <<<<<<<<<<<<<


%% >>>>>>>>>>>>> START comparison with analytical solution of <<<<<<<<<<<<< 
% %% >>>>>>>>>>>>>        species mass conservation             <<<<<<<<<<<<<
%  
% %analyrtical solution
%  t=linspace(0,time_model(end),1000);
%  x_m_i=x_m_model(1);
%  x_vs_i=x_vs_model(1);
%  x_c_i=x_c_model(1);
% 
%  K_R1=A_R1*exp(-Ta_R1/T_g);
%  K_R2=A_R2*exp(-Ta_R2/T_g);
%  K_R3=A_R3*exp(-Ta_R3/T_g);
%  if geometry=="rectangle"
%      V_i=delta_i;             %rectangle
%  elseif geometry=="cylinder"
%      V_i=pi*delta_i^2;        %cylinder
%  elseif geometry=="sphere"
%      V_i=4/3*pi*delta_i^3;    %sphere
%  end
%  alpha=x_O2_g*K_R3;
%  beta=(-K_R1+x_O2_g*K_R3)*x_m_i*V_i;
%  gamma=(x_O2_g*K_R3-(1-eta_c*rho_vs/rho_c)*K_R2)*x_vs_i*V_i;
%  C=V_i-beta/(alpha-K_R1)-gamma/(alpha-K_R2);
%  for i=1:length(t)
%      V_analytical(i)=(beta/(alpha-K_R1)*exp((alpha-K_R1)*t(i)) + gamma/(alpha-K_R2)*exp((alpha-K_R2)*t(i)) + C) ...
%          / exp(alpha*t(i));
%      
%      if geometry=="rectangle"
%      delta_analytical(i)=V_analytical(i);                 %rectangle length
%      elseif geometry=="cylinder"
%      delta_analytical(i)=sqrt(V_analytical(i)/pi);        %cylinder radius
%      elseif geometry=="sphere"
%      delta_analytical(i)=(3*V_analytical(i)/4/pi)^(1/3);  %sphere radius
%      end
%      
%      x_vs_analytical(i)=V_i*x_vs_i*exp(-K_R2*t(i)) / V_analytical(i) ;
%      x_m_analytical(i)=V_i*x_m_i*exp(-K_R1*t(i)) / V_analytical(i) ;
%      x_c_analytical(i)=1 -x_m_analytical(i)-x_vs_analytical(i) ;
% 
%  end 
% 
% figure(2) %thickness vs time
% plot(time_model,delta_model,'-b','linewidth',1.2);
% hold on
% p=plot(t,delta_analytical,'ob','linewidth',0.7,'MarkerSize',5);
% p.MarkerIndices = 1:30:length(t);
% xlabel('time (s)','fontsize',12,'Interpreter','latex');
% %xlim([0.0 0.025]);
% ylabel('sample radius or length (m)','fontsize',12,'Interpreter','latex');
% %ylim([0.0 0.014]);
% legend('Model','Anlytical','fontsize',10,'Interpreter','latex');
% box on
% set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.01 0.01],'ticklabelinterpreter','latex');
% 
% 
% figure(3)  %volume frations vs time
% AA=plot(time_model,x_m_model(1,:),'-r','linewidth',1.2);
% hold on;
% BB=plot(time_model,x_vs_model(1,:),'-k','linewidth',1.2);
% CC=plot(time_model,x_c_model(1,:),'-b','linewidth',1.2);
% 
% DD=plot(t,x_m_analytical,'ok','linewidth',1,'MarkerSize',4);
% EE=plot(t,x_vs_analytical,'sqk','linewidth',1,'MarkerSize',5);
% FF=plot(t,x_c_analytical,'xk','linewidth',1,'MarkerSize',4);
% DD.MarkerIndices = round(linspace(1,length(t),150));
% EE.MarkerIndices = round(linspace(1,length(t),150));
% FF.MarkerIndices = round(linspace(1,length(t),150));
% xlabel('time (s)','fontsize',12,'Interpreter','latex');
% %xlim([0 3]);
% ylabel('volume fraction','fontsize',12,'Interpreter','latex');
% ylim([0 1.1]);
% legend([AA BB CC DD EE FF],'moisture','virgin solid','char'...
%     ,'moisture (Analytical)','virgin solid (Analytical)','char (Analytical)','fontsize',10,'Interpreter','latex');
% box on
% set(gca,'fontsize',12,'linewidth',1,'TickLength',[0.012 0.012],'ticklabelinterpreter','latex');
% % >>>>>>>>>>>>> END comparison with analytical solution of <<<<<<<<<<<<< 
% % >>>>>>>>>>>>>        species mass conservation             <<<<<<<<<<<<<
% 
% 




return; 