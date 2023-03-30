%%  Initials
clc; close all; clear all;
%%initial of x and time variables 
m = 15; % number of observations
obs_start = 2.1; obs_end = 7; %interval of observations
time_final = 20;

time_mesh = linspace(0,time_final,m);
% x_initial = [x_T(0); x_M1(0); x_M2(0)]
% initial values that seems to fit fig 1a 
x_initial = [0; 10^3; 0]; 

% (I): (10^-6, 10^-3, 10^-3) gick bort från origo
% (II): (10^-6, 0, 0) gick mot (Beta_T, 0, 0)
% (III): Lösning 1 (+): (3*10^9, 10^(-6), 10^8) gick mot (3.1298490038×10^9, 0, 402531911.894)
%        {Beräknade i Desmos}
%        Lösning 2 (-): enligt Desmos: (160473576.808, 0,
%        −8.8025319119×10^9) inte intressant lösning pga negativ densitet
% (IV): vet inte hur vi ska göra denna!!!!

%% Steady state analysis
% We have 3 steady states with non-zero values for initial conditions on
% alpha. 
time_mesh2 = linspace(0,time_final,100);
alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh2);

%% First steady state is (0, 0, 0)
x_init = [10^(-9), 10^(-7), 10^(-5), 0;... 
          10^(-9), 10^(-7), 10^(-5), 0;...
          10^(-9), 10^(-7), 10^(-5), 0];

sch = 0;
hold on

for i=1:length(x_init)
    
    F45 = ForwardODE45(alpha,time_mesh2,x_init(:, i));
    sch = sch+1;
    subplot(2, 2, sch);
    
    plot(time_mesh2,F45(1,:),'--r',time_mesh2,F45(2,:),'--b',time_mesh2,F45(3,:),'--m','linewidth',2);
    
    title("(" + x_init(1, i) + ", " + x_init(1, i) + ", " + x_init(1, i) + ")");
    xlabel("Tid (dagar)");
    ylabel("Densitet")
    
    legend('x_T','x_{M_1}','x_{M_2}', 'location', 'northwest')
end

%% The second is (Beta_T, 0, 0)

x_init = [2*10^8, 2*10^8, 3.2*10^9, 3.2*10^9;... %Xt
          1, 10^(-7), 10^(-5), 0;... % Xm1 
          1, 10^(-7), 10^(-5), 0]; %Xm2
sch = 0;
figure
hold on
for i=1:length(x_init)
      
    time_mesh2 = linspace(0,time_final,100);
    alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh2);
    
    F45 = ForwardODE45(alpha,time_mesh2,x_init(:, i));
    sch = sch+1;
    subplot(2, 2, sch);
    
    plot(time_mesh2,F45(1,:),'--r',time_mesh2,F45(2,:),'--b',time_mesh2,F45(3,:),'--m','linewidth',2);
    
    title("(" + sprintf('%0.5g', x_init(1, i)) + ", " + sprintf('%0.5g', x_init(2, i)) + ", " + sprintf('%0.5g', x_init(3, i)) + ")");
    xlabel("Tid (dagar)");
    ylabel("Densitet")
    
    legend('x_T','x_{M_1}','x_{M_2}', 'location', 'northeast')

end

%% The third is () Are we even gonna use this?

% x_init = [10^(-9), 10^(-7), 10^(-5), 0;... 
%           10^(-9), 10^(-7), 10^(-5), 0;...
%           10^(-9), 10^(-7), 10^(-5), 0];
% sch = 0;
% for i=1:length(x_init)
%       
%     time_mesh2 = linspace(0,time_final,100);
%     alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh2);
%     
%     F45 = ForwardODE45(alpha,time_mesh2,x_init(:, i));
%     sch = sch+1;
%     subplot(2, 2, sch);
%     
%     plot(time_mesh2,F45(1,:),'--r',time_mesh2,F45(2,:),'--b',time_mesh2,F45(3,:),'--m','linewidth',2);
%     
%     
%     title("(" + x_init(1, i) + ", " + x_init(1, i) + ", " + x_init(1, i) + ")");
%     xlabel("Tid (dagar)");
%     ylabel("Densitet")
%     
%     legend('x_T','x_{M_1}','x_{M_2}', 'location', 'northwest')
% 
% end

%% Inner functions

function alpha = alpha_vec(dm1,dm2,at1,at2,k12,time_mesh)
    scaling_factor_dm1 = dm1;
    scaling_factor_dm2 = dm2;
    scaling_factor_at1 = at1;
    scaling_factor_at2 = at2;
    scaling_factor_k12 = k12;

    function_flag = 0; % constant

    exact_dm1 = ExactParameter(scaling_factor_dm1,function_flag,time_mesh); %Exact profile for dm1 to produce data.
    exact_dm2 = ExactParameter(scaling_factor_dm2,function_flag,time_mesh); %Exact profile for dm2 to produce data.
    exact_at1 = ExactParameter(scaling_factor_at1,function_flag,time_mesh); %Exact profile for at1 to produce data.
    exact_at2 = ExactParameter(scaling_factor_at2,function_flag,time_mesh); %Exact profile for at2 to produce data.
    exact_k12 = ExactParameter(scaling_factor_k12,function_flag,time_mesh); %Exact profile for k12 to produce data.

    alpha = [exact_dm1; exact_dm2; exact_at1; exact_at2; exact_k12];

end

