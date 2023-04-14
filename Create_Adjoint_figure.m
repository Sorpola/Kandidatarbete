%%  Initials
clc; close all; clear all;
%%initial of x and time variables 
m = 15; % number of observations
obs_start = 2.1; obs_end = 7; %interval of observations
time_final = 20;

time_mesh = linspace(0,time_final,m);
% x_initial = [x_T(0); x_M1(0); x_M2(0)]
% initial values that seems to fit fig 1a 
x_initial = [5*10^6; 10^3; 10^3]; 


%% Plot forward
alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh);
noise_level = 0.1;% 10% noise

%using ode45 since newton seems to have problems with low numbers for
% this particular system of odes

[g, g_brus, g_add] = ExactODE45(alpha,time_mesh,noise_level,x_initial);
figure
hold on

m2=1000;
time_mesh2 = linspace(0,time_final,m2);
alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh2);
FN = ForwardNewton(alpha,time_mesh2,x_initial);
F45 = ForwardODE45(alpha,time_mesh2,x_initial);

subplot(1,2,1);
plot(time_mesh2,FN(1,:),'-r',time_mesh2,F45(1,:),'--b','linewidth',2);
title('Framåtproblemet');
xlabel('tid (dagar)');
legend('x_T Newton', 'x_T ode45')

%% Plot Adjoint

time_obs = linspace(obs_start, obs_end, 15);
alpha_obs = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_obs);
[g, g_brus, g_add] = ExactODE45(alpha_obs,time_obs,noise_level,x_initial); % get observations
%interpolate to be able to use for other time meshes
gPol =@(t) [];
for j = 1:size(g_brus,1)
   Pol_j =@(t) interp1(time_mesh,g_brus(j,:),t);
   gPol =@(t) [gPol(t); Pol_j(t)]; 
end

time_mesh2 = linspace(0,time_final,m2);
alpha = alpha_vec(10^-9,10^-10,10^-8,10^-10,5*10^-10,time_mesh2);
g_brus = gPol(time_mesh2);
FN = ForwardNewton(alpha,time_mesh2,x_initial);  
F45 = ForwardODE45(alpha,time_mesh2,x_initial);  
AN = AdjointNewton(alpha, FN, g_brus, time_mesh2, obs_start, obs_end);
A45 = AdjointODE45(alpha, F45, g_brus, time_mesh2, obs_start, obs_end);

subplot(1,2,2);
plot(time_mesh2,AN(2,:),'-r',time_mesh2,A45(2,:),'--b','linewidth',2)
title('Bakåtproblemet'); 
xlabel('tid (dagar)');
legend('x_T Newton', 'x_T ode45')


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