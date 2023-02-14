%%  Initials
clc; close all; clear
%%initial of x and time variables
m = 100; % number of observations
obs_start = 2.1; obs_end = 7; %interval of observations
time_final = 20;

time_mesh = linspace(0,time_final,m);
% x_initial = [x_T(0); x_M1(0); x_M2(0)]
% initial values that seems to fit fig 1a
x_initial = [5*10^6; 10^3; 10^3];

%% Variations in observations
alpha = [10^-9,10^-10,10^-8,10^-10,5*10^-10];
alpha_var = 2;

if alpha_var == 1
    alpha1 = [10^-11 10^-7];
elseif alpha_var == 2
    alpha1 = [10^-12 10^-8];
elseif alpha_var == 3
    alpha1 = [10^-10 10^-6];
elseif alpha_var == 4
    alpha1 = [10^-12 10^-8];
elseif alpha_var == 5
    alpha1 = [10^-12 10^-7];
end
variation=linspace(alpha1(1),alpha1(end),5);
small_variation=2;

figure
clf
for i_big_var=1:length(variation)
    alpha(alpha_var)=variation(i_big_var);
    alpha=alpha_vec(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),time_mesh);
    F45 = ForwardODE45(alpha,time_mesh,x_initial);


    for i_small_var=1:small_variation
        line_thickness = 1/small_variation;
        if i_big_var == 1
            alpha_plus_var = alpha;
            alpha_plus_var(alpha_var,:) = alpha_plus_var(alpha_var,:)*(1+0.05*i_small_var);
            small_plus_F45 = ForwardODE45(alpha_plus_var,time_mesh,x_initial);

            subplot(3,1,1)
            plot(time_mesh,small_plus_F45(1,:),'--r','linewidth',line_thickness);
            hold on
            subplot(3,1,2)
            plot(time_mesh,small_plus_F45(2,:),'--b','linewidth',line_thickness);
            hold on
            subplot(3,1,3)
            plot(time_mesh,small_plus_F45(3,:),'--m','linewidth',line_thickness);
            hold on
        elseif i_big_var == length(variation)
            alpha_minus_var=alpha;
            alpha_minus_var(alpha_var,:) = alpha_minus_var(alpha_var,:)*(1-0.05*i_small_var);
            small_minus_F45 = ForwardODE45(alpha_minus_var,time_mesh,x_initial);

            subplot(3,1,1)
            plot(time_mesh,small_minus_F45(1,:),'--r','linewidth',line_thickness);
            hold on
            subplot(3,1,2)
            plot(time_mesh,small_minus_F45(2,:),'--b','linewidth',line_thickness);
            hold on
            subplot(3,1,3)
            plot(time_mesh,small_minus_F45(3,:),'--m','linewidth',line_thickness);
            hold on
        else
            alpha_plus_var = alpha; alpha_minus_var=alpha;
            alpha_plus_var(alpha_var,:) = alpha_plus_var(alpha_var,:)*(1+0.05*i_small_var);
            alpha_minus_var(alpha_var,:) = alpha_minus_var(alpha_var,:)*(1-0.05*i_small_var);
            small_plus_F45 = ForwardODE45(alpha_plus_var,time_mesh,x_initial);
            small_minus_F45 = ForwardODE45(alpha_minus_var,time_mesh,x_initial);

            subplot(3,1,1)
            plot(time_mesh,small_minus_F45(1,:),'--r',time_mesh,small_plus_F45(1,:),'--r','linewidth',line_thickness);
            hold on
            subplot(3,1,2)
            plot(time_mesh,small_minus_F45(2,:),'--b',time_mesh,small_plus_F45(2,:),'--b','linewidth',line_thickness);
            hold on
            subplot(3,1,3)
            plot(time_mesh,small_minus_F45(3,:),'--m',time_mesh,small_plus_F45(3,:),'--m','linewidth',line_thickness);
            hold on
        end
    end


    subplot(3,1,1)
    plot(time_mesh,F45(1,:),'r','linewidth',1.5);
    hold on
    subplot(3,1,2)
    plot(time_mesh,F45(2,:),'b','linewidth',1.5)
    hold on
    subplot(3,1,3)
    plot(time_mesh,F45(3,:),'m','linewidth',1.5)
    hold on
end
subplot(3,1,1)
xlabel('Dagar')
ylabel('Tumörstorlek')
subplot(3,1,2)
xlabel('Dagar')
ylabel('Densitet av M1 makrofager')
subplot(3,1,3)
xlabel('Dagar')
ylabel('Densitet av M2 makrofager')

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

