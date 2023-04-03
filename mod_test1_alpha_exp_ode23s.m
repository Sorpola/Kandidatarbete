%%
clc, clf

t_min=0;t_max=20; m=500; x_initial =[5*10^6; 10^3; 10^3];%[3.1298490038*10^9; 10^-5; 402531911.894];% [160473576.808; 0;-8.8025319119*10^9];%  % Måste vi vara i SS för att antagandet om konstanta parametrar ska gälla? 
time_mesh=linspace(t_min,t_max,m);
alpha =  [10^-11 10^-12 10^-10 10^-12 10^-12]*100;  %[1,1,1,1,1];%
alpha_1=alpha_vec(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),time_mesh);
x_23s = ForwardODE23s(alpha_1,time_mesh,x_initial); 
t_plot = [1:500];

gron = [102,194,165]/255;
orange = [252,141,98]/255;
lila = [141,160,203]/255;

alpha_unknown=4;
alpha_exp_23s=mod_calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);

%orginal
alpha_exp_23s_org=calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);


figure('name','Vanlig')
subplot(2,1,1)
plot(time_mesh(2:end-1),alpha_exp_23s,'color', gron,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','True value of parameter')
title('Modifierad')


subplot(2,1,2)
plot(time_mesh(2:end-1),alpha_exp_23s_org, 'color',gron,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','True value of parameter')
title('Orginal')



%%
figure


plot(time_mesh,x_45(1,:), 'r');
hold on
plot( time_mesh,x_45(2,:),'g');
hold on
plot(time_mesh,x_45(3,:),  'b');
hold on
legend('x_T', 'x_{M1}', 'x_{M2}');
xlabel('Mätpunkt');
ylabel('Densitet')




%% 
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
