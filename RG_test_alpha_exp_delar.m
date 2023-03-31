%%
clc, clf

t_min=0;t_max=20; m=500; x_initial =[5*10^6; 10^3; 10^3];%[3.1298490038*10^9; 10^-5; 402531911.894];% [160473576.808; 0;-8.8025319119*10^9];%  % Måste vi vara i SS för att antagandet om konstanta parametrar ska gälla? 
time_mesh=linspace(t_min,t_max,m);
alpha =  [10^-11 10^-12 10^-10 10^-12 10^-12]*100;  %[1,1,1,1,1];%
alpha_1=alpha_vec(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),time_mesh);
x_23s = ForwardODE23s(alpha_1,time_mesh,x_initial); 
x_Newton = ForwardNewton(alpha_1,time_mesh,x_initial);  
x_45 = ForwardODE45(alpha_1,time_mesh,x_initial);
t_plot = [1:500];


alpha_unknown=5;
%modifierad, hela
alpha_exp_23s=RG_calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton=RG_calculate_alpha_exp(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45=RG_calculate_alpha_exp(alpha,alpha_unknown,x_45,t_min,t_max);

%Orginal, hela
alpha_exp_23s_org=calculate_alpha_exp(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton_org=calculate_alpha_exp(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45_org=calculate_alpha_exp(alpha,alpha_unknown,x_45,t_min,t_max);

%modifierad, del 1
alpha_exp_23s_del1=RG_calculate_alpha_exp_del1(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton_del1=RG_calculate_alpha_exp_del1(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45_del1=RG_calculate_alpha_exp_del1(alpha,alpha_unknown,x_45,t_min,t_max);

%modifierad, del 2
alpha_exp_23s_del2=RG_calculate_alpha_exp_del2(alpha,alpha_unknown,x_23s,t_min,t_max);
alpha_exp_Newton_del2=RG_calculate_alpha_exp_del2(alpha,alpha_unknown,x_Newton,t_min,t_max);
alpha_exp_45_del2=RG_calculate_alpha_exp_del2(alpha,alpha_unknown,x_45,t_min,t_max);



figure
subplot(2,2,1)
plot(time_mesh(2:end-1),alpha_exp_23s,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_Newton,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_45,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','Explicit calculation, Newton' , 'Explicit calculation, ode45','True value of parameter')
title('Modifierad, hela')

subplot(2,2,2)
plot(time_mesh(2:end-1),alpha_exp_23s_org,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_Newton_org,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_45_org,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','Explicit calculation, Newton' , 'Explicit calculation, ode45','True value of parameter')
title('Orginal, hela')

subplot(2,2,3)
plot(time_mesh(2:end-1),alpha_exp_23s_del1,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_Newton_del1,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_45_del1,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','Explicit calculation, Newton' , 'Explicit calculation, ode45','True value of parameter')
title('Modifierad, del 1 borta')

subplot(2,2,4)
plot(time_mesh(2:end-1),alpha_exp_23s_del2,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_Newton_del2,LineWidth=1.5)
hold on
plot(time_mesh(2:end-1),alpha_exp_45_del2,LineWidth=1.5)
hold on
plot([0 20],[alpha(alpha_unknown) alpha(alpha_unknown)], 'r--')
legend('Explicit calculation, ode23s','Explicit calculation, Newton' , 'Explicit calculation, ode45','True value of parameter')
title('Orginal, del 2  borta')

hold on


%%
figure
%subplot(2,1,1)
%t = tiledlayout(1,1);

%at = axes;% first axes, save handle
%pos = get(time_mesh, 'position');% get the position vector
%pos1 = pos(2);% save the original bottom position
%pos(2) = pos(2)+pos(1); pos(4) = pos(4) -pos1; % raise bottom/reduce height->same overall upper position
%set(time_mesh, 'position', pos); % and resize first axes
%pos(2) = pos1; pos(4)=0.01; % reset bottom to original and small height
%at(2)=axes('position',pos,'color','none');  % and create the second



plot(time_mesh,x_45(1,:), 'r');
%addaxis(t_plot,x_45(1,:),'kh');
hold on
plot( time_mesh,x_45(2,:),'g');
%addaxis(t_plot,x_45(2,:),'kh');
hold on
plot(time_mesh,x_45(3,:),  'b');
%addaxis(t_plot,x_45(3,:),'kh');
hold on
legend('x_T', 'x_{M1}', 'x_{M2}');
xlabel('Mätpunkt');
ylabel('Densitet')
%time_mesh.XAxisLocation = 'top';
%time_mesh.Box = 'off';



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

end