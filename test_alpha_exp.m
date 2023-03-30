%%
clc, close all

t_min=0;t_max=20; m=[250 500 1000 5000]; 
x_initial = [5*10^6; 10^3; 10^3]; 
alpha = [10^-11 10^-12 10^-10 10^-12 10^-12]*100;
alpha_unknown=1;



alpha_exp=zeros(max(m)-2,4);
for i=1:4
    m_i=m(i); h=(t_max-t_min)/(m_i-1);
    time_mesh=linspace(t_min,t_max,m_i);
    alpha_1=alpha_vec(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),time_mesh);
    
    [x, df] = ForwardODE45(alpha_1, time_mesh, x_initial);

    alpha_exp(1:m_i-2,i)=calculate_alpha_exp(alpha,alpha_unknown,x,t_min,t_max, df);
end

for i=1:4
    m_i=m(i); h=(t_max-t_min)/(m_i-1);
    time_mesh2=t_min+h:h:t_max-h;

    subplot(2,2,i)
    plot(time_mesh2,(alpha_exp(1:m_i-2,i)),LineWidth=1.5)
    hold on
    plot([t_min t_max],[(alpha(alpha_unknown)) (alpha(alpha_unknown))], 'r--')
    legend('Explicit beräknat värde','Korrekt värde')
    title("Jämförelse mellan explicit beräkning och det korrekta värdet av en parameter", m_i + " tidspunkter")
    xlabel('Dagar'); ylabel('parametern')
end

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