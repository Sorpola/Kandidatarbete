function alpha_exp=calculate_alpha_exp(alpha,alpha_unknown,x,t_min,t_max,df)

if alpha_unknown==1 || alpha_unknown==2 || alpha_unknown==3 || alpha_unknown==4 || alpha_unknown==5

else
    disp('---Inkompatibelt värde på variablen alpha_unknown---')
end

x_T=x(1,:);
x_M1=x(2,:);
x_M2=x(3,:);
m=length(x_T);
time_step=(t_max-t_min)/m;
alpha_exp=zeros(m-2,1);

for k = 2:m-1
        x_k = x(:, k);

        dxT_dt = df(1, k);
        dxM1_dt = df(2, k);
        dxM2_dt = df(3, k);
        
        alpha_exp(k - 1) = calc_explicit(alpha, alpha_unknown, x_k, [dxT_dt, dxM1_dt, dxM2_dt]);
end
end


function alpha_exp_k=calc_explicit(alpha,alpha_unknown,x,dx_dt)
dm1=alpha(1);
dm2=alpha(2);
at1=alpha(3);
k12=alpha(5);
r=0.93;
deltaM1=0.173;
deltaM2=0.173;
betaT=3*10^9;
betaM=9*10^8;

xT=x(1);
xM1=x(2);
xM2=x(3);

dxT_dt=dx_dt(1);
dxM1_dt=dx_dt(2); 
dxM2_dt=dx_dt(3);

if alpha_unknown==1
    alpha_exp_k=(-dxT_dt+r*xT*(1-xT/betaT)+dm2*xM2*xT)/(xT*xM1);
elseif alpha_unknown==2
    alpha_exp_k=(dxT_dt-r*xT*(1-xT/betaT)+dm1*xM1*xT)/(xT*xM2);
elseif alpha_unknown==3
    alpha_exp_k=(dxM1_dt+deltaM1*xM1+k12*xM1*xT)/(xT*xM1*(1-(xM1+xM2)/betaM));
elseif alpha_unknown==4
    alpha_exp_k=(dxM2_dt+deltaM2*xM2-k12*xM1*xT)/(xT*xM2*(1-(xM1+xM2)/betaM));
else
    alpha_exp_k=(-dxM1_dt+at1*xT*xM1*(1-(xM1+xM2)/betaM)-deltaM1*xM1)/(xT*xM1);
end
end