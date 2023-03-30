function [f, df] = Forwardfunc(x,alpha)
    % alpha = (d_m1, d_m2, a_t1, a_t2, k_12)
     r = 0.93;
     beta_T = 3*10^9;
     d_m1 = alpha(1);
     d_m2 = alpha(2);
     a_t1 = alpha(3);
     a_t2 = alpha(4);
     beta_M = 9*10^8;
     sigma_m1 = 0.173;
     sigma_m2 = 0.173;
     k_12 = alpha(5);
     
     x_T = x(1); x_M1 = x(2); x_M2 = x(3);
     
     f1 = r*x_T*(1 - x_T/beta_T) - d_m1*x_M1*x_T + d_m2*x_M2*x_T;
     f2 = a_t1*x_T*x_M1*(1 - (x_M1+x_M2)/beta_M) - sigma_m1*x_M1 - k_12*x_M1*x_T;
     f3 = a_t2*x_T*x_M2*(1 - (x_M1+x_M2)/beta_M) - sigma_m2*x_M2 + k_12*x_M1*x_T;
     f = [f1; f2; f3];

     % Jacobian df
    df1_dxT = r*(1 - 2*x_T/beta_T) - d_m1*x_M1 + d_m2*x_M2;
    df1_dxM1 = -d_m1*x_T;
    df1_dxM2 = d_m2*x_T;

    df2_dxT = a_t1*x_M1*(1 - (x_M1+x_M2)/beta_M) - k_12*x_M1;
    df2_dxM1 = a_t1*x_T*(1 - (x_M1+x_M2)/beta_M) - sigma_m1 - k_12*x_T;
    df2_dxM2 = -a_t1*x_T*x_M1/beta_M;

    df3_dxT = a_t2*x_M2*(1 - (x_M1+x_M2)/beta_M) + k_12*x_M1;
    df3_dxM1 = -a_t2*x_T*x_M2/beta_M;
    df3_dxM2 = a_t2*x_T*(1 - (x_M1+x_M2)/beta_M) - sigma_m2;

    df = [df1_dxT, df1_dxM1, df1_dxM2;
          df2_dxT, df2_dxM1, df2_dxM2;
          df3_dxT, df3_dxM1, df3_dxM2];
end