function [f, df] = ForwardODE45(alpha, time_mesh, x_initial)
    [time, f] = ode45(@(t, x) Forfunc(t, x, alpha, time_mesh), time_mesh, x_initial);
    f = f';

    % Calculate the derivatives
    df = zeros(size(f));
    for i = 1:length(time)
        [~, df_temp] = Forfunc(time(i), f(:, i), alpha, time_mesh);
        df(:, i) = df_temp(:, 1);
    end
    function [func, df_func] = Forfunc(t, x, alpha, time_mesh)
        alphaPol =@(t) [];
        for i = 1:size(alpha,1)
           Pol_i =@(t) interp1(time_mesh, alpha(i,:), t, 'linear', 'extrap');
           alphaPol =@(t) [alphaPol(t); Pol_i(t)];
        end
        alphat = real(alphaPol(t));
        [func, df_func] = Forwardfunc(x, alphat);
    end
    
end

