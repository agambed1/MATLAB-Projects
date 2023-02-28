%% Function Newton
function [roots, n_iter] = Newton (x0,f,d1F,d2F,method)

% Initialization
maxIter = 20;                   % Maximum number of iterations
n = 1; % Iteration counter (count. starts from 1)
x = zeros(maxIter+1,1);         % x_n-values
x(n) = x0;                      % Initial guess x_0 (note that n=1 for x_0)
eps = 1e-12;                       % Tolerance for Newton's method
t1 = 0.001;                        % Tolerance t_1 to correct x_n
t2 = 3;                        % Tolerance t_2 to correct x_n
res = -f(x0);                   % Determine residual for initial guess x_0

% Iterate over n until |res| <= eps or n >= maxIter
while abs(res) > eps && n < maxIter+1
    if method == 2           % Correct x_n for extended Newton method
        % Iterate until one condition turns false=
        while ( f(x(n)) > t1 && (d1F(x(n)) < 0 && d1F(x(n)) > -t2 && ...
                d2F(x(n)) < 0) || ...
                (f(x(n)) > t1 && d1F(x(n)) < -t2 && d2F(x(n)) > 0) || ...
                (f(x(n)) < -t1 && (d1F(x(n)) < t2 && d1F(x(n)) > 0) && ...
                d2F(x(n)) > 0) || ...
                (f(x(n)) < -t1 && d1F(x(n)) > t2 && d2F(x(n)) < 0) )
            x(n) = x(n) + t1;   % Move x_n by t_1 to the right
        end
        % Iterate until one condition turns false
        while (f(x(n)) > t1 && (d1F(x(n)) > 0 && d1F(x(n)) < t2) && ...
                d2F(x(n)) < 0) || ...
                (f(x(n)) > t1 && d1F(x(n)) > t2 && d2F(x(n)) > 0) || ...
                (f(x(n)) < -t1 && (d1F(x(n)) < 0 && d1F(x(n)) > -t2) && ...
                d2F(x(n)) > 0) || ...
                (f(x(n)) < -t1 && d1F(x(n)) < -t2 && d2F(x(n)) < 0)
            x(n) = x(n) - t1;   % Move x_n by t_1 to the left
        end
    end

    % Newton's method
    res = -f(x(n));               % Determine residual for x_n
    DeltaX = (res/d1F(x(n)));                % Determine Delta x
    x(n+1) = x(n)+DeltaX;                % Update x_n+1 according to x_n + DeltaX
    n = n + 1;                 % Increase iteration counter
    res = -f(x(n));               % Determine residual for x_n
end

x = x(1:n);                     % Shrink x-array to n elements
roots = x(end);                  % Determine approx. root as last x_n-value
n = n - 1;                      % Decrease n (since n starts from 1, not 0)
n_iter = n;

end