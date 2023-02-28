%% computeBEulerSol Function
function [y,method] = computeBEulerSol(y0,f,dfdy,x)
% Initialization
method = 'Backward Euler Method';   % Define method for plot
h = x(2) - x(1);                    % Step size h
nSI = (x(end) - x(1))/h;           % Number of subintervals
y = zeros(nSI+1, 1);                      % y-array (approximate solution)
y(1) = y0;                         % Initial condition
eps = 1e-12;                        % Tolerance
maxIt = 150;                        % Maximum number of iterations
whichMethod = 1;                   % 1 - Newton's method
                                    % 2 - Predictor-corrector method
% Loop over all subintervals
for n = 1:nSI
    j = 1;                          % Iteration index
    ynp1 = zeros(1,maxIt);          % y_n+1^(j) initialization
    
    % Initial guess/predictor based on the Forward Euler method
    ynp1(j) = y(n)+h*f(x(n),y(n));
    
    if whichMethod == 1
        % Compute initial residual res = -F
        res = -(ynp1(j) - y(n) - h*f(x(n+1),ynp1(j)));
        
        % Loop until |res| <= eps or j >= maxIt
        while abs(res) > eps && j < maxIt
            % Compute first-order derivative of F = -res
            dFdy = 1 - h*dfdy(x(n+1),ynp1(j));
            
            % Compute delta from F' delta = res
            delta = res/dFdy;
            
            % Update ynp1, j, and residual
            ynp1(j+1) = ynp1(j) + delta;
            j = j + 1;
            res = -(ynp1(j) - y(n) - h*f(x(n+1),ynp1(j)));
        end
    else
        % Compute initial corrector
       ynp1(j+1) = y(n) + h*f(x(n+1), ynp1(j));
        
        % Loop until |y_n+1^(j+1) - y_n+1^(j)| <= eps or j >= maxIt
       while abs(ynp1(j+1)-ynp1(j)) > eps && j < maxIt
            % Update j
            j = j + 1;
            
            % Compute corrector
            ynp1(j+1) = y(n) + h*f(x(n+1), ynp1(j));
        end
    end
    % Define y_n+1 as last Newton/Predictor-corrector iteration
    y(n+1) = ynp1(j);
end
end