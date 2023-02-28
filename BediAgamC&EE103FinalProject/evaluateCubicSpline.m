%% Function evaluateCubicSpline
function [yEP,interpolant] = evaluateCubicSpline(x,y,xEP)

% Initialization
interpolant = "Cubic Spline";
h = x(2) - x(1);
nSI = (x(end)-x(1))/h;
n = nSI+1;
M = zeros(n,1);                     % Vector of M_i-values
yEP = zeros(1, length(xEP));
f = zeros(n,1);                     % Right-hand side vector f
g = zeros(n-2,1);                   % Right-hand side vector g
a = h/6*ones(n-2,1);                % Vector a of lower diagonal elements
b = 2*h/3*ones(n-2,1);              % Vector b of diagonal elements
c = h/6*ones(n-3,1);                % Vector c of upper diagonal elements
alpha = zeros(n-2,1);               % Vector alpha
beta = zeros(n-2,1);                % Vector beta
nLU = n-2;                          % Define nLU for n in LU-factorization

for i = 2:n-1                       % Loop over all rows of f from 2 to n-1
    f(i) = ((y(i+1)-y(i))/(x(i+1)-x(i))) - ((y(i)-y(i-1))/(x(i)-x(i-1))); % Calculate right-hand side vector f
end

f = f(2:n-1);       %Shrink f-vector
beta(1) = b(1);     %Define beta(1)
%LU Factorization:
for i = 2:nLU                                % Loop over remaining rows
    alpha(i) = a(i)/beta(i-1);                         % Compute alpha_i
    beta(i) = b(i)-(alpha(i)*c(i-1));     % Compute beta_i
end

% Employ forward substitution to solve Lg=f:
% Define g_1-value
g(1) = f(1);

for i = 2:nLU                               % Loop over remaining rows
    g(i) = f(i)-(alpha(i)*g(i-1));        % Compute g_i
end

% Employ back substitution to solve UM=g:
M(n-1) = g(nLU)/beta(nLU);                       % Define M_n-value

for i = (nLU-1):-1:1                          % Loop over remaining rows
    M(i+1) = (g(i)-(c(i)*M(i+2)))/(beta(i));     % Compute x_i
end

% Determine natural cubic spline:
for i=2:n                              % Loop over all end points of subinter.
    % Find indices of eval. points x_EP that are inside subinterval i-1
    indices = find(xEP >= x(i-1) & xEP <= x(i));

    % Evaluate natural cubic spline at x_EP inside subinterval i-1
    yEP(indices) = ((x(i)-xEP(indices)).^3 * M(i-1) + (xEP(indices)-x(i-1)).^3*M(i))/(6*(x(i)-x(i-1)))+(((x(i)-xEP(indices))...
        *y(i-1))+((xEP(indices)-x(i-1))*y(i)))/(x(i)-x(i-1)) - (1/6)*(x(i)-x(i-1))*((x(i)-xEP(indices))*M(i-1)+(xEP(indices)-x(i-1))*M(i));
        
end

end 