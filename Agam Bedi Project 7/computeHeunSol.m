%% computeHeunSol Function
function [y,method] = computeHeunSol(y0,f,x)
% Initialization
method = 'Heun''s Method';          % Define method for plot
h = x(2) - x(1);                    % Step size h
nSI = (x(end) - x(1))/h;           % Number of subintervals
y = zeros(nSI, 1);                      % y-array (approximate solution)
y(1) = y0;                          % Initial condition
% Loop over all subintervals
for n = 1:nSI
    % y_n+1 based on the Forward Euler method
    yFEnp1 = y(n)+h*f(x(n),y(n));
    
    % Compute y_n+1 according to Heun's method
    y(n+1) = y(n) + (h/2)*(f(x(n),y(n))+f(x(n+1), yFEnp1));
end
end