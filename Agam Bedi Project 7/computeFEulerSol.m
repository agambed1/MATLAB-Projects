%% computeFEulerSol Function
function [y,method] = computeFEulerSol(y0,f,x)
% Initialization
method = 'Forward Euler Method';    % Define method for plot
h = x(2) - x(1);                    % Step size h
nSI = (x(end) - x(1))/h;           % Number of subintervals
y = zeros(nSI + 1, 1);                      % y-array (approximate solution)
y(1) = y0;                          % Initial condition
% Loop over all subintervals
for n = 1:nSI
    % Compute y_n+1 according to the Forward Euler method
    y(n+1) = y(n)+h*f(x(n),y(n));
end
end