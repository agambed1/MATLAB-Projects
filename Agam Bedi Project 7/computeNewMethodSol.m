%% computeNewMethodSol Function
function [y,method] = computeNewMethodSol(y0,f,x)
% Initialization
method = 'New Method';              % Define method for plot
h = x(2) - x(1);                    % Step size h
nSI = (x(end) - x(1))/h;           % Number of subintervals
y = zeros(nSI+1, 1);                      % y-array (approximate solution)
y(1) = y0;                          % Initial condition
% Loop over all subintervals
for n = 1:nSI
    % Compute y_n+1 according to the new method
    y(n+1) = y(n)+h*f(x(n)+h/2,y(n)+(h/2)*f(x(n),y(n)));
end
end