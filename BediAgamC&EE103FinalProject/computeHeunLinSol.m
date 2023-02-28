%% Function computeHeunLinSol
function [y,method] = computeHeunLinSol(y0,f,theta,x)

method = 'Heun''s Method';          % Define method for plot
h = x(2) - x(1);                    % Step size h
nSI = (x(end) - x(1))/h;           % Number of subintervals
y = zeros(nSI, 1);                      % y-array (approximate solution)
y(1) = y0;                          % Initial condition
% Loop over all subintervals
for n = 1:nSI
    y(n+1) = y(n) + (h/2)*(f(x(n),theta)+f(x(n+1), theta));
end
end