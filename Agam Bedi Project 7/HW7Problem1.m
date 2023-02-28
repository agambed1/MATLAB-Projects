%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HW#7 P1 Template
%   Numerical schemes to solve IVPs
%
%   Authors: Elias Gueidon and Marcus Rüter
%   Date: 07/2019 and 08/03/2022
%   Edited by: Agam Bedi
%   Date: 08/09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please fill in the missing parts as indicated by "..", or write your own
% Matlab script
%% Clear Cache
clc;
clearvars;
close all;
    
%% Initialization
h = [.5 .25 .125 .0625];                                       % Array of step sizes
X = [1.5 2 3 4 5 6 7];                                       % x-values of interest
x0 = 1;                                        % Start x-value
b = 7;                                         % End x-value
y0 = 2;                                        % Initial value y_0
set(0,'DefaultFigureRenderer','painters');      % Create vector graphics
whichMethod = 1;                               % 1 - Forward Euler method
                                                % 2 - Backward Euler method
                                                % 3 - Trapezoidal method
                                                % 4 - Heun's method
                                                % 5 - New method

results = zeros(length(h),length(X));           % Results array
dataPlot1 = zeros(length(x0:h(1):b),1);         % Plot1 array
dataPlot2 = zeros(length(x0:h(end):b),1);       % Plot2 array


% Anonymous functions
f = @(x,y) (sin(2*x)-(2*x)*y)/(x^2);                                  % f(x,y)-function
dfdy = @(x,y) -2/x;                               % Derivative df(x,y)/dy
Y = @(x) (4+cos(2)-cos(2*x))./(2*x.^2);                                    % True solution Y(x)

%% Compute Numerical Results
% Loop over h-array to compute all y_n+1 for different step sizes h
for i = 1:length(h)
    x = x0:h(i):b;                              % Nodes array
    nSI = (b-x0)/h(i);                          % Number of subintervals
    
    switch whichMethod                          % Numerical methods
        case 1
            [y,method] = computeFEulerSol(y0,f,x);
        case 2
            [y,method] = computeBEulerSol(y0,f,dfdy,x);
        case 3
            [y,method] = computeTrapezoidalSol(y0,f,dfdy,x);
        case 4
            [y,method] = computeHeunSol(y0,f,x);
        case 5
            [y,method] = computeNewMethodSol(y0,f,x);
        otherwise
            error('Unknown numerical method!');
    end
    
    % Loop over x-values of interest to store results
    for j = 1:length(X)
        id = (X(j) - x0) / (h(i)) + 1;
        results(i,j) = y(id);
    end
    
    % Store results for plots of first and last element in h-array
    if i == 1
        dataPlot1 = y;
    elseif i == length(h)
        dataPlot2 = y;
    end
end
% Compute (true) errors at x-values of interest
errors = Y(X) - results;

%% Plot Results
figure(1);
xPlot = linspace(x0,b,2000);
hold on;
grid on;
title(sprintf('HW #7, Problem 1, %s',method));
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
xlim([x0 b]);
ylim([-0.2 2]);
plot(xPlot,Y(xPlot), 'LineWidth',4);
plot(x0:h(1):b,dataPlot1,'.-', 'LineWidth',2, 'MarkerSize',18);
plot(x0:h(end):b,dataPlot2,'.-', 'LineWidth',2, 'MarkerSize',18);
legend('True solution $Y(x)$','Approx. solution $y(x)$ w/$h = 0.5$',...
    'Approx. solution $y(x)$ w/$h = 0.0625$','Location','NE',...
    'interpreter','latex');
set(gcf,'Position',[30 350 1200 750]);
set(gca,'LineWidth',2,'FontSize',20);

%% Output Results
fprintf('y-values at x-values of interest:\n');
fprintf('-----------------------------------------------------------\n');
fprintf('    h      %2d         %2d         %2d        %2d        %2d\n',X'); 
%2d\n',X');
fprintf('-----------------------------------------------------------\n');
fprintf(' %.4f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f\n',...
    [h',results]');
fprintf('\nTrue errors at x-values of interest:\n');
fprintf('----------------------------------------------------------\n');
fprintf('    h      %2d         %2d         %2d        %2d        %2d\n',X'); 
%2d\n',X');
fprintf('----------------------------------------------------------\n');
fprintf(' %.4f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f  %11.8f\n',...
    [h',errors]');
fprintf('\nRatio at which the error decreases:\n');
fprintf('-------------------------------------------------------------\n');
fprintf(' x   %.4f  -->  %.4f  -->  %.4f  -->  %.4f\n', h');
fprintf('-------------------------------------------------------------\n');
fprintf('%2d         %.5f      %.5f      %.5f\n',...
    [X',(errors(1:end-1,:) ./ errors(2:end,:))']');