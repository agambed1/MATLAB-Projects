%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HW#4 P3 Template
%   Adaptive Newton Cotes Rules
%
%   Author: Marcus Rüter
%   Date: 07/13/2022
%   Edited by: Agam Bedi
%   Date: 7/21/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please fill in the missing parts as indicated by "..", or write your own
% Matlab script

%% Clear Cache
clc;
clearvars;
close all;

%% Main Script

% Initialization
NewtonCotesRule = 2;           % 1: Trapezoidal rule, 2: Simpson's rule
eps = 1e-7;                    % Tolerance epsilon
a = -1;                         % Lower integration limit
b = 3;                         % Upper integration limit
h_i = 0.5;                       % Length of each subinterval [x_i,x_i+1]
x_i = a:h_i:b;                       % Evaluation points x_i in interval [a,b]
n_i = (b-a)/h_i;                       % Number of subintervals [x_i,x_i+1]
n_tot = 0;                      % Total number of subsubintervals
I_n_tot = 0;                    % Value of T_n_tot(f) or S_n_tot(f)
IOfF = 3.5939652851277285;                      % True value of I(f)

% Anonymous functions
f = @(x) sin(x.^4) - ((x.^2)/10)+1;        % f(x)
d2F = @(x) ((-16*x.^6).*sin(x.^4))+((12*x.^2).*cos(x.^4))-.20;                  % f"(x)
d4F = @(x) (((256*x.^12)-(816*x.^4)).*sin(x.^4))+((24-(1152*x.^8)).*cos(x.^4));                  % f^(4)(x)

figure(1)                                   % Open Figure 1
hold on;

for i = 1:n_i                             % Loop over all subintervals
    % For each subinterval, determine its evaluation points x_j, the number
    % of subsubintervals n_j, and the value of T_n_j(f) or S_n_j(f)
    [x_j,n_j,I_n_j] = adaptiveNewtonCotesRule(x_i(i),x_i(i+1),eps,f,d2F,...
        d4F,NewtonCotesRule);
    n_tot = n_tot+n_j;                             % Update n_tot
    I_n_tot = I_n_tot+I_n_j;                           % Update I_n_tot
    xPlot = linspace(x_i(i),x_i(i+1),500);            % x-values in subint. for plot.
    plot(xPlot,f(xPlot),'b','LineWidth',2);       % Plot the function f(x)
    plot(x_j,zeros(length(x_j),1),'r.');     % Plot evaluation points x_j
end

title('HW#4, Problem 3','FontSize',20);     % Set title
xlabel('$x$','Interpreter','latex');        % Set x-label
ylabel('$y$','Interpreter','latex');        % Set y-label
legend('Function $y = f(x)$',...
    'Evaluation points $x_j$ in each $[x_i,x_{i+1}]$','Location','SW',...
    'Interpreter','latex');                 % Plot legend (using LaTeX)
grid on;                                    % Turn on grid
set(gcf,'Position',[30 350 1050 650]);      % Change position and size
set(gca,'LineWidth',2,'FontSize',20);       % Change linewidth of axes

E_n_tot = abs(IOfF-I_n_tot);                               % Definition of true error

% Output results
if NewtonCotesRule == 1                     % Trapezoidal rule
    fprintf('Total number of subsubintervals for the trapezoidal rule n_tot = %d\n', n_tot);
    fprintf('Numerical approximation T_n_tot(f) = %12.10f\n',I_n_tot);
    fprintf('Integration error |E^T_n_tot(f)| = %12.10f\n',abs(E_n_tot));
elseif NewtonCotesRule == 2                 % Simpson's rule
    fprintf('Total number of subsubintervals for Simpson''s rule n_tot = %d\n',...
        n_tot);
    fprintf('Numerical approximation S_n_tot(f) = %12.10f\n',I_n_tot);
    fprintf('Integration error |E^S_n_tot(f)| = %12.10f\n',abs(E_n_tot));
else
    error('Integration rule not implemented');
end
