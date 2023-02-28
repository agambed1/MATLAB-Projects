%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HW#2 P3 Template
%   Extended Newton method and (the original version of) Newton's method
%
%   Author: Marcus Rüter
%   Date: 06/30/2022
%   Edited by: Agam Bedi
%   Date: 07/08/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please fill in the missing parts as indicated by "..", or write your own
% Matlab script

%% Clear Cache
clc;
clearvars;
close all;

%% Main Script

% Initialization
n_inGuess = 1000;                 % Number of initial guesses
x = linspace(0,10, n_inGuess);               % x-values of initial guesses
roots = zeros(n_inGuess,1);            % Approximated root for each initial guess
n_iter = zeros(n_inGuess,1);           % Number of iterations for each init. guess
method = 2;                    % 1: Newton's method, 2: ext. Newton method 

% Anonymous functions
f = @(x) 2*sin(x).*cos(2*x)+0.25;                     % f(x)
d1F = @(x) 2*cos(x).*cos(2*x)-4*sin(x).*sin(2*x);                                      % f'(x)
d2F = @(x) -8*cos(x).*sin(2*x)-10*sin(x).*cos(2*x);                                      % f"(x)

for i = 1:n_inGuess                    % Loop over all initial guesses
    % Determine the approximated root and n_iter for each initial guess
    [roots(i),n_iter(i)] = Newton(x(i),f,d1F,d2F,method);
end

% Output total number of iterations
fprintf('Total number of iterations for method %d: %d\n',method,...
    sum(n_iter));

figure(1);                      % Open Figure 1
plot(x,f(x),'LineWidth',2);       % Plot the function f(x)
hold on;                        % Set hold to on
plot(x,roots,'LineWidth',2);       % Plot the roots for each initial guess
ylim([-4 12]);                  % Set the y-limits
title('HW#2, Problem 3, Figure 1','FontSize',20);   % Set title
xlabel('$x$','Interpreter','latex');                % Set x-label
ylabel('$y$','Interpreter','latex');                % Set y-label
legend('Function $y = f(x)$','Roots $y = \alpha$','Location','NW',...
    'Interpreter','latex');     % Plot legend (with LaTeX formulas)
grid on;                        % Turn on grid
set(gcf,'Position',[30 350 850 450]);   % Change position and size
set(gca,'LineWidth',2,'FontSize',20);   % Change linewidth of axes

figure(2);                      % Open Figure 2
bar(x,n_iter,'LineWidth',2);        % Bar plot of the number of iterations
title('HW#2, Problem 3, Figure 2','FontSize',20);   % Set title
xlabel('$x_0$','Interpreter','latex');              % Set x-label
ylabel('$n_{iter}$','Interpreter','latex');         % Set y-label
grid on;                        % Turn on grid
set(gcf,'Position',[30 350 850 450]);   % Change position and size
set(gca,'LineWidth',2,'FontSize',20);   % Change linewidth of axes