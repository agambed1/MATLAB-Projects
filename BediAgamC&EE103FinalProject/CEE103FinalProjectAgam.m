%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Final Project Structure with Code
%   A Boat Traversing a River
%
%   Author: Marcus Rüter
%   Date: 08/11/2022
%   Edited By: Agam Bedi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clc;
clearvars;
close all;

%% Initialization
whichProblem = 1;           % 1 = Lin
                            % 2 = Non-Lin
whichMethod = 1;           % 1 = Heun
                            % 2 = Backward Euler
whichInterpolant = 1;      % 1 = Radial Basis
                            % 2 = Cubic Spline
whichIntegration = 1;      % 1 = Trapezoidal Rule

color = lines(6);           % Default Matlab colors for piecewise plots
B = 431;                     % Point B [m]
delta = 1e-10;                 % delta [m] to create B+delta
vB = 5;                    % Boat velocity [m/s]
x0 = 0;                    % Start x-value [m]
y0 = 0;                    % Initial condition [m]
theta = pi/6;                 % Start value of theta
nSI = [4 16 64 256];                 % Array of numbers of subintervals
nSIMax = nSI(end);          % Maximum number of subintervals
nN = nSI + 1;                % Array of numbers of nodes
h = (B-x0) ./ nSI;            % Array of step sizes
nRun = length(nSI);         % Number of runs
y = zeros(nRun,nN(end));    % Array of y-values at nodes x_i
w = zeros(nRun,nN(end));    % Array of w-values at nodes x_j
T = zeros(nRun, 1);              % Array of travel times
nPP = nN(end)*10+1;         % Number of points for plotting
xPlot = linspace(x0,B,nPP); % x-values for plotting (1D array)
yPlot = zeros(nRun,nPP);    % y-values for plotting for all runs (2D array)
xRiver = [3.65 15.13 28.24 41.98 55.67 70.82 87.71 99.67 112.33 135.05 151.73 164.54 181.71 197.92 207.61 223.42 236.51 249.84 259.48 277.89 295.58 307.47 319.18 339.29 357.63 369.74 387.13 395.62 407.87 415.21 427.47];              % x-values of river data points
vRiver = [0.09 0.61 0.82 1.32 1.77 1.84 1.87 2.19 2.42 2.04 2.27 1.97 1.79 1.38 1.33 1.27 0.85 0.83 0.53 0.72 0.91 1.14 0.99 1.36 1.02 1.07 0.98 0.47 0.39 0.11 0.03];              % v-values of river data points
xRiverScaled = xRiver/B;    % Scaled xRiver-values for least squares
nDP = length(xRiver);           % Number of river data points
polyOrd = 7;               % Polynomial order for least squares
p = zeros(polyOrd,nDP);     % Vector of basis functions
M = zeros(polyOrd,polyOrd); % M-matrix
b = zeros(polyOrd,1);       % Right-hand side b-vector

set(0,'DefaultFigureRenderer','painters');  % Create vector graphics

%% Least-squares Data Fitting

for i = 1:polyOrd
       p(i,:) = xRiverScaled.^i; % Evaluate vector of basis functions
end
for i = 1:nDP % Loop over all river data points 
       M = M + p(:,i)*p(:,i)';
       b = b + vRiver(i)*p(:,i);
end
aLS = M\b;
vRPlot = computeLSRiverVel(B+delta, polyOrd, aLS, xPlot);% Call computeLSRiverVel to generate data points and to plot the
%least-squares function (river velocity)

%% Anonymous Functions:
% f(x) for linear IVP
fLin = @(x,theta) computeLSRiverVel(B,polyOrd,aLS,x')/(vB*cos(theta))-(tan(theta));
% f(x,y(x)) for non-linear IVP
fNonLin = @(x,y) computeLSRiverVel(B,polyOrd,aLS,x')/vB * (sqrt((B+delta-x)^2 + y^2)/(B+delta-x)) - (y/(B+delta-x));
% Derivative of f(x,y(x)) w.r.t. y for non-linear IVP
dfNonLinDy = @(x,y) computeLSRiverVel(B,polyOrd,aLS,x') .* ((y)/(vB*(B+delta-x)*(sqrt(y^2+(B+delta-x)^2))))-(1/(B+delta-x));%derivative with respect to y of line 65
% Integrand for numerical integration
f = @(x,y) sqrt((B+delta-x).^2 + (y.^2))./(vB*(B+delta-x));

%% Compute Numerical Results for Linear or Non-linear IVP
if whichProblem == 1            % Find (fixed) theta only for linear IVP
    x = x0:h(end):B; % Nodes array
    theta = pi/6;
    y(end,end) = 10;            % Large value at y(B) to start while loop
    yPlot(end,end) = 10;        % Large val. at yPlot(B) to st. while loop

    % Loop until y(B) <= 0.1 or yPlot(B) <= 0.1 (for collocation method)
    while abs(y(end, end)) >= 0.1
        % Decrease theta
        theta = theta - (1/10000);
        [y(end,:),~] = computeHeunLinSol(y0,fLin,theta,x);
    end
end

% Solve linear or non-linear IVP
for i = 1:nRun
    x = x0:h(i):B; % Nodes array

    if whichProblem == 1        % Linear IVP
        problem = "Linear Problem";         % Define problem for plot

        [y(i,1:nN(i)),method] = computeHeunLinSol(y0,fLin,theta,x);
    else                        % Non-linear IVP
        problem = "Non-linear Problem";     % Define problem for plot

        [y(i,1:nN(i)),method] = computeBEulerNonLinSol(y0,fNonLin,dfNonLinDy,x);
    end
end

%% (Numerical) Integration (and Interpolation for Non-linear IVP)
for i = 1:nRun
    x = x0:h(i):B; % Nodes array

    if whichProblem == 1                % Linear IVP
        T(i) = B/(vB*cos(theta));
        integration = "Exact Integration";
    else                                % Non-linear IVP
        xEP = x;
        [yEP, interpolant] = evaluateCubicSpline(x,y(i,1:nN(i)),xEP);

        %or
        % Initialize y_EP
        % Determine y_EP based on w-values and phi-functions
        [T(i),integration] = computeTrapezoidalRuleSol(f,xEP,yEP);
   end
end

%% Interpolation for Plotting (only for non-collocation methods)
for i = 1:nRun
    x = x0:h(i):B; % Nodes array

    [yPlot(i,:),interpolant] = evaluateInterpolant(whichInterpolant,x,y(i,1:nN(i)),xPlot);
end

%% Plot Results
figure(1);                                  % Open Figure 1

% Title for all other methods
titleText = ...
        sprintf('%s, %s, %s, %s',problem,method,interpolant,integration);

title(titleText);                           % Create title
xlabel('$x$ [m]','interpreter','latex');    % x-label
xlim([x0 B]);                               % x-limits
set(gcf,'Position',[30 350 1250 750]);      % Plot window size and position
set(gca,'LineWidth',2,'FontSize',18);       % Axes line width and font size
grid on;                                    % Turn on grid
hold on;                                    % Set hold to on

colororder({'k','k'});                      % Color order for two y-axes
yyaxis right;                               % Define right y-axis
ylabel('$v_R$ [m/s]','interpreter','latex');% y-label for right y-axis
ylim([-0.5 5]);                             % y-limits for right y-axis

% Plot LS fitting function plus horizontal line (optional) and data points
plot(xPlot,vRPlot,'-','LineWidth',2,'Color',color(1,:));
plot([x0 B],[0 0],'-','LineWidth',2,'Color',color(1,:));
plot(xRiver,vRiver,'o','MarkerSize',8,'Linewidth',2,'Color',color(2,:));

% Plot arrows to indicate direction of current (optional)
arrowsX = x0+10:(B-20)/40:B-10;             % x-positions of arrows
arrowsY = zeros(length(arrowsX),1);         % y-positions of arrows
vR = computeLSRiverVel(B,polyOrd,aLS,arrowsX);  % River velocity at x-pos.
quiver(arrowsX',vR'-0.05,arrowsY,vR'-0.05,'-^','LineWidth',2,...
    'AutoScale','off','ShowArrowHead','off','Alignment','head',...
    'MarkerSize',4,'Color',color(1,:));     % Plot arrows

yyaxis left;                                % Define left y-axis
ylabel('$y$ [m]','interpreter','latex');    % y-label for left y-axis
% y-limits for left y-axis based on maximum values of numerical results
ylim([floor(min(yPlot(2,:))/20)*20 ceil(max(yPlot(2,:))/20)*20]);

for i = 1:nRun                              % Loop over all runs
    % Plot numerical solution y(x)
    p(i) = plot(xPlot, yPlot(i,:),'-','LineWidth',3,'Color',color(i+2,:));

    % Add data points only if needed
    plot(x0:h(i):B, y(i,1:nN(i)),'o','LineWidth',3,'Color',...
        color(i+2,:),'MarkerSize',10/i,MarkerFaceColor=color(i+2,:));
end

% Create plot legend
legend([p(1) p(2) p(3) p(4)],...
    sprintf('$y_1(x)$, $y_1(B) = $ %3.2em, $T_1 = $ %5.3fs',...
    yPlot(1,end),T(1)),...
    sprintf('$y_2(x)$, $y_2(B) = $ %3.2em, $T_2 = $ %5.3fs',...
    yPlot(2,end),T(2)),...
    sprintf('$y_3(x)$, $y_3(B) = $ %3.2em, $T_3 = $ %5.3fs',...
    yPlot(3,end),T(3)),...
    sprintf('$y_4(x)$, $y_4(B) = $ %3.2em, $T_4 = $ %5.3fs',...
    yPlot(4,end),T(4)),...
    'Location','NW','interpreter','latex');