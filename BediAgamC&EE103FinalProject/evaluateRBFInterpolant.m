%% Function evaluateRBFInterpolant
function [yEP,interpolant] = evaluateRBFInterpolant(x,y,xEP)
% Initialization
h = x(2)-x(1);
nSI = (x(end)-x(1))/h;
nN = nSI+1;
A = zeros(nN);
eps = 1/(2*h);                      % Shape parameter (depending on h)
yEP = zeros(length(xEP), 1);
interpolant = 'Radial Basis Functions';

for i = 1:size(A, 1)                              % Loop over all rows of A
    for j = 1:size(A, 2)                          % Loop over all columns of A
        % Determine A(i,j) by evaluating phi_j(x_i)
        phiVal = evaluatePhiJatX(x(j),x(i),eps);
        A(i,j) = phiVal;
    end
end
% Solve linear system for weights w
w = A\y';
for i = 1:length(xEP)                      %phi = sum(        % Loop over all evaluation points
    for j = 1:nN                         % Loop over all rows of w
        % Construct solution y(x_EP) at all evaluation points x_EP
        phiVal = evaluatePhiJatX(x(j),xEP(i),eps);
        yEP(i) = yEP(i) + w(j)*phiVal;
    end
end
end