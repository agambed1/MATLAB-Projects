%% Function evaluatePhiJatX
function phiVal = evaluatePhiJatX(nodeJ,x,eps)

% Initialization
% Define radius r(x - x_j)
r = x - nodeJ;
phiVal = sqrt(1+(eps*r)^2);
end