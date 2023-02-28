%% Function evaluateInterpolant
function [yEP,interpolant] = evaluateInterpolant(whichInterpolant,x,y,xEP)
switch whichInterpolant                          % Numerical methods
        case 1
            [yEP,interpolant] = evaluateRBFInterpolant(x,y,xEP);
        case 2
            [yEP,interpolant] = evaluateCubicSpline(x,y,xEP);
        otherwise
            error('Unknown numerical method!');
end
end