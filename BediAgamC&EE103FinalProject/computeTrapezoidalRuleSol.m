%% Function computeTrapezoidalRuleSol
function [TnSI,integration] = computeTrapezoidalRuleSol(f,x,y)
% Initialization
h = x(2)-x(1);                                                   
integrand = f(x,y);
integration = "Trapezoidal Rule";
% Calculate T_nSI according to the trapezoidal rule
TnSI = h/2 * (integrand(1) + 2*sum(integrand(2:end-1)) + integrand(end));
end