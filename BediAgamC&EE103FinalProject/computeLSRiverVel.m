%% Function computeLSRiverVel
function [vRAtX] = computeLSRiverVel(B,polyOrd,a,x)
vRAtX = zeros(1, length(x));
for i = 1:polyOrd
    vRAtX = vRAtX + a(i)*(x/B).^i;%Construct polynomial, sum of that, coefficients of monomial multiplied by
    %the polynomial order of x's.
end