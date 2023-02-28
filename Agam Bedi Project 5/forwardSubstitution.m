%% Function forwardSubstitution
function [g,AS,MD] = forwardSubstitution(alpha,f)

% Initialization
n = length(alpha);                             % Determine n from alpha or f
g = zeros(n,1);                            % Solution vector g
AS = 0;                                    % Counter for A and S
MD = 0;                                    % Counter for M and D

g(1:2) = f(1:2);                             % Define g_1,2-values

for j = 3:n                                % Loop over remaining rows
    g(j) = f(j)-alpha(j)*g(j-2);                              % Compute g_j
    AS = AS + 1;                           % Increase AS counter
    MD = MD + 1;                           % Increase MD counter
end

end