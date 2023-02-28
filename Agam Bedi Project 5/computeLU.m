%% Function computeLU
function [alpha, beta, AS, MD] = computeLU(a,b,c)

% Initialization
n = length(a);                             % Determine n from a or b
alpha = zeros(n,1);                         % Vector alpha
beta = zeros(n,1);                         % Vector beta
AS = 0;                                    % Counter for A and S
MD = 0;                                    % Counter for M and D

beta(1:2) = b(1:2);                          % Define beta_1,2-values

for j = 3:n                                 % Loop over remaining rows
    alpha(j) = a(j)/beta(j-2);                          % Compute alpha_j
    MD = MD + 1;                           % Increase MD counter
    
    beta(j) = b(j)-c(j-2)*alpha(j);                           % Compute beta_j
    AS = AS + 1;                           % Increase AS counter
    MD = MD + 1;                           % Increase MD counter
end

end
