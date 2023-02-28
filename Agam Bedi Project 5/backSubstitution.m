%% Function backSubstitution
function [x,AS,MD] = backSubstitution(beta,c,g)

% Initialization
n = length(beta);                             % Determine n from beta or g 
x = zeros(n,1);                            % Solution vector x
AS = 0;                                    % Counter for A and S
MD = 0;                                    % Counter for M and D
x(n-1:n) = g(n-1:n) ./ beta(n-1:n);            % Define x_n-1,n-values
MD = MD + 2;                               % Increase MD counter

for j = (n-2):-1:1                          % Loop over remaining rows
    x(j) = (g(j)-c(j)*x(j+2))/beta(j);                              % Compute x_j
    AS = AS + 1;                           % Increase AS counter
    MD = MD + 2;                           % Increase MD counter
end

end