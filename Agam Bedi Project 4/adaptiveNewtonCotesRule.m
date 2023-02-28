%% Function adaptiveTrapezoidalRule
function [x_j,n_j,I_n_j] = adaptiveNewtonCotesRule(x_i,x_ip1,eps,f,d2F,...
    d4F,NewtonCotesRule)

% Initialize x-val. to find max |f"| or |f^(4)| in subinterval [x_i,x_i+1]
xToFindMaxDeriv = linspace(x_i,x_ip1,500);

if NewtonCotesRule == 1             % Trapezoidal rule
    y = abs(d2F(xToFindMaxDeriv));                    % Determine y-array of |f"(x)|-values
    maxD2F = max(y);               % Determine max |f"(x)| in [x_i,x_i+1]
    
    % Determine the number of subsubintervals [x_j,x_j+1] (round up result)
     n_j = ceil(sqrt(maxD2F*(x_ip1-x_i).^3)/(12*eps));
elseif NewtonCotesRule == 2         % Simpson's rule
    y = abs(d4F(xToFindMaxDeriv));                    % Determine y-array of |f^(4)(x)|-val.
    maxD4F = max(y);               % Deter. max |f^(4)(x)| in [x_i,x_i+1]
    
    % Determine the number of subsubintervals [x_j,x_j+1] (round up result)
    n_j = ceil( nthroot((((maxD4F*(x_ip1-x_i).^5)/(180*eps)).^0.25), 4) );
    
    if mod(n_j,2) == 1              % Check if n_j is an odd number
        n_j = n_j + 1;             % Correct n_j to next higher even numb.
    end
else
    error('Integration rule not implemented');
end

h_j = (x_ip1-x_i)/n_j;                           % Length of each subsubint. [x_j,x_j+1]
x_j = x_i:h_j:x_ip1;                           % Evaluation points x_j in [x_i,x_i+1]
y = f(x_j);                             % Determine y-array of f(x_j)-values

if NewtonCotesRule == 1             % Trapezoidal rule
    % Determine T_n_j(f|_[x_i,x_i+1]) in the subinterval [x_i,x_i+1]
    I_n_j = h_j/2 * (y(1) + 2*sum(y(2:end-1)) + y(end));
elseif NewtonCotesRule == 2         % Simpson's rule
    % Determine S_n_j(f|_[x_i,x_i+1]) in the subinterval [x_i,x_i+1]
    I_n_j = h_j/3 * (y(1) + 4*sum(y(2:2:end-1)) + ...
        2*sum(y(3:2:end-2)) + y(end));
else
    error('Integration rule not implemented');
end