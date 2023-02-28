%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HW#5 P3 Template
%   LU-factorization of a pentadiagonal matrix and operation counts
%
%   Author: Marcus Rüter
%   Date: 07/20/2022
%   Edited by: Agam Bedi
%   Date: 07/27/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please fill in the missing parts as indicated by "..", or write your own
% Matlab script

%% Clear Cache
clc;
clearvars;
close all;
    
%% Main Script

% Initialization
A = [2 0 9 0 0 0; 0 -5 0 2 0 0; 6 0 4 0 -7 0; 0 -2 0 7 0 1; 0 0 1 0 -3 0; 0 0 0 5 0 8];                                   % System matrix A
f = [-2; 3; -1; 5; 9; -4];                                   % Right-hand side vector f
n = length(A);                             % Dimension of A or f
a = zeros(n,1);                            % Vector a
c = zeros(n-2,1);                           % Vector c

% Define the vectors a and c
for i = 1:(n-2)                                  % Loop from 1 to n-2
    a(i+2) = A(i+2,i);                        % Def. a (use A from row 2 on)
    c(i) = A(i,i+2);                       % Def. c (use A until row n-2)
end

b = diag(A);                               % Define b (diag. elem. of A)

% Compute alpha and beta from LU-factorization of A, and count operations
[alpha,beta,AS1,MD1] = computeLU(a,b,c);

% Compute g from Lg=f using forward substitution, and count operations
[g,AS2,MD2] = forwardSubstitution(alpha,f);

% Compute x from Ux=g using back substitution, and count operations
[x,AS3,MD3] = backSubstitution(beta,c,g);

% Use Matlab to double check if x is the correct solution
xMatlab = A \ f;

% Compute operation counts for standard Gaussian elimination
ASGE = ((2*n.^3)+(3*n.^2)-(5*n))/6;                                % Additions and subtractions
MDGE = ((n.^3)+(3*n.^2)-n)/3;                                % Multiplications and divisions
TotalGE = ASGE + MDGE;                      % Total number of operations

% Compute operation counts for LU-factorization and forward/back subst.
ASLU = AS1 + AS2 + AS3;                            % Additions and subtractions
MDLU = MD1 + MD2 + MD3;                            % Multiplications and divisions
TotalLU = ASLU + MDLU;                      % Total number of operations

% Output results
sol = sprintf('%g ', x);
fprintf('Solution vector:\nx = [ %s]^T\n', sol);
fprintf('Error (norm) of the solution: %g\n', norm(xMatlab-x));
fprintf('Operation counts for Gaussian elimination:\n');
fprintf('AS: %d, MD: %d, total: %d\n', ASGE, MDGE, TotalGE);
fprintf('Operation counts for LU-factorization:\n');
fprintf('AS: %d, MD: %d, total: %d\n', ASLU, MDLU, TotalLU);
fprintf('Gaussian elimination needs %gx more operations!\n',...
    TotalGE/TotalLU);