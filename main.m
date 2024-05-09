close all;
clear;
clc;

%%
addpath('./functions')

%% Variables
N = 1;                          % Hamming weight of error in v(x)
q = 16;                          % Alphabet size
m = log2(q);                    %
d = 6;                          % Minimum distance wh(c1-c2)
j0 = 1;

n = q-1;                        % CW length
k = n-d+1;                      % Message length
t = floor((d-1)/2);             % max number of errors that can be corrected
alpha = gf(2, m);               % 2 = 010 = alpha

%% Encoding
[a, c] = codewort_generator(alpha, q, m, n, k, t);
[v, e] = received_cw_generator(c, n, q, m, N);

%% Decoding
[S, tm_e] = syndrome_evaluation(alpha, t, m, n, v);            % Syndrome evaluation

if tm_e     
    [Lambda, no_solution] = error_locator_polynomial(S, m, t); % Get coefficients of error locator polynomial
else
    e_red = 0;      % If there was no error while transmitting --> e(x) = 0
end

if tm_e && ~no_solution
    e_red = calc_error_poly(Lambda, S, alpha, m, n);     % Calculate error polynomial
end

% Reconstruct cw and extract information bits
c_red = v-e_red;        
a_red = c_red.x(1:k);   

%% Output

% Polynomial output
% Alpha output
% Binary coefficient output