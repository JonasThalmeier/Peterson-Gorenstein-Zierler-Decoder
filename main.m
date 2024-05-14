close all;
clear;
clc;

%%
% addpath('./Functions')

%% Variables
N = 4;                          % Hamming weight of error in v(x)
q = 256;                          % Alphabet size
m = ceil(log2(q));                    %
d = 10;                          % Minimum distance wh(c1-c2)
j0 = 1;
disp_mode = 3;
seedc = 42;
seede = 42;

n = q-1;                        % CW length
k = n-d+1;                      % Message length
t = floor((d-1)/2);             % max number of errors that can be corrected
alpha = gf(2, m);               % 2 = 010 = alpha

%% Encoding
[a, c] = codewort_generator(alpha, q, m, n, k, t, seedc);
[v, e] = received_cw_generator(c, n, q, m, N, seede);

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
if e==e_red
    fprintf('\nError was recovered');
end

%% Output
fprintf('\nTX Information vector a: %s', gf_to_string(gf(a,m),alpha,m,disp_mode));
fprintf('\nTX Codeword c: %s', gf_to_string(c,alpha,m,disp_mode));
fprintf('\nError e: %s', gf_to_string(e,alpha,m,disp_mode));
fprintf('\nRX Codeword v: %s', gf_to_string(v,alpha,m,disp_mode));
fprintf('\nRecovered Error e_rec: %s', gf_to_string(e_red,alpha,m,disp_mode));
fprintf('\nRecovered Codeword c_rec: %s', gf_to_string(c_red,alpha,m,disp_mode));
fprintf('\nInformation vector a_rec: %s\n', gf_to_string(gf(a_red,m),alpha,m,disp_mode));


% Polynomial output
% Alpha output
% Binary coefficient output