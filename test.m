clear all,
close all;
clc;

N = 2;          % Hamming weight of error in v(x)
q = 8;          % Alphabet size
m = log2(q);    %
d = 6;          % Minimum distance wh(c1-c2)
j0 = 1;

n = q-1;        % CW length
k = n-d+1;      % Message length
t = floor((d-1)/2); % max number of errors that can be corrected
 

% Generate random information vector
seedc = 42;
rng(seedc);
a = randi([0 q-1], 1, k); % a(x)=alpha^a(1)x^k+alpha^a(2)x^k-1+....+alpha^a(end)
% Create info poly
gf_a = gf(a,m);

% create gen poly
alpha = gf(2, m);
g = gf([1], m); %Initialise g(x)
for exp=1:2*t
    g = conv(g, gf([1 -(alpha^exp)], m)); % g(x)=(x-alpha^1)...(x-alpha^d-1)
end

a_prime = [a zeros(1,n-k)]; % shift a to the left
[quot,p] = deconv(a_prime, g);  % calculate the remainder
c = a_prime+p;  % c is the systematic codeword

% Create the error vector and add it to c

seede = 1;
rng(seede);
errorIndices = randperm(n, N); % vector of length N with random numbers between 1 and n
errorVector = zeros(1, n);
errorVector(errorIndices) = randi([1 q-1], length(errorIndices), 1);
e = gf(errorVector,m);
v = c+e;

Sj = gf(zeros(1,2*t),m);
% S_j=v(alpha^j)
for j = 1:2*t
    for exp = 1:n
% to evaluate v(alpha^j) all positions have to be summed together
% Every position is the product of the polynomials coefficient v.x(idx) and
% the position (j) to the power of n-idx
        Sj(j) = Sj(j)+gf(v.x(exp),m)*alpha^(j*(n-exp));
    end
end

% Build the S-matrix for every t,...,1
for ind = 1:-1:1
    S = gf(zeros(ind),m);
    S_vec = gf(zeros(ind,1),m);
    for l_ind = 1:ind
        S_vec(l_ind) = -Sj(ind+l_ind);
        for k_ind = 1:ind
            S(l_ind,k_ind) = Sj(l_ind+k_ind-1);
        end
    end
    % if all S_j == 0: break
    try
        Lambda = S\S_vec;  % Tries to solve equation system
        % Checks if there is a valid solution
        if ~isempty(Lambda) && all(~isnan(Lambda.x))
            break;
        end
    catch ME
        % continous the loop
    end
end


Lambda = gf([Lambda.x 1],m);    % Add coefficient for position x^0
Lambda_roots = roots(Lambda);   % Find roots
X = 1/Lambda_roots;
X_mat = gf(zeros(length(S)),m);
for idx = 1:length(S)
    X_mat(idx,:) = X.^idx;
end
Y = X\S;
X_exp = log(X)./log(alpha); % Find the exponents of alpha in the X array (X=alpha.^X_exp)
e_red = gf(zeros(1,length(v)),m);
e_red(n-X_exp) = Y;

c_red = v-e_red;

a_red = c_red.x(1:k);