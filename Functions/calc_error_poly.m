function e = calc_error_poly(Lambda, S, alpha, m, n)
% Function calculates error polynomial

v = length(Lambda.x);               % Matrix dimensions of X_mat
S_vec = S(1:v);

Lambda = gf([Lambda.x' 1],m);       % Add coefficient for position x^0
Lambda_roots = roots(Lambda);       % Find roots
X = 1./Lambda_roots;
if length(X.x) ~= v
    disp('Lambda does not have enough roots');
    e = gf(zeros(1,n),m);
else
    X_mat = gf(zeros(v,v),m);
    for idx = 1:v
        X_mat(idx,:) = X.^idx;
    end

    Y = X_mat\S_vec;

    X_exp = log(X)./log(alpha); % Find the exponents of alpha in the X array (X=alpha.^X_exp)
    e = gf(zeros(1,n),m);
    e(n-X_exp) = Y;
end
end