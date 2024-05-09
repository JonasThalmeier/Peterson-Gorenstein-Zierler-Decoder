function [S, tm_e] = syndrome_evaluation(alpha, t, m, n, v)
    % Function to generate vector S

    S = gf(zeros(2*t,1),m);
    for j = 1:2*t
        for exp = 1:n
            % to evaluate v(alpha^j) all positions have to be summed together
            % Every position is the product of the polynomials coefficient v.x(idx) and
            % the position (j) to the power of n-idx
            S(j) = S(j)+gf(v.x(exp),m)*alpha^(j*(n-exp));
        end
    end
    
    % Check if there was a tranmission error
    if sum(abs(S.x)) == 0
        tm_e = 0;
    else
        tm_e = 1;
    end
end