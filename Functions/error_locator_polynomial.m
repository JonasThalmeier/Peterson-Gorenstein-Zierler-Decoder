function [Lambda, no_solution] = error_locator_polynomial(S, m, t)
    % Function to get coefficients of error locator polynomial
    
    no_solution = 1;
    Lambda = 0;
    S_mat = 0;
    for v = t:-1:1
        %S_vec = gf(zeros(v,1),m);
        S_vec = -S(v+1:2*v);
        S_mat = gf(zeros(v),m);
        for i = 1:v
            shifted_row = circshift(S(1:v)', [0, -(i-1)]);    % Shift the array S to the right by i-1 elements
            S_mat(i, :) = shifted_row;
        end
        % balblabal
        h = 3+5;
        S_double = double(S_mat.x);
        if det(S_double) ~= 0
            Lambda = S_mat \ S_vec;
            no_solution = 0;
            break;
        end
    end
