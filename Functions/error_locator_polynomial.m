function [Lambda, no_solution] = error_locator_polynomial(S, m, t)
% Function to get coefficients of error locator polynomial

no_solution = 1;
Lambda = 0;
S_mat = gf(zeros(t),m);
for i = 1:t
    shifted_row = circshift(S', [0, -i+1]);
    % shifted_row = circshift(S', [0, -(i-1)]);    % Shift the array S to the right by i-1 elements
    S_mat(i, :) = shifted_row(1:t);
end
for v = t:-1:1
    %S_vec = gf(zeros(v,1),m);
    S_vec = -S(v+1:2*v);
    S_mat_small = S_mat(1:v,1:v);
    S_double = double(S_mat.x);
    if det(S_mat_small) ~= 0
        Lambda = S_mat_small \ S_vec;
        no_solution = 0;
        break;
    end
end
