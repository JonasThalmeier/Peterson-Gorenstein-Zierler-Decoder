function [a, c] = codewort_generator(alpha, q, m, n, k, t, seedc)

    % Generate random information vector
    rng(seedc);
    a = randi([0 q-1], 1, k);           % a(x)=alpha^a(1)x^k+alpha^a(2)x^k-1+....+alpha^a(end)
    
    % Generate generator polynomial
    g = gf([1], m);                     % Initiate g(x)
    for exp=1:2*t
        g = conv(g, gf([1 -(alpha^exp)], m));       % g(x)=(x-alpha^1)...(x-alpha^d-1)
    end
    
    a_shift = [a zeros(1,n-k)];         % shift a to the left by n-k bits
    [quot,p] = deconv(a_shift, g);      % calculate remainder
    c = a_shift+p;                      % c is the cw
end