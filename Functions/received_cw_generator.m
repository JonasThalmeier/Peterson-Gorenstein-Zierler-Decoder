function [v, e] = received_cw_generator(c, n, q, m, N)
    seede = 1;
    rng(seede);
    errorIndices = randperm(n, N); 
    errorVector = zeros(1, n);
    errorVector(errorIndices) = randi([1 q-1], length(errorIndices), 1);
    e = gf(errorVector,m);
    v = c+e;
end