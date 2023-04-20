function r = boundUniformRand(bound,n)
    a = bound(1);
    b = bound(2);
    r = a + (b-a).*rand(n,1);
end