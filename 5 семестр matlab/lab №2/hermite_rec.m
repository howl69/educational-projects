function a=hermite_rec(n, x)
    if n < 2
        a = hermite(n, x);
    else
        a = 2*hermite_rec(n - 1, x) .* x - (n - 1) * 2 * hermite_rec(n - 2, x);
    end
end