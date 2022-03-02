function res = f2(x)
    res = 0;
    if abs(x) > 0.0001
        res = x.*cos(log(abs(x)));
    end
end