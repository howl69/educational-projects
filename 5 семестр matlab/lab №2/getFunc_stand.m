function [f, g] = getFunc_stand(n)
    f = @(x, L) sin(pi * n * x ./ L);
    g = @(x, L) cos(pi * n * x ./ L);
end