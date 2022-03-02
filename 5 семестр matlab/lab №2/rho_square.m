function [point, val] = rho_square(x)
    a = 2;
    c = [0 0];
    val = norm(x, 1) * a / 2 + x * c';
    point = zeros(1, 2);
    for i = 1:2
        if x(i) > 0
            point(i) = c(i) + a / 2;
        else
            point(i) = c(i) - a / 2;
        end
    end
    %disp(abs(x*point' - val) < 0.001)
end