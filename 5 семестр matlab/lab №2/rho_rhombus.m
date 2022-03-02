function [point, val] = rho_rhombus(x)
    a = 2;
    b = 1;
    c = [0 0];%[2 2];
    val = max(abs(x).*[b a]) + x * c' ;
    point = c;
    if b*abs(x(1)) >= a*abs(x(2))
        point(1) = point(1) + b*sign(x(1));
    else
        point(2) = point(2) + a*sign(x(2));
    end
    %disp(abs(x*point' - val) < 0.1)
end