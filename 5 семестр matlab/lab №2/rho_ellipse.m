function [point, val] = rho_ellipse(x)
    a = 4;
    b = 3;
    c = [0 0];
    C = [1 0; 0 1];
    T = [a 0; 0 b];
    P = C'*T.^2*C;
    val = sqrt(x*P*x') + x*c';
    point = c + (P*x')'./(sqrt(x*P*x'));
    %disp(abs(x*point' - val) < 0.001)
end
