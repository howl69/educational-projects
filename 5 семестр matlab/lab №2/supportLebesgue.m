function [point, val] = supportLebesgue(l)
    f = @(x) x(1)^2 + x(2)^2 - 1;
    x0 = [l(1)/100, l(2)/100];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    sc_func = @(x) -(x(1)*l(1) + x(2)*l(2));
    function [c, ceq] = nlcon(x)
        c = f(x);
        ceq = [];
    end
    [point, val] = fmincon(sc_func, x0, A, b, Aeq, beq, lb, ub, @nlcon);
    val = -val;
end

