function myfmin(f, x0)
    n = size(x0, 2);
    if n == 2
        x = linspace(-5, 5, 1000);
        y = linspace(-5, 5, 1000);
        [X,Y] = meshgrid(x,y);
        Z = zeros(1000,1000);
        for i = 1:1000
            for j = 1:1000
                Z(i,j) = f([x(i) y(j)]);
            end
        end
        contour(X, Y, Z, 10);
        hold on
        old_x0 = x0;
    end
    for i = 1:n
        g = @(x) f([x0(1:(i-1)) x x0(i+1:n)]);
        x0(i) = fminunc(g,x0(i));
        clc
        if n == 2
            plot([old_x0(1) x0(1)], [old_x0(2) x0(2)], '--')
        end
        old_x0 = x0;
    end
    disp("my res = " + string(f(x0)))
    if n == 1
        disp("fminbnd res = " + string(f(fminbnd(f, -1000, 1000))))
    end
end