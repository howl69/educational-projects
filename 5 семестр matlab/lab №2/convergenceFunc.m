function mov = convergenceFunc(fn,f,a,b, n, convType) 
    mov(1:n) = struct('cdata', [],'colormap', []);
    x = linspace(a,b,100);
    y_min = min(fn(x, 1));
    y_max = max(fn(x, 1));
    for i = 2:n
        y_min = min(y_min, min(fn(x, i)));
        y_max = max(y_max, max(fn(x, i)));
    end
    for i = 1:n
        plot(x, fn(x, i), x, f(x));
        xlabel('x')
        ylabel('f(x)')
        legend("fn(x)", "f(x)")
        axis([a b y_min y_max])
        if convType ~= 2
            if convType == 0
                title("среднеквадратичная " + string(Metric(convType, x, f, fn, i)));
            else
                title("равномерная " + string(Metric(convType, x, f, fn, i)));
            end
        end
        mov(i) = getframe();
    end
end