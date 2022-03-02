function mov = fourierApprox(f,a,b, n, meth)
    mov(1:n) = struct('cdata', [],'colormap', []);
    if meth == "standart"
        L = (b - a) / 2;
        x = linspace(a, b, 1000);
        y_min = +inf;
        y_max = -inf;
        res = 1/(2*L) * trapz(x, f(x));
        for i = 1:n
            [g1, g2] = getFunc_stand(i);
            res = res + 1./L .* g1(x, L) * trapz(x, f(x).*g1(x, L)) + 1./L .* g2(x, L) .* trapz(x, f(x).*g2(x, L));
            y_min = min(y_min, min(res));
            y_max = max(y_max, max(res));
        end
        res = 1/(2*L) * trapz(x, f(x));
        for i = 1:n
            [g1, g2] = getFunc_stand(i);
            res = res + 1./L .* g1(x, L) * trapz(x, f(x).*g1(x, L)) + 1./L .* g2(x, L) .* trapz(x, f(x).*g2(x, L));
            plot(x, res, x, f(x));
            axis([a b y_min y_max])
            title(i)
            legend("standart", "f(x)")
            mov(i) = getframe();
        end
    elseif meth == "hermite"
        x = linspace(a, b, 1000);
        H_n = hermite(0, x);
        psi = (sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
        res = psi * trapz(x, f(x) .* psi);
        y_min = +inf;
        y_max = -inf;
        for i = 1:n
            H_n = hermite(i, x);
            psi = (2^i*factorial(i)*sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
            res = res + psi * trapz(x, f(x) .* psi);
            y_min = min(y_min, min(res));
            y_max = max(y_max, max(res));
        end
        H_n = hermite(0, x);
        psi = (sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
        res = psi * trapz(x, f(x) .* psi);
        for i = 1:n
            H_n = hermite(i, x);
            psi = (2^i*factorial(i)*sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
            res = res + psi * trapz(x, f(x) .* psi);
            plot(x, res, x, f(x));
            xlabel('x')
            ylabel('f(x)')
            legend("hermite", "f(x)")
            axis([a b y_min y_max])
            title(i)
            mov(i) = getframe();
        end
    elseif meth == "hermite_rec"
        x = linspace(a, b, 1000);
        H_n = hermite_rec(0, x);
        psi = (sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
        res = psi * trapz(x, f(x) .* psi);
        y_min = +inf;
        y_max = -inf;
        for i = 1:n
            H_n = hermite_rec(i, x);
            psi = (2^i*factorial(i)*sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
            res = res + psi * trapz(x, f(x) .* psi);
            y_min = min(y_min, min(res));
            y_max = max(y_max, max(res));
        end
        H_n = hermite_rec(0, x);
        psi = (sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
        res = psi * trapz(x, f(x) .* psi);
        for i = 1:n
            H_n = hermite_rec(i, x);
            psi = (2^i*factorial(i)*sqrt(pi))^(-0.5)*exp(-x.^2./2).*H_n; 
            res = res + psi * trapz(x, f(x) .* psi);
            plot(x, res, x, f(x));
            xlabel('x')
            ylabel('f(x)')
            legend("hermite_rec", "f(x)")
            axis([a b y_min y_max])
            title(i)
            mov(i) = getframe();
        end
    elseif meth == "legendre"
        x = linspace(a, b, 1000);
        res = 0;
        y_min = +inf;
        y_max = -inf;
        for i = 1:n
            P_n = legendre(i - 1, x);
            P_n = P_n(1, :);
            res = res + (2*i - 1) / 2 * P_n .* trapz(x, f(x) .* P_n);
            y_min = min(y_min, min(res));
            y_max = max(y_max, max(res));
        end
        res = 0;
        for i = 1:n
            P_n = legendre(i - 1, x);
            P_n = P_n(1, :);
            res = res + (2*i - 1) / 2 * P_n .* trapz(x, f(x) .* P_n);
            plot(x, res, x, f(x));
            xlabel('x')
            ylabel('f(x)')
            legend("legendre", "f(x)")
            axis([a b y_min y_max])
            title(i)
            mov(i) = getframe();
        end
    end
    %movie(mov);
end