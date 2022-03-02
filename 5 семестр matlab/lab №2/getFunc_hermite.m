function H = getFunc_hermite(n, x)
    H = 0;%@(x) sum((-1)^j .* factorial(n) ./ (factorial(j) .* factorial(n - 2.*j)) * (2.*x).^
    %for i = 0:(floor(n / 2))
     %   H = H + (-1)^i*(factorial(n)/(factorial(i)*factorial(n-2*i)))*(2*x).^(n-2*i);
    %end
    %disp(H);
    %disp(n);
    if n == 0
        H = 1;
    elseif n > 0
        H = 2*x.*getFunc_hermite(n - 1, x) + 2*n*getFunc_hermite(n - 2, x);
    end
end