%% 1
clc
a = -10;
b = 10;
n = 1000;
x = linspace(a, b, n);
y = f(x);
[ymax, max_n] = max(y);
xmax = x(max_n);
[ymin, min_n] = min(y);
xmin = x(min_n);
plot(x, y, xmax, ymax, 'ro', xmin, ymin, 'ro')
title('f(x) = |x|cos(2x^2 - 1)');
xlabel('x');
ylabel('f(x)');
text(xmax, ymax,'  max');
text(xmin, ymin,'  min');
disp(max(y));
disp(min(y));
%% 2
n = input("Enter n:\n");
flag = 0;
if n > 1 && n ~= inf
    k = 2;
    while (k * k <= n)
        if mod(n, k) == 0
            break
        end
        k = k + 1;
    end
    if mod(n, k)
        flag = 1;
    end
end
if flag == 0
    disp("Not prime");
else
    v = 7:14:n;
    mat = repmat(2:(n + 1), n, 1);
    B = ones(n + 1).*reshape((1:(n + 1)^2), [n + 1,n+1])';
    c = reshape(B, 1, (n + 1)^2);
    D = B(1:(n + 1), (n:(n + 1)));
end
%% 3
A = 5 + randn(13,7)*sqrt(0.001);
diagMax = max(abs(diag(A)));
maxRatio = max(prod(A) ./ sum(A));
minRatio = min(prod(A) ./ sum(A));
sortedA = sortrows(A, 'descend');
%% 4.1
n = 7;
A = zeros(2*n + 1);
A(n + 1, n + 1) = 2;
i = [n n+2];
A(i, i) = -1;
A([1, 2*n + 1],2:2:(2*n + 1)) = 1;
A(2:2:(2*n + 1), [1, 2*n + 1]) = 1;
%% 4.2
n = 7;
A = zeros(2*n + 1);
A(n + 1, n + 1) = 2;
for i = 1:4
    A(n, n + 2) = -1;
    A(1, 2:2:(2*n + 1)) = 1;
    A = rot90(A);
end
%% 4.3
n = 5;
L = diag([0 0 ones(1, 2*n - 1)]);
L(2:2:(2*n+1), 1) = 1;
L([1 2*n+1], 2) = 1;
L(n + 2, n) = 1;
U = zeros(2*n + 1);
U(1, [1 2*n+1]) = 1;
U(2, (2:2:(2*n+1))) = 1;
U(n, [n n+2]) = 2;
U(n+1, n+1) = -1;
A = L*U;
%% 5
a1 = [0 2; 4 5; 6 7];
a2 = a1(:, 2:-1:1)';
a2(1,:) = -a2(1,:);
A = a1*a2;
%% 6
a = [1 2 3];
b = [2 5 2];
res = max(max(a) - min(b), max(b) - min(a));
%% 7
A = [0 5 0; 1 1 0; 1 1 0; 0 2 1];
[row, col, v] = find(A);
[C, ia, ic] = unique(col);
disp(row(ia));
%% 8
n = 4;
k = 3;
A = [1 2 0; 4 3 1; 2 2 1; 2 0 1];
X_2 = repmat(sum(A.^2, 2), 1, n);
Y_2 = X_2';
XY = A*A';
R1 = sqrt(X_2 + Y_2 - 2*XY);
disp(R1)
%% 8 проверка
n = 4;
k = 3;
A = [1 2 0; 4 3 1; 2 2 1; 2 0 1];
B = (repmat(A, 1, n) - repmat(reshape(A', 1, n*k), n, 1)).^2;
for n = 1:n
    R(:,n) = sqrt(sum(B(:,(k*(n-1)+1):(k*n)), 2));
end
disp(R)
%% 9
numTries = 100;
vectorSize = 1000;
res = zeros(numTries,1); 
res1_median = zeros(vectorSize, 1);
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n);
        y = rand(1, n);
        tic(); 
        my_prod(x,y);
        cnt = 0;
        for j = 1:min(length(x), length(y))
            cnt = cnt + x(j)*y(j);
        end
        res(i) = toc();
    end
    res1_median(n) = median(res);
end
res2_median = zeros(vectorSize, 1);
res = zeros(numTries,1); 
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n);
        y = rand(1, n);
        tic(); 
        dot(x,y);
        res(i) = toc();
    end
    res2_median(n) = median(res);
end
res3_median = zeros(vectorSize, 1);
res = zeros(numTries,1); 
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n); 
        y = rand(1, n);
        tic(); 
        x*y';
        res(i) = toc();
    end
    res3_median(n) = median(res);
end
x = 1:vectorSize;
plot(x, res1_median, x, res2_median, x, res3_median)
legend('my prod','dot','x*y')
title('operating time comparison')
ylabel('time');
xlabel('array size');
%% 9
numTries = 100;
vectorSize = 1000;
res = zeros(numTries, vectorSize); 
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n);
        y = rand(1, n);
        tic(); 
        my_prod(x,y);
        cnt = 0;
        for j = 1:min(length(x), length(y))
            cnt = cnt + x(j)*y(j);
        end
        res(i,n) = toc();
    end
end
res1_median = median(res);
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n);
        y = rand(1, n);
        tic(); 
        dot(x,y);
        res(i,n) = toc();
    end
end
res2_median = median(res);
for n = 1:vectorSize
    for i = 1:numTries
        x = rand(1, n); 
        y = rand(1, n);
        tic(); 
        x*y';
        res(i,n) = toc();
    end
end
res3_median = median(res);
x = 1:vectorSize;
plot(x, res1_median, x, res2_median, x, res3_median)
legend('my prod','dot','x*y')
title('operating time comparison')
ylabel('time');
xlabel('array size');
%% 10
clc;
A = [2 3 6 8; 2 3 5 7; 0 0 0 0];
B = [2 3 6 8; 2 3 5 7; 2 3 6 7; 0 0 0 1];
A1 = string(int2str(A));
B1 = string(int2str(B));
[Lia, Locb] = ismember(A1, B1);
[Lia, Locb1] = ismember(A,B, 'rows');
disp(Locb1 == Locb);
%% 10
clc;
A = [2 3 6 8; 2 3 5 7];
B = [2 3 6 9; 2 3 5 6; 2 3 6 8];
Locb = zeros(2,1);
R = size(A);
for i = 1:R(1)
    C = find(sum(A(i,:) == B, R(1)) == R(2));
    if isempty(C) == 0
        Locb(i) = C(1);
    end
end
[Lia, Locb1] = ismember(A,B, 'rows');
disp(Locb1 == Locb);
%% 11
a = 100;
sgm = 100;
n = 1000;
b = 200;
numTries = 100;
r = zeros(1,numTries);
for i = 1:numTries
    x = a + sgm*randn(n, 1);
    r(i) = sum(abs(x-a)>b)./length(x);
end
disp("sgm^2/b^2 = " + string(sgm^2/b^2));
disp("median result = " + string(median(r)));
%% 11
a = 100;
sgm = 100;
n = 1000;
b = 200;
x = a + sgm*randn(n, 1);
numTries = 100;
r = zeros(1,numTries);
for i = 1:numTries
    x = a + sgm*randn(n, 1);
    r(i) = sum(abs(x-a)>b)./length(x);
end
y = 1:numTries;
plot(y, r, y, repmat(sgm^2/b^2, 1, numTries));
legend('proportion', 'sgm^2/b^2');
title('inequality check');
xlabel('attempt');
ylabel('proportion');
%% 12                  
h = 0.01;                   
X = -5:h:5;       
Y = cos(X.^2);                     
n = length(X);                  
F_trapz = zeros(1,n);                 
F_rect = zeros(1,n);                
F_simps = zeros(1,n);                
for k = 2:n
    x = X(1:k);                 
    fun_val = Y(1:k);           
    F_trapz(1, k) = trapz(x, fun_val);             
    F_rect(1, k) = rect(x, fun_val);   
    F_simps(1, k) = simps(x, fun_val);     
end
plot(X, F_trapz, X, F_rect, X, F_simps); 
xlabel('x');
ylabel('g(x)');
title('comparison of methods');
legend('trapz', 'rect', 'simps');
%% 12 скорость сходимости
h = 0.001:0.0001:0.1;                                         
n = length(h);                  
F_trapz = zeros(1,n);                 
F_rect = zeros(1,n);                
F_simps = zeros(1,n);                
for k = 1:n
    X = 1:h(k):6;                                   
    fun_val = sin(X)./X;           
    F_trapz(1, k) = trapz(X, fun_val);             
    F_rect(1, k) = rect(X, fun_val);   
    F_simps(1, k) = simps(X, fun_val);
    X = 1:(h(k)/2):6;                                   
    fun_val = sin(X)./X;           
    F_trapz(1, k) = F_trapz(1, k) - trapz(X, fun_val);             
    F_rect(1, k) = F_rect(1, k) - rect(X, fun_val);   
    F_simps(1, k) = F_simps(1, k) - simps(X, fun_val);
end
plot(h, F_trapz, h, F_rect, h, F_simps); 
xlabel('h');
ylabel('g(x)');
title('comparison of convergence rate');
legend('trapz', 'rect', 'simps');
%% 12 сравнение времени вычислений
h = 0.01;                   
X = 0:h:2;                 
Y = g(X);                       
n = length(X); 
numTries = 100;
t_trapz = zeros(1,numTries);                 
t_rect = zeros(1,numTries);                
t_simps = zeros(1,numTries);
t_trapz_median = zeros(1,n);
t_rect_median = zeros(1,n);
t_simps_median = zeros(1,n);
for k = 2:n
    tic();
    for i = 1:numTries
        x = X(1:k);                 
        fun_val = Y(1:k);
        %tic();
        rect(x, fun_val); 
        %t_rect(i) = toc();
    end
    t_rect_median(k) = toc()/numTries%median(t_rect);
end
for k = 2:n
    tic()
    for i = 1:numTries
        x = X(1:k);                 
        fun_val = Y(1:k);
        %tic();
        trapz(x, fun_val);   
        %t_trapz(i) = toc();
    end
    t_trapz_median(k) = toc()/numTries%median(t_trapz);
end
for k = 2:n
    tic()
    for i = 1:numTries
        x = X(1:k);                 
        fun_val = Y(1:k);
        %tic();
        simps(x, fun_val);  
        %t_simps(i) = toc();
    end
    t_simps_median(k) = toc()/numTries%median(t_simps);
end
plot(X, t_trapz_median, X, t_rect_median, X, t_simps_median);  
legend('trapz', 'rect', 'simps');
title('comparison of calculation time');
xlabel('x');
ylabel('time');
%% 13 f(x) = sin2x f'(x) = 2cos2x
n = 1000;
h = logspace(-11, -8, n);
x = pi;
res = 2*cos(3*x) - (sin(2*(x+h)) - sin(2*(x-h)))./(2*h);
figure(1)
title('central');
plot(h, res, 'b');
xlabel('h');
legend('f`(x) - (f(x+h) - f(x-h)/2h')
figure(2)
title('right');
xlabel('h');
legend('f`(x) - (f(x+h) - f(x)/2h')
plot(h, (sin(2*(x+h)) - sin(2*x))./(2*h) - 2*cos(2*x), 'r');
%%
n = 1000;
h = logspace(-11, -8, n);
x = pi;
centr = - 3 * sin(3 * x) - (cos(3 * (x + h)) - cos(3 * (x - h))) ./ (2 * h); % f'(x) - f' centr
figure(1)
title('central');
plot(h, centr, 'r');
xlabel('h');
legend('f`(x) - (f(x+h) - f(x-h)/2h')

right = - 3 * sin(3 * x) - (cos(3 * (x + h)) - cos(3 * (x))) ./ (2 * h); % f'(x) - f' right
figure(2)
title('right');
plot(h, right, 'r');
xlabel('h');
legend('f`(x) - (f(x+h) - f(x)/2h')
%%
function res = f(x)
    res = abs(x).*cos(2.*x.^2 - 1);
end

function res = g(x)
    res = cos(x.^2);
end

function res = my_prod(x, y)
    res = 0;
    for i = 1:min(length(x), length(y))
        res = res + x(i)*y(i);
    end
end

function res = rect(x, y)
    h = diff(x);
    res = sum(y(1:end - 1).*h, 2);
end

function res = simps(x, y)
    h = x(2) - x(1);
    res = (y(1) + y(end) + 4*sum(y(2:2:(end - 1))) + 2*sum(y(3:2:(end - 2))))*h/3;
end