%% 1
a = 0;
b = pi;
n = 10;
x = linspace(a, b, n);
xx = linspace(a, b, n*4);
f = @(x) sin(x.^2);
compareInterp(x, xx, f)
%% 2
figure(1) %-pchip +nearest, spline, linear
f = @(x) floor(x);
a = 0;
b = 5;
n = 20;
x = linspace(a, b, n);
xx = linspace(a, b, n*4);
compareInterp(x, xx, f)
figure(2) %-pchip +linear
f = @(x) round(x - floor(x));
a = 0;
b = 5;
n = 20;
x = linspace(a, b, n);
xx = linspace(a, b, n*4);
compareInterp(x, xx, f)
figure(3) %+spline
f = @(x) sin(10*x);
a = 0;
b = pi;
n = 20;
x = linspace(a, b, n);
xx = linspace(a, b, n*4);
compareInterp(x, xx, f)
figure(4) %-nearest
f = @(x) sin(x).*x;
a = 0;
b = pi;
n = 10;
x = linspace(a, b, n);
xx = linspace(a, b, n*4);
compareInterp(x, xx, f)
%% 3.1 малое отклонение
xmin = 0;
xmax = 2*pi;
n = 10;
n_xx = 100;
xx = linspace(xmin, xmax, n_xx);
x = linspace(xmin, xmax, n);
f = @(x) sin(2*x);
fmax = 16;
interp_res = interp1(x, f(x), xx, 'spline');
h = x(2) - x(1);
aprior = fmax * (h^4); % <= max(f''''(x))*h^4
apost = abs(f(xx) - interp_res);
plot(xx, aprior*ones(1, n_xx), xx, apost);
legend('aprior', 'apost');
%% 3.2 большое отклонение
xmin = 0;
xmax = 2*pi;
n = 10;
n_xx = 100;
xx = linspace(xmin, xmax, n_xx);
x = linspace(xmin, xmax, n);
f = @(x) x.^3;
fmax = 0;
interp_res = interp1(x, f(x), xx, 'spline');
h = x(2) - x(1);
aprior = fmax * (h^4); % <= max(f''''(x))*h^4
apost = abs(f(xx) - interp_res);
plot(xx, aprior*ones(1, n_xx), xx, apost);
legend('aprior', 'apost');
%% 4 поточечная+среднеквадратичная
fn = @(x,n) sin(x).^n;
f = @(x) x == pi/2;
a = 0;
b = pi/2;
n = 100;
convType = 0;
mov = convergenceFunc(fn,f,a,b, n, convType); %convType == 0 - среднеквадратическая, 1 - равномерная, 2 - поточечная
%%
writetask(VideoWriter('#4_1.avi'), mov)
%% 4 равномерная+поточечная
fn = @(x,n) x.^2 + 1 / n;
f = @(x) x.^2;
a = 0;
b = pi/4;
n = 100;
convType = 1;
mov = convergenceFunc(fn,f,a,b, n, convType); %convType == 0 - среднеквадратическая, 1 - равномерная, 2 - поточечная
%%
writetask(VideoWriter('#4_2.avi'), mov)
%% 5.1
f = @(x) x.^2;
a = 1;
b = 3;
n = 100;
meth = 'standart';
mov = fourierApprox(f,a,b, n, meth);
writetask(VideoWriter('#5_1.avi'), mov)
%% 5.2
f = @(x) sin(x);
a = 0;
b = 2*pi;
n = 20;
meth = 'hermite';
mov = fourierApprox(f,a,b, n, meth);
writetask(VideoWriter('#5_2.avi'), mov)
%% 5.3
f = @(x) sin(pi*x);
a = -1;
b = 1;
n = 20;
meth = 'legendre';
mov = fourierApprox(f,a,b, n, meth);
writetask(VideoWriter('#5_3.avi'), mov)
%% 6.1 example 1
f = @(x) -x - x.^2 + x.^3 + x.^4;
a = -2;
b = 2;
n = 100;
x = linspace(a, b, n);
y = f(x);
%% 6.1.2 example 2
f = @(x) -x - x.^2 + x.^3 + x.^4;
a = -3;
b = 2;
n = 100;
x = linspace(a, b, n);
y = f(x);
%% 6.1.3 example 3
f = @(x) abs(x-3);
a = 0;
b = 4;
n = 100;
x = linspace(a, b, n);
y = f(x);
%% 6.1.3 example 4
f = @(x) x.*sin(x);
a = 0;
b = 100;
n = 100;
x = linspace(a, b, n);
y = f(x);
%% 6.2
hold on
plot(x, y)
local_min = find((y(2:end-1) <= y(1:end - 2) & y(2:end-1) <= y(3:end)) == 1) + 1;
[ymax, max_n] = max(y);
if ismember(max_n + min(abs(max_n - local_min)), local_min)
    end_pos = max_n + min(abs(max_n - local_min));
else
    end_pos = max_n - min(abs(max_n - local_min));
end
plot(x(local_min), y(local_min), 'bo')
title('comet from max to min')
xlabel('x')
ylabel('f(x)')
if max_n > end_pos
    comet(flip(x(end_pos:max_n)), flip(y(end_pos:max_n)))
else
    comet(x(max_n:end_pos), y(max_n:end_pos))
end 
%% 7 example 1
f = @(t) t;
g = @(t) t.^2;
t0 = 0;
t1 = 9;
N = 20;
[x, y] = getEqual(f,g,t0,t1,N);
disp("max :" + string(max(sqrt(diff(x).^2 + diff(y).^2)) - min(sqrt(diff(x).^2 + diff(y).^2))));
disp("getEqual res: " + string(mean(sqrt(diff(x).^2 + diff(y).^2))))
disp("uniform res: " + string(mean(sqrt(diff(f(linspace(t0,t1,N+1))).^2 + diff(g(linspace(t0,t1,N+1))).^2))))
n = 10000;
disp("length: " + string(sum(sqrt(diff(f(linspace(t0,t1,n))).^2 + diff(g(linspace(t0,t1,n))).^2))))
%plot(f(t),g(t),f(res),g(res))
%% 7 example 2
A = 1;
B = 1;
alpha =  1;
delta = 1;
b = 1;
f = @(t) A*sin(alpha*t + delta);
g = @(t) B*sin(b*t);
t0 = 0;
t1 = 2*pi;
N = 20;
t = linspace(0, pi, 100);
[x y] = getEqual(f,g,t0,t1,N);
disp("max :" + string(max(sqrt(diff(x).^2 + diff(y).^2)) - min(sqrt(diff(x).^2 + diff(y).^2))));
disp("getEqual res: " + string(mean(sqrt(diff(x).^2 + diff(y).^2))))
disp("uniform res: " + string(mean(sqrt(diff(f(linspace(t0,t1,N+1))).^2 + diff(g(linspace(t0,t1,N+1))).^2))))
n = 10000;
disp("length: " + string(sum(sqrt(diff(f(linspace(t0,t1,n))).^2 + diff(g(linspace(t0,t1,n))).^2))))
%% 7 example 3
A = 2;
B = 4;
alpha =  4;
delta = 1;
b = 1;
f = @(t) A*sin(alpha*t + delta);
g = @(t) B*sin(b*t);
t0 = 0;
t1 = pi;
N = 20;
t = linspace(0, pi, 100);
[x y] = getEqual(f,g,t0,t1,N);
disp("max :" + string(max(sqrt(diff(x).^2 + diff(y).^2)) - min(sqrt(diff(x).^2 + diff(y).^2))));
disp("getEqual res: " + string(mean(sqrt(diff(x).^2 + diff(y).^2))))
disp("uniform res: " + string(mean(sqrt(diff(f(linspace(t0,t1,N+1))).^2 + diff(g(linspace(t0,t1,N+1))).^2))))
n = 10000;
disp("length: " + string(sum(sqrt(diff(f(linspace(t0,t1,n))).^2 + diff(g(linspace(t0,t1,n))).^2))))
%% 7 example 4
A = 4;
B = 2;
alpha =  5;
delta = 1;
b = 4;
f = @(t) A*sin(alpha*t + delta);
g = @(t) B*sin(b*t);
t0 = 0;
t1 = 2*pi;
N = 20;
t = linspace(0, pi, 100);
[x y] = getEqual(f,g,t0,t1,N);
disp("max :" + string(max(sqrt(diff(x).^2 + diff(y).^2)) - min(sqrt(diff(x).^2 + diff(y).^2))));
disp("getEqual res: " + string(mean(sqrt(diff(x).^2 + diff(y).^2))))
disp("uniform res: " + string(mean(sqrt(diff(f(linspace(t0,t1,N+1))).^2 + diff(g(linspace(t0,t1,N+1))).^2))))
n = 10000;
disp("length: " + string(sum(sqrt(diff(f(linspace(t0,t1,n))).^2 + diff(g(linspace(t0,t1,n))).^2))))
%% 7 example 5
A = 4;
B = 2;
alpha =  4;
delta = 1;
b = 3;
f = @(t) A*sin(alpha*t + delta);
g = @(t) B*sin(b*t);
t0 = 0;
t1 = 2*pi;
N = 20;
t = linspace(0, pi, 100);
[x y] = getEqual(f,g,t0,t1,N);
disp("max :" + string(max(sqrt(diff(x).^2 + diff(y).^2)) - min(sqrt(diff(x).^2 + diff(y).^2))));
disp("getEqual res: " + string(mean(sqrt(diff(x).^2 + diff(y).^2))))
disp("uniform res: " + string(mean(sqrt(diff(f(linspace(t0,t1,N+1))).^2 + diff(g(linspace(t0,t1,N+1))).^2))))
%% 11.1
ax = 1;
bx = 10;
ay = 1;
by = 20;
[X,Y] = meshgrid(ax:0.5:bx,ay:by);
nFrames = 100;
mov(1:nFrames) = struct('cdata', [],'colormap', []);
az = inf;
bz = -inf;
for num = 1:nFrames
    Z = sin(X*num/10) + cos(Y);
    az = min(az, min(Z(:)));
    bz = max(bz, max(Z(:)));
end
for num = 1:nFrames
    Z = sin(X*num/10) + cos(Y);
    new_Z = [ones(1, size(Z, 2) + 2)*(-inf); ones(size(Z, 1), 1)*(-inf) Z ones(size(Z, 1), 1)*(-inf); ones(1, size(Z, 2) + 2)*(-inf)];
    z_max = ones(size(Z, 1), size(Z, 2));
    for i = 1:3
        for j = 1:3
            if i ~= 2 || j ~= 2
                z_max = z_max & (new_Z(i:end-(3-i),j:end-(3-j)) < new_Z(2:end-1,2:end-1));
            end
        end
    end
    z_max = find(z_max == 1);
    new_Z = [ones(1, size(Z, 2) + 2)*(+inf); ones(size(Z, 1), 1)*(+inf) Z ones(size(Z, 1), 1)*(+inf); ones(1, size(Z, 2) + 2)*(+inf)];
    z_min = ones(size(Z, 1), size(Z, 2));
    for i = 1:3
        for j = 1:3
            if i ~= 2 || j ~= 2
                z_min = z_min & (new_Z(i:end-(3-i),j:end-(3-j)) > new_Z(2:end-1,2:end-1));
            end
        end
    end
    z_min = find(z_min == 1);
    surf(X,Y,Z)
    axis([ax bx ay by az bz])
    hold on
    [row,col] = ind2sub(size(Z), z_max);
    plot3(X(1, col), Y(row, 1), Z(z_max),'ro');
    xlabel('X')
    ylabel('Y')
    zlabel('f(X,Y)')
    title('surface evolution')
    [row,col] = ind2sub(size(Z), z_min);
    plot3(X(1, col), Y(row, 1), Z(z_min),'b*');
    hold off
    mov(num) = getframe();
end
%% 11.2
movie(mov)
%% 11.3
[X,Y] = meshgrid(1:0.5:10,1:20);
num = 1;
Z = sin(X*num/10) + cos(Y);
contour(X,Y,Z)
title('contour f(X,Y) = sin(X*num/10) + cos(Y)')
xlabel('x')
ylabel('y')
%% 12
v = VideoWriter('task11.avi');
v.FrameRate = 10;
open(v);
writeVideo(v,mov);
close(v);
%% 13
clc
%points = [0 0; 1 1; 4 3; 2 2; 3 3; 4 4];
points = [0 0; 1 1; 2 2; -1 -1; -2 -2; 0 1; 0 2; 0 -1; 0 -2];
%points = [0 0; -1 0; 1 1; 2 0; 1 -1];
%points = [0 0; 0 1; 2 1; 2 2];
L = 0.75;
x = points(:,1);
y = points(:,2);
N = length(x);
max_x = max(x) + L;
min_x = min(x) - L;
max_y = max(y) + L;
min_y = min(y) - L;
n = 100;
[X, Y] = meshgrid(linspace(min_x, max_x, n), linspace(min_y, max_y, n));
Z = ones(size(X))*(inf);
for i = 1:N
    Z = min(Z, sqrt((X - x(i)).^2 + (Y - y(i)).^2));
end
Z = -Z;
[C, h] = contourf(X,Y,Z,[-L -L]);
X = repmat(x, [1, N]);
Y = repmat(y, [1, N]);
dist_neighb = sqrt((X-X').^2+(Y-Y').^2);
dist_neighb(dist_neighb > 2*L) = inf;
S = floyd(dist_neighb);
if size(find(S == Inf), 1) == 0
    disp('Обойдет')
    min_dist = Inf;
    p = perms(2:size(points, 1));
    V = [ones(size(p, 1), 1), p];
    iter_cnt = size(V,1);
    for i = 1:iter_cnt
        dist = 0;
        for j = 1:(N-1)
            dist = dist + S(V(i,j+1), V(i,j)); 
        end
        min_dist = min(dist, min_dist);
    end
    disp(min_dist)
else
    disp('Не обойдёт :(')
end
%% 14
s.ax = 1;
s.bx = 1;
s.ay = 1;
s.by = 1;
s.az = 1;
s.bz = 1;
s.N = 10;
s.name = {'FaceColor', 'EdgeColor', 'FaceAlpha'};
s.value = {'blue', 'none', 0.75};
figure(1)
drawBall(3, 5, s)
figure(2)
s.value = {'blue', 'none', 0.75};
drawBall(Inf, 4, s)
%% 15
alphas = [0.5 1 2 5 Inf];
colors = ["blue" "green" "red" "yellow" "cyan"];
edges = ["red" "blue" "none" "black" "yellow"];
drawManyBalls(alphas, colors, edges)
%% 8
N = 10;
figure(1)
drawSet(@rho_ellipse, N);
figure(2)
drawSet(@rho_square, N);
figure(3)
drawSet(@rho_rhombus, N);
%% 10
clc
N = 40;
figure(1)
%drawPolar(@rho_ellipse, N);
%drawPolar(@rho_rhombus, N);
drawPolar(@rho_square, N);
figure(2)
%[X Y] = meshgrid(linspace(-5,5,100),linspace(-5,5,100));
X = linspace(-5,5,1000);
Y = linspace(-5,5,1000);
Z = ones(1000);
for i = 1:1000
    for j = 1:1000
        [point, Z(i,j)] = rho_square([X(j), Y(i)]);
    end
end
[X Y] = meshgrid(X,Y);
contourf (X, Y, Z < 1)
%% 9
clc
s.c = [0 0];
N = 10;
drawSet(@supportLebesgue, N);
%%
function writetask(v, mov)
    v.FrameRate = 10;
    open(v);
    writeVideo(v,mov);
    close(v);
end