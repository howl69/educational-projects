%% 1
clc
x = linspace(0, 2, 1000); 
f = @(x) x.^3 + x - 1;
g = @(x) sin(x);
h = @(x) f(x) - g(x);
plot(x, f(x), x, g(x))
[x , ~] = ginput;
for i = 1:size(x)
    x_i = fzero(h, x(i));
    disp("x_" + string(i) + " = " + string(x_i))
    disp("|h(x_" + string(i) + ")| = " + h(x_i))
end
%% 2
a = 0.1;
n = 100;
x = linspace(-a, a, n);
y = f2(x);
plot(x, y)
for i = 1:n
    y(i) = fzero(@f2, x(i));
end
hold on
plot(x, y)
%% 3
clc
A = [1 2; 3 4];
E = [1 0; 0 1];
n = 100;
res_pow = E;
for i = 1:n
    E = E*A/i;
    res_pow = res_pow + E;
end
B = [A, zeros(2,2); zeros(2,2), A];
dydt = @(t,y) B * y;
[t, res_ode] = ode45(dydt,linspace(0, 1, 1000),[1, 0, 0, 1]);
res_ode = reshape(res_ode(end,:), 2, 2);
disp(res_pow)
disp(expm(A))
disp(res_ode)
%% 4
clc
a = 0;
b = 1;
c = 1.5;
x0 = 1;
y0 = 1;
x1 = 3;
y1 = 1;
t = linspace(0, 10, 1000);
figure(1)
dxdt = @(t, x) [x(3); x(4); -a*x(1); -a*x(2)];
[t,x] = ode45(dxdt,t,[x0; y0; x1; y1]);
plot(x(:,1), x(:,2))
xmax = max(abs(x(:,1)));
ymax = max(abs(x(:,2)));
xymax = max(xmax, ymax);
axis([-xymax xymax -xymax xymax])
hold on
title('normal movement')
xlabel('x');
ylabel('y');

figure(2)
dxdt = @(t, x) [x(3); x(4); -a*x(1); -a*x(2)];
[t,x] = ode45(dxdt,t,[x0; y0; x1; y1]);
flag1 = 1;
flag2 = 1;
for i = 1:size(x,1)
    if abs(x(i, 1)) > b && flag1
        x(i:end, 1) = 2*x(i,1)-x(i:end, 1);
        flag1 = 0;
    else
        flag1 = 1;
    end
    if abs(x(i, 2)) > c && flag2
        x(i:end, 2) = 2*x(i,2)-x(i:end, 2);
        flag2 = 0;
    else
        flag2 = 1;
    end
end
plot(x(:,1), x(:,2))
plot([-b b], [c c], [-b b], [-c -c], [b b], [-c c], [-b -b], [-c c])
xymax = max(b + 0.5,c + 0.5);
axis([-xymax xymax -xymax xymax])
hold on
title('bouncing ball')
xlabel('x');
ylabel('y');
comet(x(:,1), x(:,2))
%% 5
clc
t = linspace(0, 10, 100);
[t,x] = ode45(@func,t,[0; 0; 10; 0; 0; 10; 0; -10]); %восьмёрка

%[t,x] = ode45(@func,t,[0; 0; 10; 0; 0; 10; 0; -10]);
n = size(x,1);
minx = min([x(:,1); x(:,3)]);
maxx = max([x(:,1); x(:,3)]);
miny = min([x(:,2); x(:,4)]);
maxy = max([x(:,2); x(:,4)]);
mov(1:n) = struct('cdata', [],'colormap', []);
for i = 1:n
    plot(x(1:i,1), x(1:i,2), x(1:i,3), x(1:i,4))
    axis([minx maxx miny maxy])
    mov(i) = getframe();
end
%% 6
SavePrivateRyan(1);
%% 8
clc
analytsol = @(x) exp(x)./(exp(1)-exp(-1)) - exp(-x)/(exp(1) - exp(-1)) - 2*x;
dydx = @(x,y) [y(2), y(1) + 2*x];
res = @(ya,yb) [ya(1), yb(1)+1];
guess = @(x) [sin(x), cos(x)];
xmesh = linspace(0,1,10);
solinit = bvpinit(xmesh, guess);
sol = bvp4c(dydx, res, solinit);
plot(sol.x, sol.y(1,:), '-o')
hold on
plot(sol.x, analytsol(sol.x), '-o')
title('Решение краевой задачи')
xlabel('x')
ylabel('y')
disp(sqrt(sum((sol.y(1,:) - analytsol(sol.x).^2))))
disp(max(abs(sol.y(1,:) - analytsol(sol.x))))
%% 9
clc
f = @(x) sin(x(1));
x0 = 1;
myfmin(f, x0)
%%
f = @(x) (2*x(1) + x(2)).^2 + x(2).^2 + 1;
x0 = [3 2];
myfmin(f, x0)
%%
f = @(x) (x(1) + x(2)).^2 + (x(3)-1).^2 + 1;
x0 = [1 2 3];
myfmin(f, x0)
%% 7 
clc;
V = @(x, y) x.^2 + y.^2;
n = 10;
phi = linspace(0, 2*pi, n);
x0 = [cos(phi); sin(phi)] / 2;
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
T = 0:0.1:2;
[X, Y] = meshgrid(x, y);
contour(X, Y, V(X, Y), 40, '--k');
hold on;
f = @(t, y) [y(1)^3 - y(2); y(1) + y(2)^3];
for i = 1:size(x0, 2)
    [t, y] = ode45(f, T, x0(1:2, i));
    vmax = max(V(y(:,1), y(:,2)));
    for j = 1:size(y, 1)-2
        quiver(y(j,1), y(j,2), y(j+2,1)-y(j,1), y(j+2,2)-y(j,2), 'Color', [V(y(j,1), y(j,2))/vmax 0 0]);
    end
end
xlim([-1 1])
ylim([-1 1])
%% 7.2
clc;
V = @(x, y) x.^2 + y.^4;
n = 5;
phi = linspace(0, 2*pi, n);
x0 = [cos(phi); sin(phi)] ;
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
T = 0:0.2:20;
[X, Y] = meshgrid(x, y);
contour(X, Y, V(X, Y), 20, '--k');
hold on;
f = @(t, y) [y(2)^3 - y(1)^5; -y(1) - y(2)^3 + y(2)^5];
for i = 1:size(x0, 2)
    [t, y] = ode45(f, T, x0(1:2, i));
    vmax = max(V(y(:,1), y(:,2)));
    for j = 1:size(y, 1)-2
        quiver(y(j,1), y(j,2), y(j+2,1)-y(j,1), y(j+2,2)-y(j,2), 'Color', [V(y(j,1), y(j,2))/vmax 0 0]);
    end
end
ylim([-1 1])
xlim([-1 1])
%% 10
clc
func1 = @(t) t.*exp((-2*t.^2));
func2 = @(t) atan(3*t) - atan(2*t);
func3 = @(t) (1 - cos(t.^3))./t.^5;
func4 = @(t) exp((-2*abs(t)))./(1+sin(t).^2);
ftfunc1 = @(t) exp(-t.^2/8).*(-1i*t*sqrt(2*pi)/8);
ftfunc2 = @(t) -i*pi*(exp(-abs(t)/3) - exp(-abs(t)/2))./t;
%%
fighandle1 = figure('Name','func1');
res = plotFT(fighandle1, func1, ftfunc1, 0.01, [5 10], [-5 10]);
%%
fighandle1 = figure('Name','func2');
res = plotFT(fighandle1, func2, ftfunc2, 0.01, [-10 -5], [0 10]);
%%
fighandle1 = figure('Name','func3');
res = plotFT(fighandle1, func3, [], 0.1, [-10 5], [-100 100]);
%%
fighandle1 = figure('Name','func4');
res = plotFT(fighandle1, func4, [], 0.001, [-100 100], [-20 20]);
%%
fighandle = figure('Name','func2');
subplot(2, 1, 1);
ax = gca;
ax.XLim = [-4 3];
subplot(2, 1, 2);
ax = gca;
ax.XLim = [-4 3];
SPlotInfo = struct('realXLim', 1, 'imagXLim', 1);
set(fighandle,'UserData',SPlotInfo);
res = plotFT(fighandle, func2, ftfunc2, 0.01, [-10 20], [-5 5]);
%%
fighandle = figure('Name','func2');
SPlotInfo = struct('outLimVec', [0 5], 'inpLimVec', [-10 10]);
set(fighandle,'UserData',SPlotInfo);
res = plotFT(fighandle, func2, ftfunc2, 0.01, [-10 20], [-3 5]);
