%% Параметры 1
clear
A = -[-2 0; 0 1];
B = -[3 0; 0 4];
%B = -[10 0; 0 4]
f = -[0; 0];
x_1 = [1; 2];
%x_1 = [2;0.1];
t_0 = 0;
t_1 = 1;
p = 0;
q = 0;
s = 3;
r = 1;
a = 1;
b = 2;
%%
A = -[-2 0; 0 1];
B = -[3 0; 0 4];
%B = -[10 0; 0 4];
f = -[0; 0];
x_1 = [1; 2];
%x_1 = [2;0.1];
t_0 = 0;
t_1 = 1;
p = 0;
q = 0;
s = 3;
r = 1;
a = 1;
b = 2;
%%
clear
A = -[1 0; 0 1];
B = -[1 1; 0 0];
f = -[1; 1];
x_1 = [1; 2];
t_0 = 0;
t_1 = 2;
p = 1;
q = 1;
s = 6;
r = 3;
a = 2;
b = 5;
%% 2
clear
A = -[0 -3; 2 0];
B = -[1 2; 0 4];
f = -[0; 1];
x_1 = [-2; 0];
t_0 = 0;
t_1 = 2;
p = 1;
q = 1;
s = 6;
r = 3;
a = 2;
b = 5;
%% 3
clear
A = -[0 -10; 7 0];
B = -[2 0; 0 2];
f = -[0; 0];
x_1 = [0; 1.7];
t_0 = 0;
t_1 = 1;
p = 0;
q = 0;
s = 1;
r = 1;
a = 2;
b = 1;
%% Основной алгоритм
% множества X_0 X_1
hold on
plot(x_1(1),x_1(2),'r*')
axis([-5 5 -5 5])
t = linspace(0, 2*pi, 100);
x_ell = p + sqrt(r/s)*cos(t);
y_ell = q + sqrt(r)*sin(t);
plot(x_ell, y_ell, 'k') 
x_rh = [-r, 0, r, 0, -r];
y_rh = [0, r/s, 0, -r/s, 0];
plot(x_rh, y_rh, 'k')

n = 12;
min_t = inf;
arr_x = [];
arr_y = [];
arr_u1 = [];
arr_u2 = [];
arr_psi1 = [];
arr_psi2 = [];
arr_t = [];
for i = 0:n-1
    psi_0 = [cos(i*2*pi/n); sin(i*2*pi/n)];
    psi = @(t) cell2mat(arrayfun(@(x) expm(-A'*x)*psi_0, t, 'UniformOutput', false));
    l_1 = @(t) B(1:2)*psi(t);
    l_2 = @(t) B(3:4)*psi(t);
    u = @(t) u_find(l_1, l_2, a, b, t); %u из условия максимума
    func = @(t,x) A*x + B*u(t) + f;
    options = odeset('Events', @(t, x) stop(t, x, p, q, r, s),'MaxStep', 0.01);
    [t,x] = ode45(func,[t_0 t_1],x_1, options);
    %if t_1 - t(end) > 0.001 %раскомментировать если нужны только траектории, дошедшие до множества
        %запоминаем все вычисленные параметры
        res_u = u(t)';
        res_psi = psi(t')';
        arr_psi1 = [[ones(max(0, size(res_psi(:,1),1) - size(arr_psi1,1)), size(arr_psi1,2))*NaN;arr_psi1] [ones(max(0, size(arr_psi1,1) - size(res_psi(:,1),1)), size(res_psi(:,1),2))*NaN;res_psi(:,1)]]; 
        arr_psi2 = [[ones(max(0, size(res_psi(:,2),1) - size(arr_psi2,1)), size(arr_psi2,2))*NaN;arr_psi2] [ones(max(0, size(arr_psi2,1) - size(res_psi(:,2),1)), size(res_psi(:,2),2))*NaN;res_psi(:,2)]];
        arr_u1 = [[ones(max(0, size(res_u(:,1),1) - size(arr_u1,1)), size(arr_u1,2))*NaN;arr_u1] [ones(max(0, size(arr_u1,1) - size(res_u(:,1),1)), size(res_u(:,1),2))*NaN;res_u(:,1)]]; 
        arr_u2 = [[ones(max(0, size(res_u(:,2),1) - size(arr_u2,1)), size(arr_u2,2))*NaN;arr_u2] [ones(max(0, size(arr_u2,1) - size(res_u(:,2),1)), size(res_u(:,2),2))*NaN;res_u(:,2)]]; 
        arr_x = [[ones(max(0, size(x(:,1),1) - size(arr_x,1)), size(arr_x,2))*x_1(1);arr_x] [ones(max(0, size(arr_x,1) - size(x(:,1),1)), size(x(:,1),2))*x_1(1);x(:,1)]]; 
        arr_y = [[ones(max(0, size(x(:,2),1) - size(arr_y,1)), size(arr_y,2))*x_1(2);arr_y] [ones(max(0, size(arr_y,1) - size(x(:,2),1)), size(x(:,2),2))*x_1(2);x(:,2)]]; 
        arr_t = [[zeros(max(0, size(t, 1) - size(arr_t, 1)), size(arr_t, 2)); arr_t] [zeros(max(0, size(arr_t,1) - size(t,1)), size(t,2));t]];
        if t(end) < min_t
            %запоминаем отдельно результаты на оптимальной траектории
            min_t = t(end);
            arr_min_t = t;
            min_x = x;
            min_u = res_u;
            min_psi = res_psi;
            min_psi0 = psi_0;
            min_phi = i*2*pi/n;
        end
    %end
end
plot(arr_x, arr_y, 'b')
plot(min_x(:,1),min_x(:,2), 'r*')
xlabel('x1')
ylabel('x2')
disp(min_t)
%% Глобальное улучшение точности
hold on
plot(x_1(1),x_1(2),'r*')
axis([-5 5 -5 5])
plot(x_ell, y_ell,'k')
plot(x_rh, y_rh,'k')

k = 2;
n = k*n;
for i = 0:n-1
    if mod(i,k) == 0 %пропускаем уже рассмотренные направления psi_0 
        continue;
    end
    psi_0 = [cos(i*2*pi/n); sin(i*2*pi/n)];
    psi = @(t) cell2mat(arrayfun(@(x) expm(-A'*x)*psi_0, t, 'UniformOutput', false));
    l_1 = @(t) B(1:2)*psi(t);
    l_2 = @(t) B(3:4)*psi(t);
    u = @(t) u_find(l_1, l_2, a, b, t); 
    func = @(t,x) A*x + B*u(t) + f;
    options = odeset('Events', @(t, x) stop(t, x, p, q, r, s),'MaxStep', 0.01);
    [t,x] = ode45(func,[t_0 t_1],x_1, options);
    %if t_1 - t(end) > 0.001
        res_u = u(t)';
        res_psi = psi(t')';
        arr_psi1 = [[ones(max(0, size(res_psi(:,1),1) - size(arr_psi1,1)), size(arr_psi1,2))*NaN;arr_psi1] [ones(max(0, size(arr_psi1,1) - size(res_psi(:,1),1)), size(res_psi(:,1),2))*NaN;res_psi(:,1)]]; 
        arr_psi2 = [[ones(max(0, size(res_psi(:,2),1) - size(arr_psi2,1)), size(arr_psi2,2))*NaN;arr_psi2] [ones(max(0, size(arr_psi2,1) - size(res_psi(:,2),1)), size(res_psi(:,2),2))*NaN;res_psi(:,2)]];
        arr_u1 = [[ones(max(0, size(res_u(:,1),1) - size(arr_u1,1)), size(arr_u1,2))*NaN;arr_u1] [ones(max(0, size(arr_u1,1) - size(res_u(:,1),1)), size(res_u(:,1),2))*NaN;res_u(:,1)]]; 
        arr_u2 = [[ones(max(0, size(res_u(:,2),1) - size(arr_u2,1)), size(arr_u2,2))*NaN;arr_u2] [ones(max(0, size(arr_u2,1) - size(res_u(:,2),1)), size(res_u(:,2),2))*NaN;res_u(:,2)]]; 
        arr_x = [[ones(max(0, size(x(:,1),1) - size(arr_x,1)), size(arr_x,2))*x_1(1);arr_x] [ones(max(0, size(arr_x,1) - size(x(:,1),1)), size(x(:,1),2))*x_1(1);x(:,1)]]; 
        arr_y = [[ones(max(0, size(x(:,2),1) - size(arr_y,1)), size(arr_y,2))*x_1(2);arr_y] [ones(max(0, size(arr_y,1) - size(x(:,2),1)), size(x(:,2),2))*x_1(2);x(:,2)]]; 
        arr_t = [[zeros(max(0, size(t, 1)) - size(arr_t, 1), size(arr_t, 2)); arr_t] [zeros(max(0, size(arr_t,1) - size(t,1)), size(t,2));t]];
        if t(end) < min_t
            min_t = t(end);
            arr_min_t = t;
            min_x = x;
            min_psi = res_psi;
            min_u = res_u;
            min_psi0 = psi_0;
            min_phi = i*2*pi/n;
        end
    %end
end
plot(arr_x, arr_y,'b')
plot(min_x(:,1),min_x(:,2), 'r*')
xlabel('x1')
ylabel('x2')
disp(min_t)
%% Локальное улучшение точности
hold on
plot(x_1(1),x_1(2),'r*')
axis([-5 5 -5 5])
plot(x_ell, y_ell,'k')
plot(x_rh, y_rh,'k')

delta = 0.1; 
% delta = 0.01; %для еще большей точности 
n = 50;
arr_phi = linspace(min_phi - delta, min_phi + delta, n);
for i = 1:n
    psi_0 = [cos(arr_phi(i)); sin(arr_phi(i))];
    psi = @(t) cell2mat(arrayfun(@(x) expm(-A'*x)*psi_0, t, 'UniformOutput', false));
    l_1 = @(t) B(1:2)*psi(t);
    l_2 = @(t) B(3:4)*psi(t);
    u = @(t) u_find(l_1, l_2, a, b, t); 
    func = @(t,x) A*x + B*u(t) + f;
    options = odeset('Events', @(t, x) stop(t, x, p, q, r, s),'MaxStep', 0.01);
    [t,x] = ode45(func,[t_0 t_1],x_1, options);
    if t_1 - t(end) > 0.001 
        res_u = u(t)';
        res_psi = psi(t')';
        arr_psi1 = [[ones(max(0, size(res_psi(:,1),1) - size(arr_psi1,1)), size(arr_psi1,2))*NaN;arr_psi1] [ones(max(0, size(arr_psi1,1) - size(res_psi(:,1),1)), size(res_psi(:,1),2))*NaN;res_psi(:,1)]]; 
        arr_psi2 = [[ones(max(0, size(res_psi(:,2),1) - size(arr_psi2,1)), size(arr_psi2,2))*NaN;arr_psi2] [ones(max(0, size(arr_psi2,1) - size(res_psi(:,2),1)), size(res_psi(:,2),2))*NaN;res_psi(:,2)]];
        arr_u1 = [[ones(max(0, size(res_u(:,1),1) - size(arr_u1,1)), size(arr_u1,2))*NaN;arr_u1] [ones(max(0, size(arr_u1,1) - size(res_u(:,1),1)), size(res_u(:,1),2))*NaN;res_u(:,1)]]; 
        arr_u2 = [[ones(max(0, size(res_u(:,2),1) - size(arr_u2,1)), size(arr_u2,2))*NaN;arr_u2] [ones(max(0, size(arr_u2,1) - size(res_u(:,2),1)), size(res_u(:,2),2))*NaN;res_u(:,2)]]; 
        arr_x = [[ones(max(0, size(x(:,1),1) - size(arr_x,1)), size(arr_x,2))*x_1(1);arr_x] [ones(max(0, size(arr_x,1) - size(x(:,1),1)), size(x(:,1),2))*x_1(1);x(:,1)]]; 
        arr_y = [[ones(max(0, size(x(:,2),1) - size(arr_y,1)), size(arr_y,2))*x_1(2);arr_y] [ones(max(0, size(arr_y,1) - size(x(:,2),1)), size(x(:,2),2))*x_1(2);x(:,2)]]; 
        arr_t = [[zeros(max(0, size(t, 1)) - size(t, 1), size(t, 2)); arr_t] [zeros(max(0, size(arr_t,1) - size(t,1)), size(t,2));t]];
        if max(t) < min_t
            min_t = t(end);
            arr_min_t = t;
            min_x = x;
            min_u = res_u;
            min_psi = res_psi;
            min_psi0 = psi_0;
            min_phi = arr_phi(i);
        end
    end
end
plot(arr_x, arr_y, 'b')
plot(min_x(:,1),min_x(:,2), 'r*')
xlabel('x1')
ylabel('x2')
disp(min_t)
%% Результаты
hold on
p1 = plot(x_ell, y_ell,'k');
plot(x_rh, y_rh,'k')
p2 = plot(min_x(:,1),min_x(:,2), 'r*');
fprintf("max_t = %f\n",min_t)
fprintf("max_psi0 = %f %f\n", min_psi0)
fprintf("max_phi = %f\n", min_phi)
% проверка выполнения второго условия трансверсальности
rate = @(v1,v2) 1 - (-v1'*v2)/(norm(v1)*norm(v2)); %функция 1-cos
psi = @(t) expm(-A'*t)*min_psi0;
v1 = psi(min_t);
last = min_x(end, :);
p3 = quiver(last(1), last(2), -v1(1)/norm(v1), -v1(2)/norm(v1), 'g'); %вектор -psi(t1)
%весь блок if - поиск нормали в точке (x(t1),y(t1)) для множества X_1(или X_0)
%в моём случае 3 варианта: попадение на границу ромба, границу эллипса, и
%в точку пересечения ромба и эллипса - соответственно для каждого варианта
%своя нормаль
if max(abs(s*(last(1)-p)^2+(last(2)-q)^2 - r), abs(last(1)) + s*abs(last(2)) - r) < 0.02
    alpha_1 = toangle(s*(last(1) - p)/r,(last(2) - q)/r);
    alpha_2 = toangle(sign(last(1))*r/s, sign(last(2))*r);
    p4 = quiver(last(1), last(2), cos(alpha_1), sin(alpha_1), 'c');
    quiver(last(1), last(2), cos(alpha_2), sin(alpha_2), 'c')
    alpha = toangle(-v1(1), -v1(2));
    if abs(alpha_2 - alpha_1) > pi
        if alpha_1 > alpha_2
            alpha_2 = alpha_2 + 2*pi;
        else
            alpha_1 = alpha_1 + 2*pi;
        end
    end
    if (alpha < max(alpha_1, alpha_2) && alpha > min(alpha_1, alpha_2)) || (2*pi + alpha < max(alpha_1, alpha_2) && 2*pi + alpha > min(alpha_1, alpha_2))
        v2 = -v1;
    elseif rate(v1, [s*(last(1) - p)/r;(last(2) - q)/r]) < rate(v1, [sign(last(1))*r/s; sign(last(2))*r])
        v2 = [s*(last(1) - p)/r;(last(2) - q)/r];
    else
        v2 = [sign(last(1))*r/s; sign(last(2))*r];
    end
elseif abs(s*(last(1)-p)^2+(last(2)-q)^2 - r) < 0.02
    v2 = [s*(last(1) - p)/r; (last(2) - q)/r];
    p4 = quiver(last(1), last(2), v2(1)/norm(v2), v2(2)/norm(v2), 'c');
else
    v2 = [sign(last(1))*r/s; sign(last(2))*r];
    p4 = quiver(last(1), last(2), v2(1)/norm(v2), v2(2)/norm(v2), 'c');
end
legend([p1, p2, p3, p4], {'X1', 'Optimal', '-psi(t1}', 'n'})
xlabel('x1')
ylabel('x2')
fprintf("1 - cos(psi(t1), n) = %f\n", rate(v1,v2))
%% u1 u2
t = linspace(-sqrt(b/a), sqrt(b/a), n);
hold on
axis([-sqrt(b/a) sqrt(b/a) -5 5])
plot(t, a*t.^2-b, t, -a*t.^2+b)
plot(arr_u1, arr_u2, '*')
xlabel('u1')
ylabel('u2')
%% x1 x2
hold on
plot(x_1(1),x_1(2),'*')
axis([-5 5 -5 5])
t = linspace(0, 2*pi, 100);
x_ell = p + sqrt(r/s)*cos(t);
y_ell = q + sqrt(r)*sin(t);
plot(x_ell, y_ell, 'k')
x_rh = [-r, 0, r, 0, -r];
y_rh = [0, r/s, 0, -r/s, 0];
plot(x_rh, y_rh, 'k')
plot(arr_x, arr_y, 'b')
plot(min_x(:,1),min_x(:,2), 'r*')
xlabel('x1')
ylabel('x2')
%% t x1, t x2
figure(1)
hold on
plot(max(arr_t) - arr_t, arr_x)
plot(max(arr_min_t) - arr_min_t, min_x(:,1), '*')
xlabel('t')
ylabel('x1')
figure(2)
hold on
plot(max(arr_t) - arr_t, arr_y)
plot(max(arr_min_t) - arr_min_t, min_x(:,2), '*')
xlabel('t')
ylabel('x2')
%% t u1, t u2
figure(1)
hold on
plot(max(arr_t) - arr_t, arr_u1)
plot(max(arr_min_t) - arr_min_t, min_u(:,1), '*')
xlabel('t')
ylabel('u1')
figure(2)
hold on
plot(max(arr_t) - arr_t, arr_u2)
plot(max(arr_min_t) - arr_min_t, min_u(:,2), '*')
xlabel('t')
ylabel('u2')
%% t psi1, t psi2
figure(1)
hold on
plot(max(arr_t) - arr_t, arr_psi1)
plot(max(arr_min_t) - arr_min_t, min_psi(:,1), '*')
xlabel('t')
ylabel('psi1')
figure(2)
hold on
plot(max(arr_t) - arr_t, arr_psi2)
plot(max(arr_min_t) - arr_min_t, min_psi(:,2), '*')
xlabel('t')
ylabel('psi2')
%%
function u = u_find(l_1, l_2, a, b, t)
    u = [];
    for i = 1:size(t, 1) %цикл для того чтобы t мог быть вектором, мб и не нужно:)
        if abs(l_1(t(i))) < 0.01 % случай (0, l_2) и (0,0) одновременно
            u1 = [0; sign(l_2(t(i)))*b]; %формулка из опорной функции P
        else 
            %аналогично формулка из опорной функции P
            u1 = [1;1];
            u1(1) = sign(l_1(t(i)))*(-abs(l_2(t(i))/l_1(t(i)))+sqrt(l_2(t(i))^2/l_1(t(i))^2+4*a*b))/(2*a);
            u1(2) = sign(l_2(t(i)))*abs(l_2(t(i))/l_1(t(i)))*(-abs(l_2(t(i))/l_1(t(i)))+sqrt(l_2(t(i))^2/l_1(t(i))^2+4*a*b))/(2*a);
        end
        u = [u u1]; 
    end
end
function alpha = toangle(x,y) %вектор -> угол
    if y > 0
        alpha = acos(x/sqrt(x^2+y^2));
    else
        alpha = 2*pi - acos(x/sqrt(x^2+y^2));
    end
end

function [value,isterminal,direction] = stop(t, x, p, q, r, s)
    if (max(s*(x(1) - p)^2 + (x(2) - q)^2, abs(x(1)) + s*abs(x(2))) <= r - 0.01)
        value = 1;
    else
        value = 0;
    end
    isterminal = 1;   
    direction = 0;  
end