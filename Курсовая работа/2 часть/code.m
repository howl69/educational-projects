%% параметрический портрет
t = linspace(0,sqrt(2)-1.05,100);
y = @(t) 4*(1-t).*(t+1).^3./((t.^2+2*t-1).^2);
plot(t,y(t), 'r')
hold on
t = linspace(sqrt(2) - 0.95, 1.1, 100);
plot(t, y(t), 'r')
plot([sqrt(2) - 1, sqrt(2) - 1], [0, 1000], 'r', [1, 1], [0, 1000], 'r')
axis([0 1.1 0 100])
xlabel('$\gamma$', 'interpreter', 'latex')
ylabel('$\alpha$', 'interpreter', 'latex')
text(0.1, 50, '1', 'FontSize',16)
text(0.3, 10, '2', 'FontSize',16)
text(0.5, 10, '3', 'FontSize',16)
text(0.7, 50, '4', 'FontSize',16)
text(1.05, 50, '5', 'FontSize',16)
%% фазовый 1
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 20;
gam = 0.2;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, gam, 1], [0, alp*gam*(1-gam)/(1+gam), 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')
%% фазовый 2
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 10;
gam = 0.2;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, gam, 1], [0, alp*gam*(1-gam)/(1+gam), 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')
%% фазовый 3
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 20;
gam = 0.45;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, gam, 1], [0, alp*gam*(1-gam)/(1+gam), 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')
%% фазовый 4
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 30;
gam = 0.7;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, gam, 1], [0, alp*gam*(1-gam)/(1+gam), 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')
%% фазовый 5
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 20;
gam = 1.1;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, 1], [0, 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')

%%
x = sym('x');
y = sym('y');
a = sym('a', 'real');
assume(a > 0);
g = sym('g');
dotx = a*x^2*(1-x)/(1+x)-x*y;
doty = -g*y+x*y;
dotx = subs(dotx, g, sqrt(2) - 1);
doty = subs(doty, g, sqrt(2) - 1);
%% 1 var
m_1 = sym('m_1');
n_1 = sym('n_1');
k_1 = sym('k_1');
m_2 = sym('m_2');
n_2 = sym('n_2');
k_2 = sym('k_2');
B1 = diff(diff(dotx, 'x'),'x')*m_1*n_1 + diff(diff(dotx, 'x'),'y')*m_1*n_2 + diff(diff(dotx, 'y'),'x')*m_2*n_1 + diff(diff(dotx, 'y'),'y')*m_2*n_2;
B2 = diff(diff(doty, 'x'),'x')*m_1*n_1 + diff(diff(doty, 'x'),'y')*m_1*n_2 + diff(diff(doty, 'y'),'x')*m_2*n_1 + diff(diff(doty, 'y'),'y')*m_2*n_2;
C1 = diff(diff(diff(dotx, 'x'),'x'),'x')*m_1*n_1*k_1 + diff(diff(diff(dotx, 'x'),'x'),'y')*m_1*n_1*k_2 + diff(diff(diff(dotx, 'x'),'y'),'x')*m_1*n_2*k_1 + diff(diff(diff(dotx, 'x'),'y'),'y')*m_1*n_2*k_2 + diff(diff(diff(dotx, 'y'),'x'),'x')*m_2*n_1*k_1 + diff(diff(diff(dotx, 'y'),'x'),'y')*m_2*n_1*k_2 + diff(diff(diff(dotx, 'y'),'y'),'x')*m_2*n_2*k_1 + diff(diff(diff(dotx, 'y'),'y'),'y')*m_2*n_2*k_2;
C2 = diff(diff(diff(doty, 'x'),'x'),'x')*m_1*n_1*k_1 + diff(diff(diff(doty, 'x'),'x'),'y')*m_1*n_1*k_2 + diff(diff(diff(doty, 'x'),'y'),'x')*m_1*n_2*k_1 + diff(diff(diff(doty, 'x'),'y'),'y')*m_1*n_2*k_2 + diff(diff(diff(doty, 'y'),'x'),'x')*m_2*n_1*k_1 + diff(diff(diff(doty, 'y'),'x'),'y')*m_2*n_1*k_2 + diff(diff(diff(doty, 'y'),'y'),'x')*m_2*n_2*k_1 + diff(diff(diff(doty, 'y'),'y'),'y')*m_2*n_2*k_2;
p1 = 1i/(2*sqrt(a*(sqrt(2)-1)));
p2 = 0.5;
q1 = 1i*sqrt(a*(sqrt(2)-1));
q2 = 1;
B1 = subs(B1, {x,y}, [0,0]);
B2 = subs(B2, {x,y}, [0,0]);
C1 = subs(C1, {x,y}, [0,0]);
C2 = subs(C2, {x,y}, [0,0]);
g20 = conj(p1)*subs(B1, {m_1,n_1,m_2,n_2}, [q1, q1, q2, q2]) + conj(p2)*subs(B2, {m_1,n_1,m_2,n_2}, [q1, q1, q2, q2]);
g11 = conj(p1)*subs(B1, {m_1,n_1,m_2,n_2}, [q1, conj(q1), q2, conj(q2)]) + conj(p2)*subs(B2, {m_1,n_1,m_2,n_2}, [q1, conj(q1), q2, conj(q2)]);
g21 = conj(p1)*subs(C1, {m_1,n_1,k_1,m_2,n_2,k_2}, [q1, q1, conj(q1), q2, q2, conj(q2)]) + conj(p2)*subs(C2, {m_1,n_1,k_1,m_2,n_2,k_2}, [q1, q1, conj(q1), q2, q2, conj(q2)]);
w0 = sqrt(a*(sqrt(2)-1)^3);
l1 = real(i*g20*g11+w0*g21)/(2*w0);
simplify(vpa(l1))
%% 2 var
p1 = 1i/(2*sqrt(a*(sqrt(2)-1)));
p2 = 0.5;
q1 = 1i*sqrt(a*(sqrt(2)-1));
q2 = 1;
syms z w
g = sqrt(2) - 1;
newx = z*q1+w*conj(q1);
newy = z*q2+w*conj(q2);
G = conj(p1)*subs(dotx, {x, y}, [newx, newy]) + conj(p2)*subs(doty, {x,y}, [newx, newy]);
g20 = diff(diff(G, 'z'),'z');
g11 = diff(diff(G, 'z'),'w');
g21 = diff(diff(diff(G, 'z'), 'z'), 'w');
w0 = sqrt(a*(sqrt(2)-1)^3);
l1 = subs(real(i*g20*g11+w0*g21), {z,w}, [0,0])/(2*w0);
simplify(vpa(l1))
%% фазовый, бифуркация
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 60;
limx = 1;
limy = 1;
x = linspace(0,limx,n);
y = linspace(0,limy,n);
alp = 1;
gam = sqrt(2) - 1;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot(gam, alp*gam*(1-gam)/(1+gam), 'r*')
axis([0 limx 0 limy])
xlabel('x')
ylabel('y')
%% возникновение предельного цикла
alp = 3;
gam = sqrt(2) - 1;
x0 = gam + 0.1;
y0 = alp*gam*(1-gam)/(1+gam);
f = @(t, x) [alp*x(1)^2*(1-x(1))/(1+x(1))-x(1)*x(2); -gam*x(2)+x(1)*x(2)];
[t,x] = ode45(f,[0 1000],[x0; y0]);
plot(x(:,1),x(:,2), 'b', x0,y0, 'r*')
xlabel('x')
ylabel('y')