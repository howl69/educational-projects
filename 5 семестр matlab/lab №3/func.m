function res = func(t,x)
    %восьмёрка
    m1 = 1;
    m2 = 10;
    G = 300;
    
    %m1 = ;
    %m2 = 1;
    %G = 50;
    %m1 = 1;
    %m2 = 1;
    %G = 1000;
    res = zeros(8, 1);
    res(1:2) = x(5:6);
    res(3:4) = x(7:8);
    res(5) = G*m2/(norm(x(1:2) - x(3:4))^3) * (x(3) - x(1));
    res(6) = G*m2/(norm(x(1:2) - x(3:4))^3) * (x(4) - x(2));
    res(7) = G*m1/(norm(x(1:2) - x(3:4))^3) * (x(1) - x(3));
    res(8) = G*m1/(norm(x(1:2) - x(3:4))^3) * (x(2) - x(4));
end