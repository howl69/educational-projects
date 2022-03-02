function Z = drawPolar(rho, N)
    arrphi = linspace(0, 2*pi, N + 1);
    [X, Y] = meshgrid(-6:0.01:6,-6:0.01:6);
    Z = ones(size(X));
    points = zeros(N, 2);
    %mov(1:N) = struct('cdata', [],'colormap', []);
    for i = 1:N
        phi = arrphi(i);
        L = [cos(phi), sin(phi)];
        [point, val] = rho(L);
        if val > 0
            points(i,:) = L / val;
        else
            [point, val] = rho(-L);
            points(i,:) = -L / val;
        end
        %contourf(X,Y,Z)
        %mov(i) = getframe();
    end
    k = convhull(points);
    array = points(k, :);
    patch(array(:,1), array(:,2), 'red');
    title('Polar')
    xlabel('X')
    ylabel('Y')
end
