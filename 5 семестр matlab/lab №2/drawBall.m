function drawBall(alpha, level, params)
    x = linspace(params.ax, params.bx, params.N);
    y = linspace(params.ay, params.by, params.N);
    z = linspace(params.az, params.bz, params.N);
    [X, Y, Z] = meshgrid(x, y, z);
    if (alpha == Inf)
        f = max(max(abs(X), abs(Y)), abs(Z));
    else 
        f = (abs(X).^alpha + abs(Y).^alpha + abs(Z).^alpha).^(1./alpha);
    end
    i = isosurface(X,Y,Z,f,level);
    if (size(i.vertices) == 0)
        error('empty');
    end
    p = patch(i);
    isonormals(X, Y, Z, f, p);
    axis([params.ax params.bx params.ay params.by params.az params.bz]);
    set(p, params.name, params.value)
    daspect([1 1 1])
    view(3); 
    camlight('left');
    lighting gouraud
end