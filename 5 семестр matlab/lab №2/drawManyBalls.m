function drawManyBalls(alphas, colors, edges)
    a = -1.1; 
    b = 1.1;
    n = 40;
    x = linspace(a,b,n);
    y = linspace(a,b,n);
    z = linspace(a,b,n);
    [X,Y,Z] = meshgrid(x,y,z);
    for i = 1:length(alphas)
        subplot(2,floor((length(alphas)+1)/2),i);
        if (alphas(i) == Inf)
            f = max(max(abs(X), abs(Y)), abs(Z));
        else 
            f = (abs(X).^alphas(i) + abs(Y).^alphas(i) + abs(Z).^alphas(i)).^(1./alphas(i));
        end
        s = isosurface(X,Y,Z,f,1);
        p = patch(s);
        isonormals(x,y,z,f,p)
        view(3);
        title (num2str(alphas(i)));
        set(p, 'FaceColor', colors(i), 'EdgeColor', edges(i));
        lighting gouraud
    end
end