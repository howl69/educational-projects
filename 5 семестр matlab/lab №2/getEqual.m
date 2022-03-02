function [x, y] = getEqual(f, g, t0, t1, N)
    eps = 0.01;
    n = 10000;
    t = linspace(t0, t1, n+1);
    x = f(t);
    y = g(t);
    dist_matrix = pdist([x;y]', [x;y]'); 
    ind = 1:N;
    ind(end + 1) = n + 1;
    dist_points = zeros(1,N);
    for i = 1:N
        dist_points(i) = dist_matrix(ind(i), ind(i+1));
    end
    while (max(dist_points) - min(dist_points)) > eps
        max_ind = find(dist_points == max(dist_points), 1, 'last' );
        if max_ind == 1
            break
        end
        ind(max_ind) = ind(max_ind) + 1;
        dist_points(max_ind) = dist_matrix(ind(max_ind), ind(max_ind + 1));
        dist_points(max_ind - 1) = dist_matrix(ind(max_ind - 1), ind(max_ind));
    end
    x = f(t(ind));
    y = g(t(ind));
    plot(f(t), g(t), x, y, 'ro')
    title('curve splitting')
    xlabel('f(t)')
    ylabel('g(t)')
    axis equal;
    disp(dist_points(1,1));
end
