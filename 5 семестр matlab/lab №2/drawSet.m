function  drawSet(rho, N)
    arrphi = linspace(0, 2*pi, N);          
    out_points = zeros(2, N);
    in_points = zeros(2, N);
    L2 = 0;         
    c2 = 0;     
    for k = 1:N
        L1 = L2;
        phi = arrphi(k);
        L2 = [cos(phi) sin(phi)]; 
        [point, val] = rho(L2);
        in_points(:,k) = point;
        if(k == 1)
            c2 = -L2(1)*in_points(1,k) - L2(2)*in_points(2,k);
        else
            c1 = c2;                
            c2 = -L2(1)*in_points(1,k) - L2(2)*in_points(2,k);
            out_points(2, k-1) =(L1(1)*c2 - L2(1)*c1)/(-L1(1)*L2(2) + L1(2)*L2(1));
            out_points(1, k-1) = -(L1(2)*c2 - L2(2)*c1)/(-L1(1)*L2(2) + L1(2)*L2(1));            
        end
    end
    out_points(2, N) = out_points(2, 1);
    out_points(1, N) = out_points(1, 1);
    out_points(1, :) = out_points(1, :);
    plot(in_points(1,:), in_points(2,:), 'r', out_points(1,:), out_points(2,:), 'b');
    legend('internal approximation', 'external approximation')
    axis([-10 10 -10 10])
end

