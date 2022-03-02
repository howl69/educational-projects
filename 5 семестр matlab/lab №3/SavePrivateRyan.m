function SavePrivateRyan(var)
    n = 5;
    m = 7;
    GenerateTable(n,m);
    load('table.mat')
    fld = -inf*ones(n, m);
    fld(n,m) = 1;
    change = zeros(n,m);
    change(n-1, m) = 1;
    change(n, m-1) = 1;
    while 1
        new_change = zeros(n,m);
        for i = n:-1:1
            for j = m:-1:1
                if change(i,j)
                    val = min_neigh(fld, P, i, j, var);
                    if val > fld(i,j)
                        new_change = neigh(new_change, i, j, var);
                        fld(i,j) = val;
                    end
                end
            end
        end
        if min(new_change == change)
            break
        end
        change = new_change; 
    end
    disp(fld)
    figure(2)
    image(fld,'CDataMapping','scaled')
    colorbar
    disp(fld)
    hold on
     for i = 1:n
         for j = 0:(m-1)
             text(0.6 + j,i - 0.2, string(fld(i,j + 1)))
         end
     end
    direct = [1, 1];
    while min(direct == [n, m]) ~= 1
        new_direct = maxval(fld, var, direct);
        plot([direct(2), new_direct(2)], [direct(1), new_direct(1)], 'r');
        direct = new_direct;
    end
end