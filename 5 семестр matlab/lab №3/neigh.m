function mat = neigh(mat, i, j, var)
    [n,m] = size(mat);
    if var
        if i + 1 <= n && j + 1 <= m
            mat(i+1, j+1) = 1;
        end
        if i + 1 <= n && j - 1 > 0
            mat(i+1, j-1) = 1;
        end
        if i - 1 > 0 && j + 1 <= m
            mat(i-1, j+1) = 1;
        end
        if i - 1 > 0 && j - 1 > 0
            mat(i-1, j-1) = 1;
        end
    end
    if i + 1 <= n
        mat(i+1, j) = 1;
    end
    if i - 1 > 0
        mat(i-1, j) = 1;
    end
    if j + 1 <= m
        mat(i, j+1) = 1;
    end
    if j - 1 > 0
        mat(i, j-1) = 1;
    end
end