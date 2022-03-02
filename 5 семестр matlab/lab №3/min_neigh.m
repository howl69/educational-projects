function val = min_neigh(fld, P, i, j, var)
    [n, m] = size(P);
    val = -Inf;
    if var
        if i + 1 <= n && j + 1 <= m && P(i, j) * fld(i+1, j+1) > val
            val = P(i, j) * fld(i+1, j+1);
        end
        if i + 1 <= n && j - 1 > 0 && P(i, j) * fld(i+1, j-1) > val
            val = P(i, j) * fld(i+1, j-1);
        end
        if i - 1 > 0 && j + 1 <= m && P(i, j) * fld(i-1, j+1) > val
            val = P(i, j) * fld(i-1, j+1);
        end
        if i - 1 > 0 && j - 1 > 0 && P(i, j) * fld(i-1, j-1) > val
            val = P(i, j) * fld(i-1, j-1);
        end
    end
    if i + 1 <= n && P(i, j) * fld(i+1, j) > val
        val = P(i, j) * fld(i+1, j);
    end
    if i - 1 > 0 && P(i, j) * fld(i-1, j) > val
        val = P(i, j) * fld(i-1, j);
    end
    if j + 1 <= m && P(i, j) * fld(i, j+1) > val
        val = P(i, j) * fld(i, j+1);
    end
    if j - 1 > 0 && P(i, j) * fld(i, j-1) > val
        val = P(i, j) * fld(i, j-1);
    end
end