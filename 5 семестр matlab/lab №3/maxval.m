function dir = maxval(fld, var, cur)
    maxv = -Inf;
    i = cur(1);
    j = cur(2);
    [n,m] = size(fld);
    if var
        if i + 1 <= n && j + 1 <= m && maxv < fld(i+1, j+1)
            maxv = fld(i+1, j+1);
            dir = [i+1, j+1];
        end
        if i + 1 <= n && j - 1 > 0 && maxv < fld(i+1, j-1)
            maxv = fld(i+1, j-1);
            dir = [i+1, j-1];
        end
        if i - 1 > 0 && j + 1 <= m && maxv < fld(i-1, j+1)
            maxv = fld(i-1, j+1);
            dir = [i-1, j+1];
        end
        if i - 1 > 0 && j - 1 > 0 && maxv < fld(i-1, j-1)
            maxv = fld(i-1, j-1);
            dir = [i-1, j-1];
        end
    end
    if i + 1 <= n && maxv < fld(i+1, j)
        maxv = fld(i+1, j);
        dir = [i+1, j];
    end
    if i - 1 > 0 && maxv < fld(i-1, j)
        maxv = fld(i-1, j);
        dir = [i-1, j];
    end
    if j + 1 <= m && maxv < fld(i, j+1)
        maxv = fld(i, j+1);
        dir = [i, j+1];
    end
    if j - 1 > 0 && maxv < fld(i, j-1)
        dir = [i, j-1];
    end
end