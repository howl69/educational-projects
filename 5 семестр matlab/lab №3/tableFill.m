function [fld, dir] = tableFill(P, fld, dir, i, j, var)
    [n, m] = size(fld);
    if var
        if i + 1 <= n && j + 1 <= m && fld(i+1, j+1) < fld(i, j) * P(i+1, j+1)
            fld(i+1, j+1) = fld(i, j) * P(i+1, j+1);
            dir(i+1, j+1, :) = [i;j];
            [fld, dir] = tableFill(P, fld, dir, i+1, j+1, var);
        end
        if i + 1 <= n && j - 1 > 0 && fld(i+1, j-1) < fld(i, j) * P(i+1, j-1)
            fld(i+1, j-1) = fld(i, j) * P(i+1, j-1);
            dir(i+1, j-1, :) = [i;j];
            [fld, dir] = tableFill(P, fld, dir, i+1, j-1, var);
        end
        if i - 1 > 0 && j + 1 <= m && fld(i-1, j+1) < fld(i, j) * P(i-1, j+1)
            fld(i-1, j+1) = fld(i, j) * P(i-1, j+1);
            dir(i-1, j+1, :) = [i;j];
            [fld, dir] = tableFill(P, fld, dir, i-1, j+1, var);
        end
        if i - 1 > 0 && j - 1 > 0 && fld(i-1, j-1) < fld(i, j) * P(i-1, j-1)
            fld(i-1, j-1) = fld(i, j) * P(i-1, j-1);
            dir(i-1, j-1, :) = [i;j];
            [fld, dir] = tableFill(P, fld, dir, i-1, j-1, var);
        end
    end
    if i + 1 <= n && fld(i+1, j) < fld(i, j) * P(i+1, j)
        fld(i+1, j) = fld(i, j) * P(i+1, j);
        dir(i+1, j, :) = [i;j];
        [fld, dir] = tableFill(P, fld, dir, i+1, j, var);
    end
    if i - 1 > 0 && fld(i-1, j) < fld(i, j) * P(i-1, j)
        fld(i-1, j) = fld(i, j) * P(i-1, j);
        dir(i-1, j, :) = [i;j];
        [fld, dir] = tableFill(P, fld, dir, i-1, j, var);
    end
    if j + 1 <= m && fld(i, j+1) < fld(i, j) * P(i, j+1)
        fld(i, j+1) = fld(i, j) * P(i, j+1);
        dir(i, j+1, :) = [i;j];
        [fld, dir] = tableFill(P, fld, dir, i, j+1, var);
    end
    if j - 1 > 0 && fld(i, j-1) < fld(i, j) * P(i, j-1)
        fld(i, j-1) = fld(i, j) * P(i, j-1);
        dir(i, j-1, :) = [i;j];
        [fld, dir] = tableFill(P, fld, dir, i, j-1, var);
    end
end