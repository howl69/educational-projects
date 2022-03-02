function a = check(mat, a, n, j, cnt)
    for i = 1:n
        if mat(j, i) && isempty(find(a == i, 1))
            cnt = cnt + 1;
            a(cnt) = i;
            a = check(mat, a, n, i);
        end
    end
end