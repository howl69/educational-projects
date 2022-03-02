function S = floyd(dist_matrix)
    N = min(size(dist_matrix));
    S = dist_matrix;
    for k=1:N
        for i=1:N
            for j=1:N
                if S(i,k)==inf || S(k,j)==inf
                    continue;
                end
                if S(i,j)>S(i,k)+S(k,j)
                    S(i,j)=S(i,k)+S(k,j);
                end
            end
        end
    end
end
