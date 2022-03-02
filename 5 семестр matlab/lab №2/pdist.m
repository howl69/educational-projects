function R = pdist(X, Y)
    sz = size(X);
    n = sz(1);
    X_2 = repmat(sum(X.^2, 2), 1, n);
    Y_2 = repmat(sum(Y.^2, 2), 1, n)';
    XY = X*Y';
    R = sqrt(X_2 + Y_2 - 2*XY);
end