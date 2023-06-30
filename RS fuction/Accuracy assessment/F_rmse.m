function [rmse] = F_rmse(X,endmember,S)
    [L,N] = size(X);
    rmse = sum(sqrt(sum((X - endmember * S).^2)/L))/N;
end

