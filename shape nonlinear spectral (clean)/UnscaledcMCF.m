function V = UnscaledcMCF(F, V, mu, iter_num)
    L = cotmatrix(V,F);    
    
    for i = 1:iter_num
        M = massmatrix(V,F,'barycentric');

        V = (M - mu * L) \ (M * V);
        drawnow

    end    
end