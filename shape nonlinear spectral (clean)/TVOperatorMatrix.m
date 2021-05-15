function [ TVV ] = TVOperatorMatrix( V, F, x, G, D, epsilon)
    %//
    %approximating the (E.L.) variation of Total-Variation operator to be linear by deciding the magnitude of grad(V). returns operation in matrix form
    % j is the dimension index
    %//
    grads=reshape(G*x,size(F,1),size(V,2),size(x,2)); % faces X dim(3) X size(x,2)
    mag = squeeze(sqrt(sum(grads.^2,2))); % faces  X size(x,2)
    invGradMagVec = 1./(mag+epsilon);
    invGradMagMat = diag(sparse(repmat(invGradMagVec,3,1)));
    TVV = D * invGradMagMat * G;

end

