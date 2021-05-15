function [res,Phi] = ShapeSpectralTVFlow(F,V, mu, iter_num, epsilon, draw_every, view_vec)
    %% gradient and divergence (Beltrami) matrices
    M = massmatrix(V,F,'barycentric');
    G = grad(V,F); 
    D = div(V,F);
    %% Rotation Invariant TV Flow
    I = speye(size(D,1));
    Phi = cell(iter_num, 1);  % Phi is a cell array (and not a matrix) to not exceed the maximal size of a variable (Phi's size would otherwise be [vertex_num X iter_num])
    if draw_every <= iter_num
        figure;
    end
    for i=1:iter_num
        % TV matrix
        grads=reshape(G*V,size(F,1),size(V,2),size(V,2)); % faces X dim(3) X size(x,2)
        mag = squeeze(sqrt(sum(grads.^2,2))); % faces  X size(x,2)
        mag = squeeze(sqrt(sum(mag.^2,2))); % faces  X 1  % this line was proposed in []
        invGradMagVec = 1./(mag+epsilon);
        invGradMagMat = diag(sparse(repmat(invGradMagVec,3,1)));
        TVmat = D * invGradMagMat * G;

        % semi-implicit step
        V_prev = V;
        V = (M-mu*TVmat)\(M*V);
        
        % time derivatives and phis
        if i > 1
            Vt_prev = Vt;
        end
        Vt = V - V_prev;
        if i > 1
            Vtt = Vt - Vt_prev;
            t = i;   % * mu;
            phi = t * Vtt;
            Phi{i} = transpose(phi);
        end
        
        % show
         if mod(i,draw_every)==0
            cla;
            trisurf(F,V(:,1),V(:,2),V(:,3),zeros(size(V(:,3))));
            view ([view_vec]);
            shading interp;
            camlight left;
            axis equal;
            colormap jet
            title(num2str(i))
            drawnow
        end
        
    end
    % delete empty cells
    Phi = Phi(~cellfun('isempty', Phi));
    
    %% residual
    res = V - iter_num*Vt;
end

