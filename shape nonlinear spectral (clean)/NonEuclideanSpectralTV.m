function [res, Phi] = NonEuclideanSpectralTV(F, V, u, mu, iter_num, epsilon, draw_every, view_vec, big_model)
% Performs spectral TV on a non-euclidean domain.

    % Requires: 
    % V, F: NX3 and MX3 matrices describing a non-euclidean 2D domain (3D triangular mesh), where N is No. of vertices, and M is No. of triangles (faces)
    % mu: scalar step size
    % iter_num: Number of flow iterations (semi-implicit, as described in [2])
    % epsilon: denominator regularization term (for TV operator)
    % draw_every: A scalar that specifies when to show shape during flow
    % view_vec: a 2X1 vector of degrees - specifies viewing angle for shape during flow
    % big_model: If tru - Phi is a cell array. else a matrix.
        
    % Returns:
    % Phi: non-euclidean spectral TV representations
    % res: non-euclidean spectral TV residual
   
    G = grad(V,F); 
    D = div(V,F);
    M = massmatrix(V,F,'barycentric');
    
%     G = grad(V,F); 
%     D = M\div(V,F);
    
    I = speye(size(D,1));
    if big_model
        Phi = cell(iter_num, 1);  % Phi is a cell array (and not a matrix) to not exceed the maximal size of a variable (Phi's size would otherwise be [vertex_num X iter_num])
    else
        Phi = zeros(size(u, 1), size(u, 2), iter_num);  % Phi is a cell array (and not a matrix) to not exceed the maximal size of a variable (Phi's size would otherwise be [vertex_num X iter_num]);
    end
    
    if draw_every <= iter_num
        figure;
        trisurf(F,V(:,1),V(:,2),V(:,3), u);
        view ([view_vec]);
        shading interp;
        camlight left;
        axis equal;
        colormap jet
        colorbar;
        title('initial')
        drawnow
    end
    
    for i=1:iter_num
        % TV operator matrix
        TVmat = TVOperatorMatrix( V, F, u, G, D, epsilon);

        % semi-implicit step
        u_prev = u;
        u = (M-mu*TVmat)\(M*u);
        
        % time derivatives and phis
        if i > 1
            ut_prev = ut;
        end
        ut = u - u_prev;
        if i > 1
            utt = ut - ut_prev;
            t = i;   % * mu;
            phi = t * utt;
             if big_model
                Phi{i-1} = transpose(phi);
            else
                Phi(:,:,i-1) = phi;
            end
        end
        
        % show
         if mod(i,draw_every)==0
            cla;
            trisurf(F,V(:,1),V(:,2),V(:,3), u);
            view ([view_vec]);
            shading interp;
            camlight left;
            axis equal;
            colormap jet
            colorbar;
%             caxis([-2,2])
            title(num2str(i))
            drawnow
        end
        
    end    
    
    %% residual
    res = u - iter_num*ut;
end