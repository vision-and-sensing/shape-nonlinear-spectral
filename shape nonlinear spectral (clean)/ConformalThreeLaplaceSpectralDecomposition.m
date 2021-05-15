function [V_smoothed, res, Phi] = ConformalThreeLaplaceSpectralDecomposition(F,V, mu, iter_num, draw_every, view_vec, big_model)
% Performs spectral filtering via the nonlinear operator "conformalized 3-Laplace" as in Brokman, Gilboa 2021

    % Requires: 
    % V, F: NX3 and MX3 matrices describing a 3D triangular mesh, where N is No. of vertices, and M is No. of triangles (faces)
    % mu: scalar step size
    % iter_num: Number of flow iterations (semi-implicit, as described in [2])
    % draw_every: A scalar that specifies when to show shape during flow
    % view_vec: a 2X1 vector of degrees - specifies viewing angle for shape during flow
    % big_model: If tru - Phi is a cell array. else a matrix.
        
    % Returns:
    % V_smoothed: NX3 matrix of vertex position after conformalized 3-Laplace flow. 
    % (triangulation F remains unchanged during flow)
    % Phi: nonlinear spectral decomposition of the mesh by our conformalized 3-Laplace (see [3])
    % res: residual of spectral decomposition, as described in [3]
    
    % Description:
    % Our conformal 3-Laplace flow, implemented for [3].
    % Ontop of the flow, this code returns nonlinear spectral decomposition 
    % of the input signal (mesh) by the Conformalized 3-Laplace operator.
    % This code embodies - 
    % 1. Adapting ideas from [1] to the Beltrami framework
    % 2. Combining above with ideas from [2]
    % 3. Spectral decomposition as in [3]
    
    % Remarks: 
    % This code uses functions from "gptoolbox" by Alec Jacobson.
    % Inputs are not modified
    
    % references:
    %  [1] Elmoatez Lezoray Bougleux 2008 
    %  [2] Kazhdan Solomon Ben-Chen 2012
    %  [3] Brokman Gilboa 2021
    
    

    %% gradient and divergence (Beltrami) matrices
    G = grad(V,F); 
    D = div(V,F);

    %% Conformal 3-Laplace Flow
    if big_model
        Phi = cell(iter_num, 1);  % Phi is a cell array (and not a matrix) to not exceed the maximal size of a variable (Phi's size would otherwise be [vertex_num X iter_num])
    else
        Phi = zeros(size(V, 1), size(V, 2), iter_num);  % Phi is a cell array (and not a matrix) to not exceed the maximal size of a variable (Phi's size would otherwise be [vertex_num X iter_num])
    end
        
    if draw_every <= iter_num
        figure;
    end
    
    for i=1:iter_num
        M = massmatrix(V,F,'barycentric');
        
        % Operator matrices: 
        % Combined gradient magnitude, as in [1]
        grads=reshape(G*V,size(F,1),size(V,2),size(V,2)); % faces X dim(3) X size(x,2)
        mag = squeeze(sqrt(sum(grads.^2,2))); % faces  X size(x,2)
        mag = squeeze(sqrt(sum(mag.^2,2))); % faces  X 1  % this line was proposed in [Elmoataz 2008]
        magMat = diag(sparse(repmat(mag,3,1)));
        % Confomalized 3-Laplace operator matrix (the actual operator is M^-1*P, not calculated explicitly due to semi-implicit flow implementation)
        P = D * magMat * G;

        % step
        V_prev = V;
        V = (M-mu*P)\(M*V);        
        
        % time derivatives and phis
        if i > 1
            Vt_prev = Vt;
        end
        Vt = V - V_prev;
        if i > 1
            Vtt = Vt - Vt_prev;
            t = i;   % * mu;
            phi = t * Vtt;
            if big_model
                Phi{i} = transpose(phi);
            else
                Phi(:,:,i) = phi;
            end
        end
        
        % show
         if mod(i,draw_every)==0
            cla;
            trisurf(F,V(:,1),V(:,2),V(:,3),zeros(size(V(:,3))));
            xlabel('x')
            zlabel('z')
            ylabel('y')
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
    if big_model
        Phi = Phi(~cellfun('isempty', Phi));
    else
        Phi(:,:,1) = [];
    end
    
    %% residual
    res = V - iter_num*Vt;
    
    %% smoothed
    V_smoothed = V;
end


