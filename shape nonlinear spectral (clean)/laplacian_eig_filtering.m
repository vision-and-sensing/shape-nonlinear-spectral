function [V_final, eVec] = laplacian_eig_filtering(vert, faces, eigenNumber, maxEigs, eig_string, do_save)
% code from http://www.mcihanozer.com/portfolio/cc/computer-graphics/thesis/mesh-smoothing-with-manifold-harmonics/

    %% STEP 1: Assemble positive semi-definite discrete Laplacian
    L = cotmatrix(vert, faces);   % Get cotangent Laplacian
    M = massmatrix(vert, faces, 'barycentric');   % Get dual area of vertices (Hodge star 0)

    Minv = sqrt(inv(M));    % Get M-1/2 for symmetry 
    beltrami = Minv * L * Minv; % Get positive semi-definite discrete Laplacian (Eqauation 2)
    beltrami = beltrami * -1;   % For positive eigenvalues

    % Handle numerical precision issue:
    % http://stackoverflow.com/a/33259074
    beltrami = (beltrami + beltrami.') * 0.5;   % Now our Laplacian is symmetric, and
                                                % its eigenvectors are orthonormal

    %% STEP 2: Compute eigenvectors for the Laplacian for getting MHB
    % Apperantly, eig() function does not give orthogonal eigenvectors
    % even if the Laplacian is positive semi-definite:
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/29459 (First entry)
    tic;
    [~, e, eVec] = svds(beltrami, maxEigs, eig_string);
    fprintf(['SVD took ', num2str(toc), ' seconds.\n'])

    % Sort eigenvectors by increasing eigenvalues (Ascending order)
    [~, I] = sort(diag(e));
    eVec = eVec(:, I);

    %% STEP 3: Map the bases into canonical basis
    Hktemp = Minv * eVec;
    Hk = Hktemp(:,1:eigenNumber);   % Take only the bases you will use
 
    %% STEP 4: Transform the mesh into frequency space (MHT)
    % For this operation, matrix multiplication is much faster than a loop
    % This operation will give you a row vector
    Xk = vert(:,1)' * M * Hk;
    Yk = vert(:,2)' * M * Hk;
    Zk = vert(:,3)' * M * Hk;

    %% STEP 5-6: Smooth the mesh and transform it back into geometry space (MHT-1)

    % Allocate memory for vertices
    [vertexNumber, ~] = size(vert); % Get vertex number
    dumVertX = zeros(vertexNumber,1);
    dumVertY = zeros(vertexNumber,1);
    dumVertZ = zeros(vertexNumber,1);

    % MHT-1
    for k = 1:eigenNumber
         dumVertX = dumVertX + (Xk(1,k) * Hk(:,k));
         dumVertY = dumVertY + (Yk(1,k) * Hk(:,k));
         dumVertZ = dumVertZ + (Zk(1,k) * Hk(:,k));
    end

    % MAP NEW VERTEX POSITIONS
    V_final = zeros(vertexNumber,3);

    if strcmp(eig_string, 'largest')
        V_final(:,1) = vert(:,1) - dumVertX;
        V_final(:,2) = vert(:,2) - dumVertY;
        V_final(:,3) = vert(:,3) - dumVertZ;

    elseif strcmp(eig_string, 'smallest')
        V_final(:,1) = dumVertX;
        V_final(:,2) = dumVertY;
        V_final(:,3) = dumVertZ;

    else
        fprintf('not implemented error')
    end




    %% save
    if do_save
        save(['armadillo_maxEigs', num2str(maxEigs), '_e_eVec', '.mat'], 'eVec');
    end

end

