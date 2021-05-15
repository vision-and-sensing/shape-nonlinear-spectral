clear all; clc; 
close all;

addpath(genpath('gptoolbox'))
addpath(genpath('meshes'))

%% flow parameters
% time parameters
iter_num = 50;
mu = 3;

%% filter parameters
% spectral bands
bands{1} = [1];
bands{2} = [2:10];
bands{3} = [11:40];  
bands{4} = [41:iter_num-1];
% per-band amplitude
amp = [1 1 1.5 1];

%% begin
% Global filter for all poses (denoted H(t) in [3])
fltr = [ones(size(bands{1}))*amp(1), ones(size(bands{2}))*amp(2), ones(size(bands{3}))*amp(3), ones(size(bands{4}))*amp(4)];

% technical parameters
big_model = false;  % set true if spectral representations exceeds MATLAB's memory limit per-variable
draw_every = 1e10;  % time-steps between shape display during flow

figure;
michael_num_vec = [1 5 12 18];  % [1 10 12 23]  % [0 1 2 4 7 10 12 17 21 23 24 25]  % 1:19
for n = 1:numel(michael_num_vec)
    %% load & centralize initial shape
    michael_num = michael_num_vec(n);
    mesh_str = ['Michael', num2str(michael_num)];
    mesh = load(mesh_str);
    vertices = @(mesh) [mesh.surface.X, mesh.surface.Y, mesh.surface.Z];
    V_orig = vertices(mesh);
    F = mesh.surface.TRIV;
    V_orig = bsxfun(@minus,V_orig,centroid(V_orig,F));

    %% display settings
    view_vec = ([0,0]);
    zoom_factor = 1.3;
    colorvec = [144, 229, 63]/256;
    
    %% show
    subplot(2, 4, n)
    trisurf(F, V_orig(:,1),V_orig(:,2),V_orig(:,3),zeros(size(V_orig,1),1));
    shading interp
    view (view_vec)
    camlight(20, 30)
    camlight(20, 30)
    material metal
    xlabel 'X'; ylabel 'Y'; zlabel 'Z';
    axis equal
    axis off
    colormap(colorvec);
    ax = gca;
    ax.Clipping = 'off';
    zoom(zoom_factor)
    drawnow

    %% Conformalized 3-Laplace filtering as in [3]
    fprintf('\nPerforming spectral analysis (conformalized 3-Laplace).. \n')
    if draw_every <= iter_num
        fprintf('Note: plotting during execution affects execution time \n')
    end

    tic;
    % spectral decomposition (flow-based)
    [~, res, phi] = ConformalThreeLaplaceSpectralDecomposition(F, V_orig, mu, iter_num, draw_every, view_vec, big_model);
    
    % filtered reconstruction
    filtered = FilteredReconstructionMatPhi(res, phi, fltr);

    % time
    execution_time = toc;
    fprintf(['Elapsed time is ', num2str(execution_time), ' seconds.\n'])

    % show
    V = filtered;
    subplot(2, 4, n + 4)
    trisurf(F, V(:,1),V(:,2),V(:,3),zeros(size(V(:,3))));
    shading interp
    view (view_vec)
    camlight(20, 30)
    camlight(20, 30)
    material metal
    xlabel 'X'; ylabel 'Y'; zlabel 'Z';
    axis equal
    axis off
    colormap(colorvec);
    ax = gca;
    ax.Clipping = 'off';
    zoom(zoom_factor)
    drawnow

end