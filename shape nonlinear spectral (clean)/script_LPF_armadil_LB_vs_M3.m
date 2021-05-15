clear all; clc; 
close all;

addpath(genpath('gptoolbox'))
addpath(genpath('meshes'))

%% display settings
draw_every = 1e10;  % 1  % 1e10
view_vec = ([180 -90]);
colorvec = [144, 229, 63]/256;
zoom_factor = 1.5;

%% load & centralize initial shape
[vertex,face] = read_mesh('Armadillo');
V_orig = vertex;
F = face;
% centralize
V_orig = bsxfun(@minus,V_orig,centroid(V_orig,F));

%% show
figure(100);
subplot(1, 3, 1)
trisurf(F,V_orig(:,1),V_orig(:,2),V_orig(:,3), zeros(size(V_orig, 1), 1));
view (view_vec)
shading interp
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
title('Input')
drawnow

%% Directional TV spectral analysis
%% flow parameters
do_M = true;
epsilon = 5e-2;
iter_num = 10;
mu = 1;

%% filter parameters
% spectral bands
num_bands = 2;
bands = cell(num_bands, 1);
bands{1} = [1:4];  
bands{2} = [5:iter_num];  
% per-band amplitude
amps = [0 1]; 

%% begin
fprintf('\nPerforming Directional TV spectral analysis... \n')
if draw_every <= iter_num
    fprintf('note: plotting during execution affects execution time \n')
end

tic;    
% over-smoothed
V_hat = UnscaledcMCF(F, V_orig, 10, 3);

% signed magnitude and unsigned direction of deformation field between V_orig and V_hat
[f0, d] = DeformationFieldSignedMagnitudesUnsignedDirections(F, V_orig, V_hat, draw_every<=iter_num, view_vec);

% Directional TV spectral decomposition, implemented as non-euclidean spectral TV of signed deformation-field magnitude
big_model = true;
[res, Phi] = NonEuclideanSpectralTV(F, V_orig, f0, mu, iter_num, epsilon, draw_every, view_vec, big_model);
Phi{iter_num} = transpose(res);

% filtered reconstruction
filtered_signed_magnitude = FilteredReconstruction(zeros(size(V_orig,1),1), Phi, bands, amps);
filtered = V_hat+ d .* filtered_signed_magnitude;

% time
fprintf(['Done. Elapsed time is ', num2str(toc), ' seconds.\n'])

% show
V = filtered;
figure(100);
subplot(1, 3, 2)
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
title('M3 filtered')
drawnow

f=getframe(gcf);


%% Laplace-Beltrami eigenfunction filtering
%% parameters
eigenNumber = 2500; % 100 500 1000 2000 3000
maxEigs = eigenNumber;
eig_string = 'smallest';
do_save = false;
do_load = true;

%% begin
if do_load
    %%  load
    fprintf('Loading Laplace-Beltrami filtered shape (because do_load == true)... \n')
    fprintf('For (a few hours long) time comparison set do_load to false \n')
    loaded = load('BeltramiEigWorkspace.mat');
    V_linear_smoothed = loaded.V_linear_smoothed;

else
    fprintf('performing Laplace Beltrami linear filtring.. \n')
    tic;
    [V_linear_smoothed, eVec] = laplacian_eig_filtering(V_orig, F, eigenNumber, maxEigs, eig_string, do_save);
    % time
    fprintf(['Done. Elapsed time is ', num2str(toc), ' seconds.\n'])
    
    
end

figure(100);
subplot(1, 3, 3)
trisurf(F, V_linear_smoothed(:,1),V_linear_smoothed(:,2),V_linear_smoothed(:,3),zeros(size(V_orig(:,3))));
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
title('Laplace-Beltrami LPF')
