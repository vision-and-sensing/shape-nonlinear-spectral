clear all; clc; close all;

addpath(genpath('gptoolbox'))
addpath(genpath('meshes'))

%% per-band amplifications
num_filters = 8;
amps= cell(num_filters, 1);
amps{1} = [0 0 1];  
amps{2} = [0 0 5];  
amps{3} = [0 1 0];  
amps{4} = [0 3 0];  
amps{5} = [1 0 0];  
amps{6} = [2 0 0];  
amps{7} = [1 1 1];  
amps{8} = [2 3 5];  

%% load & centralize initial shape
[vertex,face] = read_mesh('vase-lion_original');
V_orig = vertex';
F = face';
% centralize
V_orig = bsxfun(@minus,V_orig,centroid(V_orig,F));

%% display settings
draw_every = 1e10;  % 5  % 1e10
view_vec = ([0 90]);
colorvec = [144, 229, 63]/256;
zoom_factor = 1;
PNG = [];  % intialize png-file (matrix of images of displayed figures)

%% show
figure;
trisurf(F,V_orig(:,1),V_orig(:,2),V_orig(:,3), zeros(size(V_orig, 1), 1));
view (view_vec)
shading interp
camlight (10, -10)
title('initial')
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
axis equal
axis off
colormap(colorvec);
ax = gca;
ax.Clipping = 'off';
zoom(zoom_factor)
drawnow

PNGrow1 = [];
PNGrow2 = [];
%% TV naive spectral analysis (per-coordinate prcessing)
% flow parameters
epsilon1 = 0.5*1e-3;
iter_num1 = 20;
mu1_x = 5e-4;
mu1_y = 5e-4;
mu1_z = 5e-5;

% filter bands
num_bands = 3;
bands1 = cell(num_bands, 1);
bands1 = cell(num_bands, 1);
bands1{1} = [1:3];  
bands1{2} = [4:13];  
bands1{3} = [14:iter_num1];

fprintf('\nPerforming per-coordinate spectral TV analysis.. \n')
if draw_every <= iter_num1
    fprintf('note: plotting during execution affects execution time \n')
end

filtered1 = cell(num_filters, 1);
    
big_model = false;

%% per-coordinate (naive) spectral decomposition
tic;
[resx, phix] = NonEuclideanSpectralTV(F, V_orig, V_orig(:,1), mu1_x, iter_num1, epsilon1, draw_every, view_vec, big_model);
[resy, phiy] = NonEuclideanSpectralTV(F, V_orig, V_orig(:,2), mu1_y, iter_num1, epsilon1, draw_every, view_vec, big_model);
[resz, phiz] = NonEuclideanSpectralTV(F, V_orig, V_orig(:,3), mu1_z, iter_num1, epsilon1, draw_every, view_vec, big_model);

res1 = [resx, resy, resz];
phi1 = [phix, phiy, phiz];

%% filtered reconstructions
for i=1:num_filters
    amp = amps{i};
    fprintf(['Filter No. ', num2str(i), '.. \n'])

    fltr = [];
    for j = 1:numel(amp)
       fltr = [fltr, ones(size(bands1{j}))*amp(j)] ;
    end
    
    filtered1{i} = FilteredReconstructionMatPhi(res1, phi1, fltr);
    
    % show filtered shape
    V = filtered1{i};
    figure;
    set(gcf,'color','w');
    trisurf(F, V(:,1),V(:,2),V(:,3),zeros(size(V(:,3))));
    material metal
    view (view_vec)
    shading interp
    camlight (0, 10)
    camlight (0, 30)
    xlabel 'X'; ylabel 'Y'; zlabel 'Z';
    axis equal
    axis off
    colormap(colorvec)
    ax = gca;
    ax.Clipping = 'off';
    zoom(zoom_factor)
    % show filter
    ax2 = axes('Position',[.8 .6 .15 .25]);
    box on;
    plot(1:iter_num1, fltr, 'b-', 'LineWidth', 2)
    set(gca, 'FontSize', 20);
    ylim([0 5])
    yticks([0 1 5])
    title('H(t)')
    grid on;
    drawnow

    f=getframe(gcf);
    if mod(i,2)
        PNGrow1 = [PNGrow1, f.cdata];
    else
%         PNGrow2 = [PNGrow2, f.cdata];
    end
        
    
    

end
PNG = [PNG; PNGrow1; PNGrow2];


PNGrow1 = [];
PNGrow2 = [];

%% shape TV spectral analysis
% flow parameters
epsilon2 = 0.5*1e-3;
iter_num2 = 20;
mu2 = 7e-5;

% filter bands
num_bands = 3;
bands2 = cell(num_bands, 1);
bands2 = cell(num_bands, 1);
bands2{1} = [1:3];  
bands2{2} = [4:13];  
bands2{3} = [14:iter_num1 - 1];

fprintf('\nPerforming shape TV spectral analysis.. \n')
if draw_every <= iter_num1
    fprintf('note: plotting during execution affects execution time \n')
end

filtered2 = cell(num_filters, 1);

%% shape TV spectral decomposition
[res2, phi2] = ShapeSpectralTVFlow(F, V_orig, mu2, iter_num2, epsilon2, draw_every, view_vec);

%% filtered reconstructions
for i=1:num_filters
    amplification = amps{i};
    fprintf(['Filter No. ', num2str(i), '.. \n'])
    filtered2{i} = FilteredReconstruction(res2, phi2, bands2, amplification);
    
    % show filtered shape
    V = filtered2{i};
    figure;
    set(gcf,'color','w');
    trisurf(F, V(:,1),V(:,2),V(:,3),zeros(size(V(:,3))));
    material metal  % ([0, 1, 0.9])
    view (view_vec)
    shading interp
    camlight (0, 30)
    camlight (0, 10)
    xlabel 'X'; ylabel 'Y'; zlabel 'Z';
    axis equal
    axis off
    colormap(colorvec)
    ax = gca;
    ax.Clipping = 'off';
    zoom(zoom_factor)
    drawnow
    % show filter
    fltr = [];
    for j = 1:numel(amplification)
       fltr = [fltr, ones(size(bands2{j}))*amplification(j)] ;
    end
    ax2 = axes('Position',[.8 .6 .15 .25]);
    box on;
    plot(1:numel(fltr), fltr, 'b-', 'LineWidth', 2)
    set(gca, 'FontSize', 20);
    ylim([0 5])
    yticks([0 1 5])
    title('H(t)')
    grid on;
    drawnow
    f=getframe(gcf);
    if mod(i,2)
        PNGrow1 = [PNGrow1, f.cdata];
    else
        PNGrow2 = [PNGrow2, f.cdata];
    end
        
    
    

end
PNG = [PNG; PNGrow1; PNGrow2];


imwrite(PNG, 'PNG.png')

figure; imshow(PNG); title('image to write')

fprintf(['Done. \nVertex reconstruction MSE: ', num2str(mean(mean((filtered1{i} - V_orig).^2))), '\n'])
