clc;clear all;
%% set up
addpath("matlab\code");

%% 1.load test data and plot
activation = load("matlab\data\test\MEEG_source_beta_14Hz_31Hz.mat");
J = activation.J_FSAve;
Cortex8K = load("matlab\data\test\tess_cortex_concat_8000V_fix.mat");
sources_iv              = sqrt(abs(J));
sources_iv              = sources_iv/max(sources_iv(:));
smoothValue             = 0.66;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Cortex8K.Vertices, smoothValue, SurfSmoothIterations, Cortex8K.VertConn, 1);
fig1 = plotBrainSurface(sources_iv, Vertices, Cortex8K.Faces);

%% 2.unsample and plot
Cortex32K = load("matlab\data\test\tess_cortex_concat.mat");

upsample_data = upsampleData(J, Cortex8K.Vertices, Cortex32K.Vertices);
sources_iv              = sqrt(abs(upsample_data));
sources_iv              = sources_iv/max(sources_iv(:));
smoothValue             = 0.66;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Cortex32K.Vertices, smoothValue, SurfSmoothIterations, Cortex32K.VertConn, 1);
fig2 = plotBrainSurface(sources_iv, Vertices, Cortex32K.Faces);

%% transformtion  and plot
target_sphere = struct('lh', 'matlab/data/fsaverage10k/tpl-fsaverage_den-10k_hemi-L_sphere.surf.gii', ...
                       'rh', 'matlab/data/fsaverage10k/tpl-fsaverage_den-10k_hemi-R_sphere.surf.gii');
source_sphere = struct('lh', 'matlab/data/test/L.sphere.32k_fs_LR.surf.gii', ...
                       'rh', 'matlab/data/test/R.sphere.32k_fs_LR.surf.gii');
resampled_data = transformationBrainData(upsample_data, source_sphere, target_sphere);

left_surf = gifti('C:\Users\Administrator\neuromaps-data\atlases\fsaverage\tpl-fsaverage_den-10k_hemi-L_pial.surf.gii');
right_surf = gifti('C:\Users\Administrator\neuromaps-data\atlases\fsaverage\tpl-fsaverage_den-10k_hemi-R_pial.surf.gii');
% 提取左右脑的顶点和面片
left_vertices = left_surf.vertices;  % 左脑顶点坐标
left_faces = left_surf.faces;        % 左脑面片索引

right_vertices = right_surf.vertices; % 右脑顶点坐标
right_faces = right_surf.faces;       % 右脑面片索引

% 合并顶点：直接拼接
Vertices10 = double([left_vertices; right_vertices]);

% 合并面片：右脑面片索引需要偏移左脑顶点数量
n_left_vertices = size(left_vertices, 1);
Faces10 = double([left_faces; right_faces + n_left_vertices]);
sources_iv              = sqrt(abs(resampled_data));
sources_iv              = sources_iv/max(sources_iv(:));
fig2 = plotBrainSurface(sources_iv, Vertices10, Faces10);
%% null model
[z_map, nulls] = generateSpatialNulls(sources_iv, 'n_nulls', 1000);
