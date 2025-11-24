function [fig, ax, patchHandle] = plotBrainSurface(sources_iv, Vertices, Faces, varargin)
%PLOTBRAINSURFACE Plot 3D brain surface map
%
% Input Parameters:
%   sources_iv    - [n¡Á1 double] Values for each vertex (e.g., activation intensity, connectivity strength)
%   Vertices      - [n¡Á3 double] Vertex coordinate matrix
%   Faces         - [m¡Á3 double] Triangle face connectivity matrix
%   
% Optional Parameters (Name-Value Pairs):
%   'Colormap'    - Color map (default: built-in brain atlas colormap)
%   'ViewAngle'   - Viewing angle [azimuth, elevation] (default: [0, 0])
%   'Alpha'       - Transparency (default: 1)
%   'FigureName'  - Figure window name (default: 'Brain Surface Visualization')
%   'LogScale'    - Whether to apply logarithmic scaling to data (default: true)
%   'Parent'      - Parent axes handle (default: create new axes)
%
% Output Parameters:
%   fig         - Figure handle
%   ax          - Axes handle  
%   patchHandle - Patch object handle
%
% Example:
%   [fig, ax, patch] = plotBrainSurface(activation, vertices, faces, ...
%       'ViewAngle', [90, 0], 'Alpha', 0.8);

% Parse input parameters
p = inputParser;
addRequired(p, 'sources_iv', @(x) validateattributes(x, {'double'}, {'vector'}));
addRequired(p, 'Vertices', @(x) validateattributes(x, {'double'}, {'ncols', 3}));
addRequired(p, 'Faces', @(x) validateattributes(x, {'double'}, {'ncols', 3}));
addParameter(p, 'Colormap', [], @(x) isempty(x) || ismatrix(x));
addParameter(p, 'ViewAngle', [0, 0], @(x) validateattributes(x, {'double'}, {'numel', 2}));
addParameter(p, 'Alpha', 1, @(x) validateattributes(x, {'double'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(p, 'FigureName', 'Brain Surface Visualization', @ischar);
addParameter(p, 'LogScale', true, @islogical);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x, 'axes'));

parse(p, sources_iv, Vertices, Faces, varargin{:});

% Data validation
nVertices = size(Vertices, 1);
if length(sources_iv) ~= nVertices
    error('Input data dimension mismatch: sources_iv(%d) vs Vertices(%d)', ...
          length(sources_iv), nVertices);
end

% Create figure and axes
fig = figure('Name', p.Results.FigureName, ...
             'Color', 'white', ...
             'NumberTitle', 'off');

if isempty(p.Results.Parent)
    ax = axes(fig);
else
    ax = p.Results.Parent;
    fig = ax.Parent;
end

% Load default colormap
if isempty(p.Results.Colormap)
    try
        colorMapData = load('matlab/data/plot/mycolormap_brain_basic_conn.mat');
        cmap = colorMapData.cmap_a;
    catch ME
        warning('Cannot load default colormap, using jet colormap: %s', ME.message);
        cmap = jet(256);
    end
else
    cmap = p.Results.Colormap;
end

% Data processing
if p.Results.LogScale
    plotData = log(1 + sources_iv);
else
    plotData = sources_iv;
end

% Create brain surface patch
patchHandle = patch(ax, ...
    'Faces', Faces, ...
    'Vertices', Vertices, ...
    'FaceVertexCData', plotData, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', p.Results.Alpha, ...
    'AlphaDataMapping', 'none', ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud', ...
    'EdgeLighting', 'gouraud');

% Set axes properties
set(ax, 'Color', 'white', ...
        'xcolor', 'white', ...
        'ycolor', 'white', ...
        'zcolor', 'white');
    
% Set view and colormap
view(ax, p.Results.ViewAngle);
colormap(ax, cmap);
colorbar(ax);

% Add lighting
light(ax, 'Position', [1 1 1], 'Style', 'infinite');
light(ax, 'Position', [-1 -1 -1], 'Style', 'infinite');

% Enable interaction
rotate3d(ax, 'on');

% Set equal axis scaling for proper proportions
axis(ax, 'equal');
axis(ax, 'tight');

% Add grid and labels
grid(ax, 'on');
xlabel(ax, 'X (mm)');
ylabel(ax, 'Y (mm)'); 
zlabel(ax, 'Y (mm)');

% Set figure title
title(ax, 'Brain Surface Map', 'Color', 'white', 'FontSize', 12);

fprintf('Brain surface plot created successfully:\n');
fprintf('  - Number of vertices: %d\n', nVertices);
fprintf('  - Number of faces: %d\n', size(Faces, 1));
fprintf('  - Data range: [%.3f, %.3f]\n', min(sources_iv), max(sources_iv));

end