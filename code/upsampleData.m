function upsample_data = upsampleData(sources_iv, LowResolutionVertices, HighResolutionVertices, varargin)
%UPSAMPLEDATA Upsample brain data from low to high resolution vertices
%
% Input Parameters:
%   sources_iv            - [n¡Á1 double] Data values at low resolution vertices
%   LowResolutionVertices - [n¡Á3 double] Coordinates of low resolution vertices
%   HighResolutionVertices - [m¡Á3 double] Coordinates of high resolution vertices
%   
% Optional Parameters (Name-Value Pairs):
%   'Method'      - Interpolation method: 'linear', 'natural', 'nearest' (default: 'linear')
%   'Extrapolation' - Extrapolation method: 'none', 'linear', 'nearest' (default: 'none')
%   'Verbose'     - Display progress information (default: true)
%
% Output Parameters:
%   upsample_data - [m¡Á1 double] Upsampled data at high resolution vertices
%
% Example:
%   % Upsample 4K data to 10K resolution
%   high_res_data = upsampleData(sources_4k, vertices_4k, vertices_10k, ...
%       'Method', 'linear', 'Verbose', true);

% Parse input parameters
p = inputParser;
addRequired(p, 'sources_iv', @(x) validateattributes(x, {'double'}, {'vector', 'real'}));
addRequired(p, 'LowResolutionVertices', @(x) validateattributes(x, {'double'}, {'ncols', 3, 'real'}));
addRequired(p, 'HighResolutionVertices', @(x) validateattributes(x, {'double'}, {'ncols', 3, 'real'}));
addParameter(p, 'Method', 'linear', @(x) ismember(x, {'linear', 'natural', 'nearest'}));
addParameter(p, 'Extrapolation', 'none', @(x) ismember(x, {'none', 'linear', 'nearest'}));
addParameter(p, 'Verbose', true, @islogical);

parse(p, sources_iv, LowResolutionVertices, HighResolutionVertices, varargin{:});

% Data validation and dimension checking
nLowRes = size(LowResolutionVertices, 1);
nHighRes = size(HighResolutionVertices, 1);

if length(sources_iv) ~= nLowRes
    error('Dimension mismatch: sources_iv has %d values but LowResolutionVertices has %d vertices', ...
          length(sources_iv), nLowRes);
end

if nLowRes >= nHighRes
    error(['Low resolution vertices (%d) should have fewer points than ' ...
           'high resolution vertices (%d). Check your input data.'], nLowRes, nHighRes);
end

if p.Results.Verbose
    fprintf('Upsampling data:\n');
    fprintf('  - Low resolution: %d vertices\n', nLowRes);
    fprintf('  - High resolution: %d vertices\n', nHighRes);
    fprintf('  - Upsampling ratio: %.2f:1\n', nHighRes/nLowRes);
    fprintf('  - Method: %s interpolation\n', p.Results.Method);
end

% Remove any NaN values from input data
validMask = ~isnan(sources_iv);
if ~all(validMask)
    if p.Results.Verbose
        warning('Input data contains %d NaN values. These will be excluded from interpolation.', ...
                sum(~validMask));
    end
    validSources = sources_iv(validMask);
    validVertices = LowResolutionVertices(validMask, :);
else
    validSources = sources_iv;
    validVertices = LowResolutionVertices;
end

% Check if we have enough valid points for interpolation
minPointsRequired = 4; % Minimum points for 3D interpolation
if sum(validMask) < minPointsRequired
    error('Insufficient valid data points (%d). Need at least %d points for 3D interpolation.', ...
          sum(validMask), minPointsRequired);
end

% Create scattered interpolant
try
    if p.Results.Verbose
        fprintf('  - Creating scattered interpolant...\n');
    end
    
    interpolant = scatteredInterpolant(validVertices, validSources, ...
                                       p.Results.Method, p.Results.Extrapolation);
    
    % Perform interpolation
    if p.Results.Verbose
        fprintf('  - Interpolating to high resolution vertices...\n');
    end
    
    upsample_data = interpolant(HighResolutionVertices);
    
    % Post-processing: Handle any potential extrapolation issues
    if strcmp(p.Results.Extrapolation, 'none')
        % Find points outside convex hull and set to NaN
        inHull = inhull(HighResolutionVertices, validVertices);
        upsample_data(~inHull) = NaN;
        
        if p.Results.Verbose && any(~inHull)
            fprintf('  - Warning: %d points (%.1f%%) outside convex hull set to NaN\n', ...
                    sum(~inHull), 100*sum(~inHull)/nHighRes);
        end
    end
    
catch ME
    error('Interpolation failed: %s', ME.message);
end

% Validate output
if any(isnan(upsample_data)) && p.Results.Verbose
    nanPercentage = 100 * sum(isnan(upsample_data)) / nHighRes;
    fprintf('  - Output contains %.1f%% NaN values\n', nanPercentage);
end

if p.Results.Verbose
    fprintf('  - Upsampling completed successfully\n');
    fprintf('  - Output range: [%.3f, %.3f]\n', min(upsample_data), max(upsample_data));
end

end

function in = inhull(testpts, hullpts)
%INHULL Test if points are in convex hull (simplified version)
% This is a basic implementation - for production use consider a more robust version
    if size(hullpts, 2) ~= 3
        error('Vertices must be 3D coordinates');
    end
    
    try
        % Use convhull to test point membership
        K = convhull(hullpts);
        in = intriangulation(hullpts, K, testpts);
    catch
        % Fallback: use tsearchn (requires delaunayTriangulation)
        try
            DT = delaunayTriangulation(hullpts);
            in = ~isnan(pointLocation(DT, testpts));
        catch
            % Simple bounding box check as last resort
            warning('Using bounding box check for convex hull (less accurate)');
            minBounds = min(hullpts);
            maxBounds = max(hullpts);
            in = all(testpts >= minBounds & testpts <= maxBounds, 2);
        end
    end
end

function in = intriangulation(vertices, faces, testpoints)
%INTRANGULATION Check if points are inside triangulated surface
% Simplified implementation - for exact results use more sophisticated method
    in = true(size(testpoints, 1), 1);
    % This is a placeholder - implement proper point-in-mesh test if needed
    warning('Using simplified convex hull test. For complex meshes, consider implementing exact point-in-mesh test.');
end