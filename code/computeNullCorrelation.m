function [corr_actual, corr_nulls, p_value, stats] = computeNullCorrelation(nulls, data, varargin)
%COMPUTENULLCORRELATION Compute correlation between data and null models
%
% This function computes the correlation between actual data and null models,
% providing both the actual correlation and the null distribution of correlations.
%
% Input Parameters:
%   nulls - Null models matrix [vertices ¡Á permutations] from generateSpatialNulls
%   data  - Actual data: can be numeric vector, GIfTI object, GIfTI file path,
%           or cell array of two files for left/right hemispheres
%
% Optional Parameters (Name-Value Pairs):
%   'corr_type'    - Correlation type: 'pearson', 'spearman' (default: 'pearson')
%   'hemi'         - Hemisphere: 'lh', 'rh', or 'both' (default: 'both')
%   'verbose'      - Display progress information (default: true)
%   'tail'         - P-value tail: 'both', 'right', 'left' (default: 'both')
%
% Output Parameters:
%   corr_actual - Actual correlation coefficient between data and null mean
%   corr_nulls  - Correlation coefficients for each null model [1 ¡Á permutations]
%   p_value     - P-value for the actual correlation
%   stats       - Additional statistics structure
%
% Example:
%   % Basic usage with vector data
%   [corr_actual, corr_nulls, p_value] = computeNullCorrelation(nulls, data_vector);
%
%   % Using GIfTI files for both hemispheres
%   [corr_actual, corr_nulls, p_value] = computeNullCorrelation(nulls, ...
%       {'lh_data.func.gii', 'rh_data.func.gii'}, 'corr_type', 'spearman');

% Parse input parameters
p = inputParser;
addRequired(p, 'nulls', @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(p, 'data', @(x) isnumeric(x) || ischar(x) || isa(x, 'gifti') || iscell(x));
addParameter(p, 'corr_type', 'pearson', @(x) ismember(x, {'pearson', 'spearman'}));
addParameter(p, 'hemi', 'both', @(x) ismember(x, {'lh', 'rh', 'both'}));
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'tail', 'both', @(x) ismember(x, {'both', 'right', 'left'}));

parse(p, nulls, data, varargin{:});

if p.Results.verbose
    fprintf('Computing correlation with null models...\n');
    fprintf('  Null models: %d vertices ¡Á %d permutations\n', ...
        size(nulls, 1), size(nulls, 2));
end

% Load and validate data
actual_data = loadAndValidateData(data, p.Results.hemi, size(nulls, 1), p.Results.verbose);

% Compute null mean (representing the null expectation)
null_mean = mean(nulls, 2);

if p.Results.verbose
    fprintf('  Data validation passed: %d vertices\n', length(actual_data));
    fprintf('  Correlation type: %s\n', p.Results.corr_type);
end

% Compute actual correlation
corr_actual = computeCorrelation(actual_data, null_mean, p.Results.corr_type);

% Compute correlations for each null model
n_nulls = size(nulls, 2);
corr_nulls = zeros(1, n_nulls);

if p.Results.verbose
    fprintf('  Computing correlations for %d null models...\n', n_nulls);
end

for i = 1:n_nulls
    if p.Results.verbose && mod(i, 100) == 0
        fprintf('    Processing null %d/%d\n', i, n_nulls);
    end
    corr_nulls(i) = computeCorrelation(actual_data, nulls(:, i), p.Results.corr_type);
end

% Compute p-value
p_value = computePValue(corr_actual, corr_nulls, p.Results.tail);

% Compute additional statistics
stats = computeAdditionalStats(corr_actual, corr_nulls, p_value);

if p.Results.verbose
    fprintf('[Done] Correlation analysis completed\n');
    fprintf('  Actual correlation: %.4f\n', corr_actual);
    fprintf('  Null correlation range: [%.4f, %.4f]\n', min(corr_nulls), max(corr_nulls));
    fprintf('  P-value (%s-tailed): %.4f\n', p.Results.tail, p_value);
    fprintf('  Z-score: %.4f\n', stats.z_score);
end

end

%% Data Loading and Validation Functions

function data = loadAndValidateData(input_data, hemi, expected_vertices, verbose)
%LOADANDVALIDATEDATA Load and validate input data
% Ensures data matches expected vertex count
    
    data = loadData(input_data, hemi, verbose);
    
    % Validate vertex count
    if length(data) ~= expected_vertices
        error(['Vertex count mismatch: Data has %d vertices, ' ...
               'but null models have %d vertices. ' ...
               'Ensure both are in the same space.'], ...
              length(data), expected_vertices);
    end
    
    % Check for NaN values
    nan_count = sum(isnan(data));
    if nan_count > 0
        if verbose
            warning('Data contains %d NaN values (%.1f%%). These will be excluded from correlation.', ...
                    nan_count, 100*nan_count/length(data));
        end
    end
end

function data = loadData(input_data, hemi, verbose)
%LOADDATA Load data from various input types
    
    if isnumeric(input_data)
        % Direct numeric input
        data = input_data(:);  % Ensure column vector
        
    elseif ischar(input_data)
        % Single file path
        if exist(input_data, 'file')
            data = loadGiftiData(input_data, verbose);
        else
            error('File not found: %s', input_data);
        end
        
    elseif iscell(input_data) && numel(input_data) == 2
        % Cell array with left/right hemisphere files
        data = loadBothHemispheres(input_data, hemi, verbose);
        
    elseif isa(input_data, 'gifti')
        % GIfTI object
        data = double(input_data.cdata);
        
    else
        error('Unsupported data type: %s', class(input_data));
    end
    
    % Ensure column vector
    data = data(:);
end

function data = loadGiftiData(filename, verbose)
%LOADGIFTIDATA Load data from GIfTI file
    
    if verbose
        fprintf('    Loading GIfTI file: %s\n', filename);
    end
    
    try
        gii_data = gifti(filename);
        data = double(gii_data.cdata);
    catch ME
        error('Failed to load GIfTI file %s: %s', filename, ME.message);
    end
end

function data = loadBothHemispheres(filepaths, hemi, verbose)
%LOADBOTHHEMISPHERES Load data from left/right hemisphere files
    
    if verbose
        fprintf('    Loading both hemispheres...\n');
    end
    
    % Load left hemisphere
    if exist(filepaths{1}, 'file')
        lh_data = loadGiftiData(filepaths{1}, verbose);
    else
        error('Left hemisphere file not found: %s', filepaths{1});
    end
    
    % Load right hemisphere  
    if exist(filepaths{2}, 'file')
        rh_data = loadGiftiData(filepaths{2}, verbose);
    else
        error('Right hemisphere file not found: %s', filepaths{2});
    end
    
    % Combine based on hemisphere specification
    switch hemi
        case 'both'
            data = [lh_data; rh_data];
        case 'lh'
            data = lh_data;
        case 'rh'
            data = rh_data;
    end
    
    if verbose
        fprintf('    Combined data: %d vertices\n', length(data));
    end
end

%% Correlation Computation Functions

function r = computeCorrelation(x, y, corr_type)
%COMPUTECORRELATION Compute correlation between two vectors
    
    % Remove NaN values (pairwise deletion)
    valid_mask = ~isnan(x) & ~isnan(y);
    x_valid = x(valid_mask);
    y_valid = y(valid_mask);
    
    if length(x_valid) < 2
        warning('Insufficient valid data points for correlation. Returning NaN.');
        r = NaN;
        return;
    end
    
    switch lower(corr_type)
        case 'pearson'
            r = corr(x_valid, y_valid, 'type', 'Pearson');
        case 'spearman'
            r = corr(x_valid, y_valid, 'type', 'Spearman');
        otherwise
            error('Unknown correlation type: %s', corr_type);
    end
end

function p_val = computePValue(corr_actual, corr_nulls, tail)
%COMPUTEPVALUE Compute p-value from null distribution
    
    % Remove any NaN values from null correlations
    valid_nulls = corr_nulls(~isnan(corr_nulls));
    
    if isempty(valid_nulls)
        p_val = NaN;
        return;
    end
    
    switch lower(tail)
        case 'both'
            % Two-tailed p-value
            extreme_count = sum(abs(valid_nulls) >= abs(corr_actual));
            p_val = (extreme_count + 1) / (length(valid_nulls) + 1);
            
        case 'right'
            % Right-tailed p-value (positive correlation)
            extreme_count = sum(valid_nulls >= corr_actual);
            p_val = (extreme_count + 1) / (length(valid_nulls) + 1);
            
        case 'left'
            % Left-tailed p-value (negative correlation)
            extreme_count = sum(valid_nulls <= corr_actual);
            p_val = (extreme_count + 1) / (length(valid_nulls) + 1);
    end
end

function stats = computeAdditionalStats(corr_actual, corr_nulls, p_value)
%COMPUTEADDITIONALSTATS Compute additional statistics
    
    % Remove NaN values
    valid_nulls = corr_nulls(~isnan(corr_nulls));
    
    stats = struct();
    stats.null_mean = mean(valid_nulls);
    stats.null_std = std(valid_nulls);
    stats.null_ci95 = prctile(valid_nulls, [2.5, 97.5]);
    
    % Compute z-score
    if stats.null_std > 0
        stats.z_score = (corr_actual - stats.null_mean) / stats.null_std;
    else
        stats.z_score = NaN;
    end
    
    stats.effect_size = corr_actual - stats.null_mean;
    stats.null_min = min(valid_nulls);
    stats.null_max = max(valid_nulls);
    stats.p_value = p_value;
end