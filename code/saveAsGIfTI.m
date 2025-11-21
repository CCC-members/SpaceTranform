function saveAsGIfTI(data, output_path, output_filename, varargin)
%SAVEASGIFTI Save data as GIfTI format files
%
% This function saves brain surface data as GIfTI files, automatically handling
% both single hemisphere and both hemispheres cases.
%
% Input Parameters:
%   data            - Input data vector (e.g., N¡Á1 for brain surface data)
%   output_path     - Directory where files will be saved
%   output_filename - Base filename (without extension)
%
% Optional Parameters (Name-Value Pairs):
%   'hemi'          - Hemisphere: 'lh', 'rh', or 'both' (default: 'both')
%   'split_point'   - Custom split point for both hemispheres (default: auto)
%   'data_type'     - Data type for GIfTI file: 'functional', 'shape', 'label' 
%                     (default: 'functional')
%   'verbose'       - Display progress information (default: true)
%   'overwrite'     - Overwrite existing files (default: false)
%
% Output:
%   Saves one or two GIfTI files depending on hemisphere specification
%
% Example:
%   % Save both hemispheres with automatic splitting
%   saveAsGIfTI(z_map, './results', 'z_score_map', 'hemi', 'both');
%
%   % Save only left hemisphere
%   saveAsGIfTI(data_lh, './results', 'left_hemi_data', 'hemi', 'lh');

% Parse input parameters
p = inputParser;
addRequired(p, 'data', @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(p, 'output_path', @ischar);
addRequired(p, 'output_filename', @ischar);
addParameter(p, 'hemi', 'both', @(x) ismember(x, {'lh', 'rh', 'both'}));
addParameter(p, 'split_point', 'auto', @(x) ischar(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'data_type', 'functional', @ischar);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'overwrite', false, @islogical);

parse(p, data, output_path, output_filename, varargin{:});

% Validate output directory
if ~exist(output_path, 'dir')
    if p.Results.verbose
        fprintf('Creating output directory: %s\n', output_path);
    end
    mkdir(output_path);
end

% Remove file extension if provided
[~, filename_base, ~] = fileparts(output_filename);

% Process data based on hemisphere
switch lower(p.Results.hemi)
    case 'both'
        saveBothHemispheres(data, output_path, filename_base, p.Results);
    case {'lh', 'rh'}
        saveSingleHemisphere(data, output_path, filename_base, p.Results.hemi, p.Results);
end

if p.Results.verbose
    fprintf('GIfTI files saved successfully to: %s\n', output_path);
end
end

%% Helper Functions

function saveBothHemispheres(data, output_path, filename_base, params)
%SAVEBOTHHEMISPHERES Save data for both hemispheres
    
    n_vertices = length(data);
    
    % Determine split point
    if strcmp(params.split_point, 'auto')
        % Auto-split: assume equal division
        split_point = floor(n_vertices / 2);
    else
        % Use custom split point
        split_point = params.split_point;
        if split_point > n_vertices
            error('Split point (%d) exceeds data length (%d)', split_point, n_vertices);
        end
    end
    
    % Split data into left and right hemispheres
    data_lh = data(1:split_point);
    data_rh = data(split_point+1:end);
    
    if params.verbose
        fprintf('Saving both hemispheres:\n');
        fprintf('  Total vertices: %d\n', n_vertices);
        fprintf('  Left hemisphere: %d vertices\n', length(data_lh));
        fprintf('  Right hemisphere: %d vertices\n', length(data_rh));
        fprintf('  Split point: %d\n', split_point);
    end
    
    % Save left hemisphere
    lh_filename = fullfile(output_path, [filename_base '_lh.func.gii']);
    saveGIfTIFile(data_lh, lh_filename, params);
    
    % Save right hemisphere
    rh_filename = fullfile(output_path, [filename_base '_rh.func.gii']);
    saveGIfTIFile(data_rh, rh_filename, params);
end

function saveSingleHemisphere(data, output_path, filename_base, hemi, params)
%SAVESINGLEHEMISPHERE Save data for single hemisphere
    
    if params.verbose
        fprintf('Saving %s hemisphere: %d vertices\n', hemi, length(data));
    end
    
    % Create filename with hemisphere suffix
    filename = fullfile(output_path, [filename_base '_' hemi '.func.gii']);
    saveGIfTIFile(data, filename, params);
end

function saveGIfTIFile(data, filename, params)
%SAVEGIFTIFILE Save single data vector as GIfTI file
    
    % Check if file exists and handle overwrite
    if exist(filename, 'file') && ~params.overwrite
        if params.verbose
            warning('File already exists and overwrite is false: %s', filename);
        end
        return;
    end
    
    try
        % Create gifti object
        g = gifti();
        
        % Set data based on data type
        switch lower(params.data_type)
            case 'functional'
                g.cdata = data;
            case 'shape'
                % For shape data, we might need to handle 3D coordinates
                if size(data, 2) == 3
                    g.vertices = data;
                else
                    warning('Shape data should be N¡Á3. Using cdata instead.');
                    g.cdata = data;
                end
            case 'label'
                % For label data, use integers
                g.cdata = int32(data);
            otherwise
                g.cdata = data;
        end
        
        % Save the gifti file
        save(g, filename);
        
        if params.verbose
            fprintf('  Saved: %s\n', filename);
        end
        
    catch ME
        error('Failed to save GIfTI file %s: %s', filename, ME.message);
    end
end

%% Additional utility function for separate hemisphere data

function saveSeparateGIfTI(data_lh, data_rh, output_path, filename_base, varargin)
%SAVESEPARATEGIFTI Save separate left and right hemisphere data as GIfTI files
%
% This function is useful when you already have separated hemisphere data
%
% Input Parameters:
%   data_lh        - Left hemisphere data vector
%   data_rh        - Right hemisphere data vector
%   output_path    - Directory where files will be saved
%   filename_base  - Base filename (without extension)
%
% Optional Parameters: Same as saveAsGIfTI

% Parse input parameters
p = inputParser;
addRequired(p, 'data_lh', @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(p, 'data_rh', @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(p, 'output_path', @ischar);
addRequired(p, 'filename_base', @ischar);
addParameter(p, 'data_type', 'functional', @ischar);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'overwrite', false, @islogical);

parse(p, data_lh, data_rh, output_path, filename_base, varargin{:});

% Validate output directory
if ~exist(output_path, 'dir')
    if p.Results.verbose
        fprintf('Creating output directory: %s\n', output_path);
    end
    mkdir(output_path);
end

if p.Results.verbose
    fprintf('Saving separate hemisphere files:\n');
    fprintf('  Left hemisphere: %d vertices\n', length(data_lh));
    fprintf('  Right hemisphere: %d vertices\n', length(data_rh));
    fprintf('  Total vertices: %d\n', length(data_lh) + length(data_rh));
end

% Save left hemisphere
lh_filename = fullfile(output_path, [filename_base '_lh.func.gii']);
saveGIfTIFile(data_lh, lh_filename, p.Results);

% Save right hemisphere
rh_filename = fullfile(output_path, [filename_base '_rh.func.gii']);
saveGIfTIFile(data_rh, rh_filename, p.Results);

if p.Results.verbose
    fprintf('Separate hemisphere GIfTI files saved successfully to: %s\n', output_path);
end
end