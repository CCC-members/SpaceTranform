function resampled_data = transformationBrainData(source_iv, source_sphere, target_sphere, varargin)
%RESAMPLEBRAINDATA Resample brain data between different coordinate spaces
%
% Input Parameters:
%   source_iv     - [n¡Á1 double] Data values at source vertices
%   source_sphere - [struct/cell] Source space sphere information
%                   Can be: file paths (string), gifti objects, or struct containing these
%   target_sphere - [struct/cell] Target space sphere information
%                   Same format as source_sphere
%
% Optional Parameters (Name-Value Pairs):
%   'Hemi'         - Hemisphere: 'lh', 'rh', or 'both' (default: 'both')
%   'Method'       - Resampling method: 'BARYCENTRIC', 'ADAP_BARY_AREA' (default: 'BARYCENTRIC')
%   'TempDir'      - Temporary directory for intermediate files (default: system temp)
%   'Verbose'      - Display progress information (default: true)
%   'Cleanup'      - Delete temporary files (default: true)
%   'SplitRatio'   - For 'both' mode: [lh_ratio, rh_ratio] or 'auto' (default: 'auto')
%
% Output Parameters:
%   resampled_data - Resampled data vector in target space

% Parse input parameters
p = inputParser;
addRequired(p, 'source_iv', @(x) validateattributes(x, {'double'}, {'vector', 'real'}));
addRequired(p, 'source_sphere', @(x) isstruct(x) || iscell(x) || ischar(x) || isa(x, 'gifti'));
addRequired(p, 'target_sphere', @(x) isstruct(x) || iscell(x) || ischar(x) || isa(x, 'gifti'));
addParameter(p, 'Hemi', 'both', @(x) ismember(x, {'lh', 'rh', 'both'}));
addParameter(p, 'Method', 'BARYCENTRIC', @ischar);
addParameter(p, 'TempDir', tempdir, @ischar);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'Cleanup', true, @islogical);
addParameter(p, 'SplitRatio', 'auto', @(x) ischar(x) || (isnumeric(x) && numel(x) == 2));

parse(p, source_iv, source_sphere, target_sphere, varargin{:});

% Validate Connectome Workbench installation
if ~isWorkbenchInstalled()
    error(['Connectome Workbench not found or not properly configured. '...
           'Please download from: https://www.humanconnectome.org/software/get-connectome-workbench '...
           'and ensure wb_command is in your system PATH.']);
end

% Normalize sphere inputs to struct format with file paths
source_sphere_norm = normalizeSphereInput(source_sphere, p.Results.Hemi, 'source');
target_sphere_norm = normalizeSphereInput(target_sphere, p.Results.Hemi, 'target');

% Create temporary directory
temp_dir = p.Results.TempDir;
if ~exist(temp_dir, 'dir')
    mkdir(temp_dir);
end

try
    % Process based on hemisphere
    switch p.Results.Hemi
        case 'both'
            resampled_data = processBothHemispheres(source_iv, source_sphere_norm, ...
                target_sphere_norm, temp_dir, p.Results);
        case {'lh', 'rh'}
            resampled_data = processSingleHemisphere(source_iv, source_sphere_norm.(p.Results.Hemi), ...
                target_sphere_norm.(p.Results.Hemi), p.Results.Hemi, temp_dir, p.Results);
    end
    
    % Cleanup temporary files
    if p.Results.Cleanup
        cleanupTempFiles(temp_dir);
    end
    
    if p.Results.Verbose
        fprintf('Spatial resampling completed successfully.\n');
        fprintf('  Input data size: %d vertices\n', length(source_iv));
        fprintf('  Output data size: %d vertices\n', length(resampled_data));
    end
    
catch ME
    % Ensure cleanup even if error occurs
    if p.Results.Cleanup
        cleanupTempFiles(temp_dir);
    end
    rethrow(ME);
end

end

%% Helper Functions

function isInstalled = isWorkbenchInstalled()
%ISWORKBENCHINSTALLED Check if Connectome Workbench is properly installed
    try
        [status, ~] = system('wb_command -version');
        isInstalled = (status == 0);
    catch
        isInstalled = false;
    end
end

function sphere_struct = normalizeSphereInput(sphere_input, hemi, space_type)
%NORMALIZESPHEREINPUT Convert various sphere inputs to standardized struct with file paths
    sphere_struct = struct();
    
    if isstruct(sphere_input)
        % Handle struct input - extract file paths from gifti objects or use strings directly
        fields = fieldnames(sphere_input);
        for i = 1:length(fields)
            field = fields{i};
            value = sphere_input.(field);
            sphere_struct.(field) = extractFilePath(value, field);
        end
    elseif iscell(sphere_input) && numel(sphere_input) == 2
        % Handle cell array input
        sphere_struct.lh = extractFilePath(sphere_input{1}, 'lh');
        sphere_struct.rh = extractFilePath(sphere_input{2}, 'rh');
    elseif ischar(sphere_input) || isa(sphere_input, 'gifti')
        % Handle single path or gifti object
        switch hemi
            case 'both'
                error(['For ''both'' hemisphere mode, provide sphere information as struct '...
                       'with ''lh'' and ''rh'' fields or as cell array {lh_info, rh_info}']);
            case {'lh', 'rh'}
                sphere_struct.(hemi) = extractFilePath(sphere_input, hemi);
        end
    else
        error('Invalid sphere input format for %s space', space_type);
    end
    
    % Validate that required hemispheres are present
    if strcmp(hemi, 'both')
        if ~isfield(sphere_struct, 'lh') || ~isfield(sphere_struct, 'rh')
            error('Both hemisphere mode requires ''lh'' and ''rh'' sphere information');
        end
    else
        if ~isfield(sphere_struct, hemi)
            error('Hemisphere ''%s'' not found in %s sphere input', hemi, space_type);
        end
    end
end

function file_path = extractFilePath(input_value, hemi)
%EXTRACTFILEPATH Extract file path from various input types
    if ischar(input_value)
        file_path = input_value;
    elseif isa(input_value, 'gifti')
        % For gifti objects, use their original filename if available
        if isprop(input_value, 'filename') && ~isempty(input_value.filename) && exist(input_value.filename, 'file')
            file_path = input_value.filename;
        else
            % Create temporary file for gifti object
            temp_dir = tempdir;
            timestamp = datestr(now, 'HHMMSSFFF');
            temp_file = fullfile(temp_dir, sprintf('temp_%s_sphere_%s.surf.gii', hemi, timestamp));
            save(input_value, temp_file);
            file_path = temp_file;
        end
    else
        error('Unsupported input type: %s', class(input_value));
    end
end

function result = processSingleHemisphere(source_iv, source_sphere_path, ...
    target_sphere_path, hemi, temp_dir, params)
%PROCESSSINGLEHEMISPHERE Process data for a single hemisphere
    
    % Validate sphere files exist and get vertex counts
    [n_vertices_source, source_sphere_gii] = getSurfaceInfo(source_sphere_path, 'source', hemi);
    [n_vertices_target, target_sphere_gii] = getSurfaceInfo(target_sphere_path, 'target', hemi);
    
    % Validate input data matches source sphere vertex count
    if length(source_iv) ~= n_vertices_source
        error(['Vertex count mismatch: Input data has %d vertices, but source sphere has %d vertices. '...
               'Please ensure your input data matches the source space.'], ...
              length(source_iv), n_vertices_source);
    end
    
    % Create temporary input and output files
    timestamp = datestr(now, 'HHMMSSFFF');
    temp_input = fullfile(temp_dir, sprintf('temp_input_%s_%s.func.gii', hemi, timestamp));
    temp_output = fullfile(temp_dir, sprintf('temp_output_%s_%s.func.gii', hemi, timestamp));
    
    % Save input data as GIFTI
    saveGiftiData(source_iv, temp_input);
    
    % Build and execute resampling command
    cmd = buildResampleCommand(temp_input, source_sphere_path, target_sphere_path, ...
        temp_output, params.Method);
    
    if params.Verbose
        fprintf('Resampling %s hemisphere...\n', hemi);
        fprintf('  Source sphere: %s (%d vertices)\n', source_sphere_path, n_vertices_source);
        fprintf('  Target sphere: %s (%d vertices)\n', target_sphere_path, n_vertices_target);
        fprintf('  Input data: %d vertices\n', length(source_iv));
        fprintf('  Method: %s\n', params.Method);
        fprintf('  Command: %s\n', cmd);
    end
    
    [status, result_msg] = system(cmd);
    
    if status ~= 0
        error('Resampling failed for %s hemisphere. System command error: %s', hemi, result_msg);
    end
    
    % Load resampled data
    if ~exist(temp_output, 'file')
        error('Output file not created: %s', temp_output);
    end
    
    gifti_data = gifti(temp_output);
    result = double(gifti_data.cdata);
    
    % Validate output vertex count
    if length(result) ~= n_vertices_target
        warning('Output vertex count (%d) does not match target sphere (%d)', ...
                length(result), n_vertices_target);
    end
end

function [n_vertices, sphere_gii] = getSurfaceInfo(sphere_path, space_type, hemi)
%GETSURFACEINFO Get vertex count and gifti object for a surface file
    if ~exist(sphere_path, 'file')
        error('%s sphere file for %s hemisphere not found: %s', space_type, hemi, sphere_path);
    end
    
    try
        sphere_gii = gifti(sphere_path);
        n_vertices = size(sphere_gii.vertices, 1);
    catch ME
        error('Failed to read %s sphere file for %s hemisphere: %s', space_type, hemi, ME.message);
    end
end

function result = processBothHemispheres(source_iv, source_sphere, ...
    target_sphere, temp_dir, params)
%PROCESSBOTHHEMISPHERES Process data for both hemispheres (DEFAULT MODE)
    
    % Get vertex counts for source spheres to help with splitting
    [n_vertices_lh_source, ~] = getSurfaceInfo(source_sphere.lh, 'source', 'lh');
    [n_vertices_rh_source, ~] = getSurfaceInfo(source_sphere.rh, 'source', 'rh');
    total_source_vertices = n_vertices_lh_source + n_vertices_rh_source;
    
    % Determine split ratio for left/right hemispheres
    if strcmp(params.SplitRatio, 'auto')
        % Use actual vertex counts from source spheres for splitting
        if length(source_iv) == total_source_vertices
            % If total matches, use the actual sphere vertex counts
            source_lh = source_iv(1:n_vertices_lh_source);
            source_rh = source_iv(n_vertices_lh_source+1:end);
        else
            % Fallback: equal division
            half_point = floor(length(source_iv) / 2);
            source_lh = source_iv(1:half_point);
            source_rh = source_iv(half_point+1:end);
            if params.Verbose
                warning('Using equal split. Input data (%d) does not match total source vertices (%d)', ...
                        length(source_iv), total_source_vertices);
            end
        end
    else
        % Use provided split ratio
        ratio = params.SplitRatio;
        split_point = floor(length(source_iv) * ratio(1));
        source_lh = source_iv(1:split_point);
        source_rh = source_iv(split_point+1:end);
    end
    
    if params.Verbose
        fprintf('Processing both hemispheres (default mode):\n');
        fprintf('  Left hemisphere: %d vertices (input), %d vertices (source sphere)\n', ...
                length(source_lh), n_vertices_lh_source);
        fprintf('  Right hemisphere: %d vertices (input), %d vertices (source sphere)\n', ...
                length(source_rh), n_vertices_rh_source);
        fprintf('  Source spheres: LH=%s, RH=%s\n', ...
            source_sphere.lh, source_sphere.rh);
        fprintf('  Target spheres: LH=%s, RH=%s\n', ...
            target_sphere.lh, target_sphere.rh);
    end
    
    % Process left hemisphere
    resampled_lh = processSingleHemisphere(source_lh, source_sphere.lh, ...
        target_sphere.lh, 'lh', temp_dir, params);
    
    % Process right hemisphere  
    resampled_rh = processSingleHemisphere(source_rh, source_sphere.rh, ...
        target_sphere.rh, 'rh', temp_dir, params);
    
    % Combine results
    result = [resampled_lh; resampled_rh];
end

function cmd = buildResampleCommand(input_file, source_sphere_path, ...
    target_sphere_path, output_file, method)
%BUILDRESAMPLECOMMAND Build wb_command resampling command
    cmd = sprintf('wb_command -metric-resample "%s" "%s" "%s" %s "%s"', ...
        input_file, source_sphere_path, target_sphere_path, method, output_file);
end

function saveGiftiData(data, filename)
%SAVEGIFTIDATA Save data as GIFTI file
    g = gifti();
    g.cdata = data;
    save(g, filename);
end

function cleanupTempFiles(temp_dir)
%CLEANUPTEMPFILES Remove temporary files
    temp_files = dir(fullfile(temp_dir, 'temp_*.gii'));
    for i = 1:length(temp_files)
        try
            delete(fullfile(temp_dir, temp_files(i).name));
        catch
            % Skip if file cannot be deleted
        end
    end
end