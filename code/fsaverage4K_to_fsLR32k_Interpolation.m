function metric_interp = fsaverage4K_to_fsLR32k_Interpolation(J, Vertices_4K, Vertices_32K)
% Interpolate 4K brain data to 32K surface mesh
%
% Input Parameters:
%   J: Metric values at each vertex of 4K mesh, vector  
%   Vertices_4K: 3D vertex coordinates of 4K surface, 
%   Vertices_32K: 3D vertex coordinates of 32K surface, 
%
% Output Parameters:
%   metric_interp: Interpolated metric values on 32K surface, vector size matching rows of Vertices_32K

% Create scattered data interpolant object
% scatteredInterpolant builds interpolation function based on 4K vertex coordinates and corresponding metric values
% Uses linear interpolation method by default
interp = scatteredInterpolant(Vertices_4K, J);

% Apply interpolation function to 32K surface vertex coordinates
% Computes interpolated metric values for each point on the 32K mesh
metric_interp = interp(Vertices_32K);

% Calculate threshold: Use 25th percentile as cutoff value
% This threshold filters out lower values, keeping the top 75% of the distribution
threshold = prctile(metric_interp, 25);

% Apply thresholding: Set data points below threshold to zero
% Implemented using logical indexing: metric_interp > threshold creates logical mask
% Only points above threshold retain their original values, others become 0
metric_interp = metric_interp .* (metric_interp > threshold);

end