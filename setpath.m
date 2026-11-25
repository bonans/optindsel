function setpath()
%SETPATH Add optindsel source and experiment directories to MATLAB path.
%   SETPATH() configures the MATLAB path for all algorithms and
%   reproducibility scripts in this repository. No external toolboxes are
%   required.
%
%   See also: REPRODUCE_ALL, RUN_OPT2, FIGURES_OPT2

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'src'));
addpath(fullfile(this_dir, 'experiments'));
end
