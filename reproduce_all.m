function reproduce_all(varargin)
%REPRODUCE_ALL Reproduce numerical experiments and figures from the paper.
%   REPRODUCE_ALL() runs all numerical experiments (~1 minute), saves
%   data/opt2_experiments.mat, and generates figures/opt2_experiments.pdf.
%
%   Optional name-value pairs:
%     'mode', M   - 'compute&plot' (default), 'compute', or 'plot'
%
%   Examples:
%     reproduce_all()                    % Compute experiments and plot figure
%     reproduce_all('mode', 'compute')   % Headless server: compute and save .mat
%     reproduce_all('mode', 'plot')      % Plot from saved data/*.mat
%
%   See also: SETPATH, RUN_OPT2, FIGURES_OPT2

setpath();

p = inputParser;
addParameter(p, 'mode', 'compute&plot', ...
    @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
parse(p, varargin{:});
mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});

root = fileparts(mfilename('fullpath'));
mat_file = fullfile(root, 'data', 'opt2_experiments.mat');

fprintf('=========================================================================\n');
fprintf('  Reproducing results for:\n');
fprintf('  "Selecting 2 x 2 submatrices with optimally bounded inverses"\n');
fprintf('=========================================================================\n');
fprintf('  mode = %s\n\n', mode);

if strcmp(mode, 'plot') && ~isfile(mat_file)
    fprintf('Saved data not found at %s; falling back to compute&plot.\n\n', mat_file);
    mode = 'compute&plot';
end

if strcmp(mode, 'compute') || strcmp(mode, 'compute&plot')
    run_opt2();
end

if strcmp(mode, 'plot') || strcmp(mode, 'compute&plot')
    figures_opt2();
end
end
