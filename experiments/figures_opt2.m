function figures_opt2()
%figures_opt2 - Generate MATLAB reproduction of Figure 1 from saved data
%   Reads data/opt2_experiments.mat and writes figures/opt2_experiments.pdf.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
data_dir = fullfile(root, 'data');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

S = load(fullfile(data_dir, 'opt2_experiments.mat'));

% Shared method styling (col index, label, color, marker, linestyle)
method_styles = {
    1, 'Algorithm 2',   [0.00 0.38 0.65], 'o', '-'
    2, 'exact maxvol',  [0.85 0.33 0.10], 's', '-'
    3, 'maxvol',        [0.80 0.18 0.42], 'p', '-.'
    4, 'QRCP',          [0.55 0.33 0.65], 'v', '-'
    5, 'BSQR',          [0.90 0.60 0.00], 'd', '--'
    6, 'LUPP',          [0.25 0.58 0.30], '^', '--'
    7, 'exact minimum', [0.12 0.12 0.12], '+', ':'
    };

W = 5.9; H = 3.20;
fig = new_figure(W, H);

%% --- Left panel: Perturbed three-block matrices vs n ---
ax1 = panel(fig, [0.66, 1.12, 2.10, 1.74], W, H);
yline(ax1, 1.0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
for k = [7, 6, 4, 5, 3, 2, 1]
    idx_c = method_styles{k, 1};
    plot(ax1, S.ns_cond, S.cond_max(:, idx_c), 'LineStyle', method_styles{k, 5}, ...
        'Marker', method_styles{k, 4}, 'Color', method_styles{k, 3}, ...
        'LineWidth', 1.25, 'MarkerSize', 4.2);
end
ymax_c = min(1.85, ceil(20 * max(S.cond_max(:))) / 20 + 0.04);
set(ax1, 'XScale', 'log', 'XLim', [8 / 1.18, 1024 * 1.18], ...
    'XTick', [8, 16, 64, 256, 1024], ...
    'XTickLabel', {'$8$', '$16$', '$64$', '$256$', '$1024$'}, ...
    'YTick', [1.0, 1.2, 1.4, 1.6], ...
    'YTickLabel', {'$1.0$', '$1.2$', '$1.4$', '$1.6$'}, ...
    'TickLabelInterpreter', 'latex', ...
    'XTickLabelRotation', 0, ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylim(ax1, [0.88, ymax_c]);
xl1 = xlabel(ax1, 'number of rows $n$', 'Interpreter', 'latex', 'FontSize', 9.5);
xl1.Units = 'normalized';
xl1.Position = [0.5, -0.125, 0];
xl1.HorizontalAlignment = 'center';
yl1 = ylabel(ax1, '$\max\,\|U_I^{-1}\|_2 / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 10);
yl1.Units = 'normalized';
yl1.Position = [-0.18, 0.5, 0];
yl1.HorizontalAlignment = 'center';
t1 = title(ax1, 'Perturbed 3-block matrices (max)', 'Interpreter', 'latex', ...
    'FontWeight', 'normal', 'FontSize', 10);
t1.Units = 'normalized';
t1.Position = [0.5, 1.045, 0];
t1.HorizontalAlignment = 'center';
t1.VerticalAlignment = 'baseline';

%% --- Right panel: Perturbed Proposition 4.2 transition across sigma ---
ax2 = panel(fig, [3.60, 1.12, 2.10, 1.74], W, H);
yline(ax2, 1.0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
for k = [7, 6, 4, 5, 3, 2, 1]
    idx_c = method_styles{k, 1};
    plot(ax2, S.sigmas, S.trans_med(:, idx_c), 'LineStyle', method_styles{k, 5}, ...
        'Marker', method_styles{k, 4}, 'Color', method_styles{k, 3}, ...
        'LineWidth', 1.25, 'MarkerSize', 4.2);
end
set(ax2, 'XScale', 'log', 'XLim', [1e-3 / 1.35, 1e1 * 1.35], ...
    'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1], ...
    'XTickLabel', {'$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$'}, ...
    'YTick', [0.4, 0.6, 0.8, 1.0, 1.2], ...
    'YTickLabel', {'$0.4$', '$0.6$', '$0.8$', '$1.0$', '$1.2$'}, ...
    'TickLabelInterpreter', 'latex', ...
    'XTickLabelRotation', 0, ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylim(ax2, [0.24, 1.24]);
xl2 = xlabel(ax2, 'noise level $\sigma$', 'Interpreter', 'latex', 'FontSize', 9.5);
xl2.Units = 'normalized';
xl2.Position = [0.5, -0.125, 0];
xl2.HorizontalAlignment = 'center';
yl2 = ylabel(ax2, 'median $\|U_I^{-1}\|_2 / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 10);
yl2.Units = 'normalized';
yl2.Position = [-0.18, 0.5, 0];
yl2.HorizontalAlignment = 'center';
t2 = title(ax2, sprintf('Perturbed Prop.~4.2 (median, $n = %d$)', S.n_trans), ...
    'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 10);
t2.Units = 'normalized';
t2.Position = [0.5, 1.045, 0];
t2.HorizontalAlignment = 'center';
t2.VerticalAlignment = 'baseline';

%% --- Unified two-row bottom legend ---
key = axes(fig, 'Units', 'inches', 'Position', [0.92, 0.04, 4.82, 0.44], ...
    'Visible', 'off', 'XLim', [0, 1], 'YLim', [0, 1]);
hold(key, 'on');
xs_leg = [0.00, 0.25, 0.53, 0.78, 0.00, 0.25, 0.53];
ys_leg = [0.76, 0.76, 0.76, 0.76, 0.24, 0.24, 0.24];
for k = 1:size(method_styles, 1)
    x0 = xs_leg(k); y0 = ys_leg(k);
    plot(key, x0 + [0, 0.050], [y0, y0], 'LineStyle', method_styles{k, 5}, ...
        'Color', method_styles{k, 3}, 'LineWidth', 1.3);
    plot(key, x0 + 0.025, y0, 'Marker', method_styles{k, 4}, ...
        'Color', method_styles{k, 3}, 'MarkerSize', 4.5, 'LineWidth', 1.2);
    text(key, x0 + 0.062, y0, method_styles{k, 2}, 'Interpreter', 'latex', ...
        'FontSize', 9.0, 'VerticalAlignment', 'middle');
end

out_file = fullfile(fig_dir, 'opt2_experiments.pdf');
save_pdf(fig, out_file);
fprintf('figures_opt2: wrote %s\n', out_file);
end

function fig = new_figure(width, height)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, width, height], ...
    'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesTickLabelInterpreter', 'none', 'DefaultTextInterpreter', 'none', ...
    'DefaultAxesFontSize', 9.0);
end

function ax = panel(fig, position, width, height)
ax = axes(fig, 'Units', 'normalized', 'Position', position ./ [width, height, width, height], ...
    'LineWidth', 0.7);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); ax.GridAlpha = 0.12;
end

function save_pdf(fig, file)
drawnow;
exportgraphics(fig, file, 'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
end
