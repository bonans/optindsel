function run_opt2()
%run_opt2 - Run numerical experiments for "Selecting 2x2 submatrices with optimally bounded inverses"
%   Writes data/opt2_experiments.mat containing:
%     1. Perturbed three-block matrices across 8 <= n <= 1024 with log-uniform
%        block sizes and uniform unit vectors z (2,000 realizations per n at
%        noise level sigma = 0.02).
%     2. Perturbed Proposition 4.2 counterexample at n = 512 across
%        noise levels 10^(-3) <= sigma <= 10^1 (2,000 realizations per sigma).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));
data_dir = fullfile(root, 'data');
if ~isfolder(data_dir)
    mkdir(data_dir);
end

methods = ["opt2row", "exact_maxvol", "maxvol", "qrcp", "bsqr", "lupp", "exact"];
m_count = numel(methods);

%% --- 1. Perturbed three-block ensemble (8 <= n <= 1024) ---
rng(314159, 'twister');
ns_cond = [8, 16, 32, 64, 128, 256, 512, 1024];
n_trials_cond = 2000;
noise_level = 0.02;

cond_max       = zeros(numel(ns_cond), m_count);
cond_med       = zeros(numel(ns_cond), m_count);
cond_ratio_max = zeros(numel(ns_cond), m_count - 1);
cond_ratio_med = zeros(numel(ns_cond), m_count - 1);

for k = 1:numel(ns_cond)
    n = ns_cond(k);
    sqn = sqrt(n);
    vals = zeros(n_trials_cond, m_count);
    for t = 1:n_trials_cond
        w3 = n.^(rand(1, 3));
        w3 = w3 / sum(w3);
        p1 = max(1, min(n - 2, round(w3(1) * n)));
        p2 = max(1, min(n - p1 - 1, round(w3(2) * n)));
        p3 = n - p1 - p2;

        z = randn(3, 1);
        z = z / norm(z);
        [Q3, ~] = qr([z, randn(3, 2)]);
        V = Q3(:, 2:3);
        U0 = threeblock2row(p1, p2, p3, z, V);

        sgn = 2 * (rand(n, 1) > 0.5) - 1;
        X = sgn .* U0(randperm(n), :) + (noise_level / sqn) * randn(n, 2);
        [U, ~] = qr(X, 0);
        vals(t, :) = eval_all_methods(U);
    end
    normed = vals / sqn;
    cond_max(k, :) = max(normed, [], 1);
    cond_med(k, :) = median(normed, 1);
    rel_opt = vals(:, 1:6) ./ vals(:, 7);
    cond_ratio_max(k, :) = max(rel_opt, [], 1);
    cond_ratio_med(k, :) = median(rel_opt, 1);
end

%% --- 2. Perturbed Proposition 4.2 transition across sigma at n = 512 ---
rng(271828, 'twister');
n_trans = 512;
sqn_trans = sqrt(n_trans);
n1 = n_trans - 1;
m_blk = (n1 - 1) / 2;
z_adv = [sqrt(1 - 4 / n1); sqrt(2 / n1); sqrt(2 / n1)];
U_star = [threeblock2row(1, m_blk, m_blk, z_adv); 0, 0];

sigmas = 10.^(-3 : 0.5 : 1);
n_trials_trans = 2000;
trans_med = zeros(numel(sigmas), m_count);
trans_max = zeros(numel(sigmas), m_count);
trans_q25 = zeros(numel(sigmas), m_count);
trans_q75 = zeros(numel(sigmas), m_count);

for s = 1:numel(sigmas)
    sig = sigmas(s);
    vals_s = zeros(n_trials_trans, m_count);
    for t = 1:n_trials_trans
        sgn = 2 * (rand(n_trans, 1) > 0.5) - 1;
        X = sgn .* U_star(randperm(n_trans), :) + (sig / sqn_trans) * randn(n_trans, 2);
        [U, ~] = qr(X, 0);
        vals_s(t, :) = eval_all_methods(U);
    end
    normed_s = vals_s / sqn_trans;
    trans_med(s, :) = median(normed_s, 1);
    trans_max(s, :) = max(normed_s, [], 1);
    trans_q25(s, :) = prctile(normed_s, 25, 1);
    trans_q75(s, :) = prctile(normed_s, 75, 1);
end

save(fullfile(data_dir, 'opt2_experiments.mat'), ...
    'methods', 'ns_cond', 'n_trials_cond', 'noise_level', ...
    'cond_max', 'cond_med', 'cond_ratio_max', 'cond_ratio_med', ...
    'n_trans', 'sigmas', 'n_trials_trans', ...
    'trans_med', 'trans_max', 'trans_q25', 'trans_q75');

fprintf('run_opt2: saved data/opt2_experiments.mat.\n');
fprintf('  1. Perturbed 3-block ensemble (max across 8 <= n <= 1024, 2000 trials):\n');
for m = 1:m_count
    fprintf('    %-13s: global_max=%.4f, med_at_512=%.4f\n', ...
        methods(m), max(cond_max(:, m)), cond_med(end - 1, m));
end
fprintf('  2. Perturbed Prop 4.2 transition (median at n = %d, %d trials):\n', n_trans, n_trials_trans);
for s = 1:numel(sigmas)
    fprintf('    sigma=%7.4f: opt=%.4f, emv=%.4f, mv=%.4f, qrcp=%.4f, bsqr=%.4f, lupp=%.4f, ex=%.4f\n', ...
        sigmas(s), trans_med(s, 1), trans_med(s, 2), trans_med(s, 3), ...
        trans_med(s, 4), trans_med(s, 5), trans_med(s, 6), trans_med(s, 7));
end
end

function vals = eval_all_methods(U)
[~, ~, v_opt] = opt2row(U);
[~, v_emv]    = exactmaxvol2row(U);
I_mv          = maxvol2(U);
v_mv          = inv2norm(U(I_mv, :));
[~, ~, p_qr]  = qr(U', 'vector');
v_qr          = inv2norm(U(p_qr(1:2), :));
[~, ~, p_bs]  = bsqr(U', 2);
v_bs          = inv2norm(U(p_bs(1:2), :));
[~, ~, p_lu]  = lu(U, 'vector');
v_lu          = inv2norm(U(p_lu(1:2), :));
[~, ~, v_ex]  = exact2row(U);
vals = [v_opt, v_emv, v_mv, v_qr, v_bs, v_lu, v_ex];
end
