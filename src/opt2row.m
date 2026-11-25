function [I, wnorm, invnorm] = opt2row(U, w)
%opt2row - Streaming 2-by-2 row selection in O(n) operations (Algorithm 2)
%   [I, wnorm, invnorm] = opt2row(U) selects a 2-element row index pair
%   I = [i, j] (1 <= i < j <= n) from an n-by-2 matrix U with orthonormal
%   columns such that ||U(I,:)^(-1)||_2 <= sqrt(n) in exact arithmetic.
%
%   [I, wnorm, invnorm] = opt2row(U, w) uses positive weights w(1:n)
%   satisfying sum(w) <= 1 and guarantees
%       wnorm = ||A(I,:)^(-1)||_2 = ||U(I,:)^(-1) * diag(sqrt(w(I)))||_2 <= 1
%   in exact arithmetic, where A = diag(w.^(-1/2)) * U. Also returns
%   invnorm = ||U(I,:)^(-1)||_2.

n = size(U, 1);
if size(U, 2) ~= 2 || n < 2
    error('opt2row:invalidInput', 'U must be an n-by-2 matrix with n >= 2.');
end
unweighted = (nargin < 2 || isempty(w));
if unweighted
    w_unif = 1 / n;
    sw_unif = sqrt(w_unif);
elseif numel(w) ~= n
    error('opt2row:invalidWeights', 'w must be a finite positive vector of length n.');
end

u    = eps(class(U)) / 2; % unit roundoff
rmin = realmin(class(U));

% Preallocated storage for the active subset P of at most 4 rows (m = |P|).
% Storing U_sub, sqrt(w_r), Ahat = u_r / ||u_r||_2, and ainv = ||a_r||_2^(-1)
% allows U and w to be read in a single sequential pass while normalizing
% intermediate quantities to reduce the risk of overflow and underflow.
P     = zeros(1, 4);
U_sub = zeros(4, 2);
Ahat  = zeros(4, 2);
sw    = zeros(4, 1);
w_sub = zeros(4, 1);
ainv  = zeros(4, 1);
m     = 0;

% Upper-triangular index pairs (k, l) with 1 <= k < l <= 4
idx_k     = [1, 1, 1, 2, 2, 3];
idx_l     = [2, 3, 4, 3, 4, 4];
upper_idx = idx_k + 4 * (idx_l - 1);

for r = 1:n
    if unweighted
        wr = w_unif;
        swr = sw_unif;
    else
        wr = w(r);
        if ~(wr > 0 && isfinite(wr))
            error('opt2row:invalidWeights', 'w must be a finite positive vector of length n.');
        end
        swr = sqrt(wr);
    end

    ur = U(r, :);
    nrm_u = hypot(ur(1), ur(2));
    if nrm_u <= rmin
        continue;
    end

    m = m + 1;
    P(m)        = r;
    U_sub(m, :) = ur;
    Ahat(m, :)  = ur / nrm_u;
    sw(m)       = swr;
    w_sub(m)    = wr;
    ainv(m)     = swr / nrm_u;

    if m == 4
        % Pairwise determinants d_{kl} for 1 <= k < l <= 4
        D = Ahat(:, 1) * Ahat(:, 2).' - Ahat(:, 2) * Ahat(:, 1).';
        [d_min, pos] = min(abs(D(upper_idx)));

        if d_min <= 4 * u
            kp = idx_k(pos);
            lp = idx_l(pos);
            s_max = max(ainv(kp), ainv(lp));
            gamma = zeros(4, 1);
            gamma(kp) =  (ainv(kp) / s_max)^2;
            gamma(lp) = -(ainv(lp) / s_max)^2;
        else
            gamma = [ D(2,3) * D(2,4) * D(3,4); ...
                     -D(1,3) * D(1,4) * D(3,4); ...
                      D(1,2) * D(1,4) * D(2,4); ...
                     -D(1,2) * D(1,3) * D(2,3)];
            gnorm = norm(gamma, inf);
            if gnorm > 0
                gamma = gamma / gnorm;
            end
            gamma = gamma .* (ainv / max(ainv)).^2;
            gnorm = norm(gamma, inf);
            if gnorm > 0
                gamma = gamma / gnorm;
            end
        end

        if sum(gamma) > 0
            gamma = -gamma;
        end

        idx = find(gamma < 0);
        if ~isempty(idx)
            [theta_star, j] = min(w_sub(idx) ./ (-gamma(idx)));
            k_star = idx(j);
            w_sub = w_sub + theta_star * gamma;
            w_sub(k_star) = 0;
        else
            [~, k_star] = min(w_sub);
            w_sub(k_star) = 0;
        end

        keep = w_sub > 0;
        m = sum(keep);
        P(1:m)       = P(keep);
        U_sub(1:m,:) = U_sub(keep, :);
        Ahat(1:m, :) = Ahat(keep, :);
        sw(1:m)      = sw(keep);
        w_sub(1:m)   = w_sub(keep);
        ainv(1:m)    = ainv(keep);
    end
end

if m < 2
    error('opt2row:rankDeficient', 'Input matrix U does not have rank 2.');
end

% Return I = {i, j} in P minimizing ||A_I^(-1)||_2 evaluated via inv2norm
A_sub = U_sub(1:m, :) ./ sw(1:m);

wnorm = Inf;
best_pair = [1, 2];
for k = 1:(m - 1)
    for l = (k + 1):m
        val = inv2norm(A_sub([k, l], :));
        if val < wnorm
            wnorm = val;
            best_pair = [k, l];
        end
    end
end

[I, ord] = sort(P(best_pair));
if nargout >= 3
    pair = best_pair(ord);
    invnorm = inv2norm(U_sub(pair, :));
end
end
