function [I, invnorm] = exactmaxvol2row(U)
%exactmaxvol2row - Exact global maximum-volume 2-by-2 row selection
%   [I, invnorm] = exactmaxvol2row(U) finds I = [i, j] maximizing
%   |det U([i, j], :)| over all 1 <= i < j <= n.

n = size(U, 1);
best_det = -1;
I = [1, 2];
u1 = U(:, 1);
u2 = U(:, 2);
for i = 1:(n - 1)
    dets = abs(u1(i) * u2((i+1):n) - u2(i) * u1((i+1):n));
    [dmax, rel] = max(dets);
    if dmax > best_det
        best_det = dmax;
        I = [i, i + rel];
    end
end
if nargout >= 2
    invnorm = inv2norm(U(I, :));
end
end
