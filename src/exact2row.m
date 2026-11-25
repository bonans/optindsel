function [I, wnorm, invnorm] = exact2row(U, w)
%exact2row - Exhaustive search for optimal 2-by-2 row selection
%   [I, wnorm, invnorm] = exact2row(U, w) enumerates all n*(n-1)/2 row
%   pairs to minimize ||U(I,:)^(-1) * diag(sqrt(w(I)))||_2 (default
%   w = ones(n, 1) / n).

n = size(U, 1);
if nargin < 2 || isempty(w)
    A = U * sqrt(n);
else
    A = U ./ sqrt(w(:));
end
wnorm = Inf;
I = [1, 2];
a1 = A(:, 1);
a2 = A(:, 2);
for i = 1:(n - 1)
    b11 = a1(i);
    b12 = a2(i);
    b21 = a1((i+1):n);
    b22 = a2((i+1):n);
    den = 2 * abs(b11 * b22 - b12 * b21);
    num = hypot(b11 + b22, b12 - b21) + hypot(b11 - b22, b12 + b21);
    vals = num ./ den;
    vals(den == 0) = Inf;
    [vmin, rel] = min(vals);
    if vmin < wnorm
        wnorm = vmin;
        I = [i, i + rel];
    end
end
if nargout >= 3
    invnorm = inv2norm(U(I, :));
end
end
