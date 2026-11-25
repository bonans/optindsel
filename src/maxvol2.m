function ind = maxvol2(a, ind, varargin)
%maxvol2 - Greedy 1-swap maximum-volume row selection
%   ind = maxvol2(a) computes row indices ind of a tall n-by-r matrix a
%   starting from partial-pivoted LU and swapping single rows until
%   max_{i,j} |(a * a(ind,:)^(-1))_{i,j}| <= 1 + tol (default tol = 5e-2).
%   ind = maxvol2(a, ind, 'eps', tol, 'niters', niters, 'qr', do_qr)
%   customizes the initial row set, swap tolerance, iteration cap, or
%   initial thin QR orthogonalization.
%
%   Adapted from TT-Toolbox:
%   https://github.com/oseledets/TT-Toolbox/blob/master/core/maxvol2.m
%   Modifications: respect user-supplied initial ind when nonempty, add the
%   missing iter increment inside the while loop, guard tol = 0 against
%   unit-roundoff self-swaps, and remove the unused external lu_full hook.
%
%   Reference:
%   S. A. Goreinov, I. V. Oseledets, D. V. Savostyanov, E. E. Tyrtyshnikov,
%   and N. L. Zamarashkin, "How to find a good submatrix," in Matrix
%   Methods: Theory, Algorithms and Applications, World Scientific, 2010,
%   pp. 247-256.
%
%   MIT License
%   Copyright (C) 2009-2012 Ivan Oseledets, Sergey Dolgov, Vladimir Kazeev,
%   Thomas Mach, Olga Lebedeva, Dmitry Savostyanov, Pavel Zhlobich, Le Song
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions: The above copyright notice and this
%   permission notice shall be included in all copies or substantial
%   portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT
%   WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
%   THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
%   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
%   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
%   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

n = size(a, 1);
r = size(a, 2);
if n <= r
    ind = 1:n;
    return;
end

do_qr = false;
niters = 100;
tol = 5e-2;

for i = 1:2:(length(varargin) - 1)
    switch lower(varargin{i})
        case 'qr'
            do_qr = varargin{i + 1};
        case 'eps'
            tol = varargin{i + 1};
        case 'niters'
            niters = varargin{i + 1};
        case 'ind'
            ind = varargin{i + 1};
        otherwise
            error('maxvol2:unrecognizedOption', 'Unrecognized option: %s', varargin{i});
    end
end

if do_qr
    [a, ~] = qr(a, 0);
end

if nargin < 2 || isempty(ind)
    [~, ~, p] = lu(a, 'vector');
    ind = p(1:r);
end

sbm = a(ind, :);
b = a / sbm;

thresh = 1 + max(tol, 10 * eps(class(a)));
iter = 0;
while iter < niters
    [mx0, big_ind] = max(abs(b(:)));
    [i0, j0] = ind2sub([n, r], big_ind);
    if mx0 <= thresh
        ind = reshape(sort(ind), 1, r);
        return;
    end
    k = ind(j0);
    b = b + b(:, j0) * ((b(k, :) - b(i0, :)) / b(i0, j0));
    ind(j0) = i0;
    iter = iter + 1;
end

ind = reshape(sort(ind), 1, r);
end
