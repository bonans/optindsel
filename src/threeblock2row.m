function [U, z, V] = threeblock2row(p1, p2, p3, z, V)
%threeblock2row - Explicit three-block n-by-2 orthonormal-column basis
%   [U, z, V] = threeblock2row(p1, p2, p3, z, V) constructs
%   U_{p_1,p_2,p_3}(z) in Section 4 with positive integer block sizes
%   p1, p2, p3 (n = p1 + p2 + p3), unit null vector z in R^3, and 3-by-2
%   orthonormal basis V of z^perp.
%
%   When z is omitted or empty, uses the GTZ equality vector
%   z_* = [sqrt((n-p1)/(2n)); sqrt((n-p2)/(2n)); sqrt((n-p3)/(2n))]
%   (Theorem 4.1).
%   When V is omitted or empty, uses the explicit basis V(z) with V(1,2)=0
%   whenever |z(1)| < 1.

if any([p1, p2, p3] < 1) || any(floor([p1, p2, p3]) ~= [p1, p2, p3])
    error('threeblock2row:invalidBlocks', 'Block sizes p1, p2, p3 must be positive integers.');
end
n = p1 + p2 + p3;

if nargin < 4 || isempty(z)
    z = [sqrt((n - p1) / (2 * n)); sqrt((n - p2) / (2 * n)); sqrt((n - p3) / (2 * n))];
else
    z = z(:) / norm(z(:));
end

if nargin < 5 || isempty(V)
    s1 = sqrt(max(0, 1 - z(1)^2));
    if s1 > 1e-14
        V = [ s1,                0; ...
             -z(1)*z(2) / s1,  z(3) / s1; ...
             -z(1)*z(3) / s1, -z(2) / s1];
    else
        V = [0, 0; 1, 0; 0, 1];
    end
end

U = [repmat(V(1, :) / sqrt(p1), p1, 1); ...
     repmat(V(2, :) / sqrt(p2), p2, 1); ...
     repmat(V(3, :) / sqrt(p3), p3, 1)];
end
