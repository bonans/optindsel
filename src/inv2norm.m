function nrm = inv2norm(B)
%inv2norm - Spectral norm ||B^(-1)||_2 of a 2-by-2 real matrix B
%   nrm = inv2norm(B) evaluates ||B^(-1)||_2 = sigma_1(B) / |det B| via
%   Proposition 3.2 without subtractive cancellation in sigma_2(B),
%   returning Inf when det B == 0.

detB = abs(B(1,1)*B(2,2) - B(1,2)*B(2,1));
if detB == 0
    nrm = Inf;
    return;
end
s_plus  = hypot(B(1,1) + B(2,2), B(1,2) - B(2,1));
s_minus = hypot(B(1,1) - B(2,2), B(1,2) + B(2,1));
nrm     = (s_plus + s_minus) / (2 * detB);
end
