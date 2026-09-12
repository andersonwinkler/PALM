function w = palm_lambertw(b,z,varargin)
% Lambert W function.
%
% Usage:
% w = palm_lambertw(z)
% w = palm_lambertw(k,z)
%
% Inputs:
% z     : Argument of W. W(z).*exp(W(z)) = z.
% k     : Integer branch. Default is 0 (principal branch).
%         Branches 0 and -1 are the only ones that can be real.
%         K and Z must be the same size, or be scalars.
%
% Outputs:
% w     : Value of the Lambert W function.
%
% This is a MATLAB/Octave port of lambertw by Nicol N. Schraudolph.
%
% Reference:
% * Corless RM, Gonnet GH, Hare DEG, Jeffrey DJ, Knuth DE.
%   On the Lambert W function. Adv Comput Math. 1996;5(4):329-359.
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Sep/2026
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 1998 Nicol N. Schraudolph
% Copyright (C) 2026 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin == 1
    z = b;
    b = 0;
elseif nargin ~= 2
    if nargin < 1
        error('Too few input arguments.');
    else
        error('Too many input arguments');
    end
else
    if any(round(real(b)) ~= b)
        error('Branch number must be integer');
    end
end

if isempty(b)
    w = b;
    return;
end
if isempty(z)
    w = z;
    return;
end
if all(isnan(z(:)))
    w = NaN(size(z));
    return;
end

% Series expansion about -1/e
ee = exp(1);
w = (1 - 2*abs(b)) .* sqrt(2*ee*z + 2) - 1;

% Asymptotic expansion at 0 and Inf
v = log(z + (z==0 & b==0)) + 2*pi*1i*b;
v = v - log(v + (v==0));

% Initial guess
c = abs(z + 1/ee);
c = (c > 1.45 - 1.1*abs(b));
c = c | (b.*imag(z) > 0) | ((imag(z)==0) & (b==1));
w = (1-c).*w + c.*v;

% Halley iteration
converged = false;
for n = 1:10
    p = exp(w);
    t = w.*p - z;
    f = ones(size(w));
    f(w==-1) = 0;
    t = f .* t ./ (p.*(w+f) - 0.5*(w+2.0).*t./(w+f));
    w = w - t;
    if all((abs(real(t)) < (2.48*eps)*(1.0+abs(real(w)))) ...
        &  (abs(imag(t)) < (2.48*eps)*(1.0+abs(imag(w)))))
        converged = true;
        break;
    end
end

% Special cases aligned with MATLAB lambertw
w(isnan(w)) = NaN;
w(isinf(z) & z>0 & b==0) = Inf;
w(z==0) = -Inf;
w(z==0 & b==0) = 0;
idx = (b==0 & imag(z)==0 & z>-exp(-1)) | ...
      (b==-1 & imag(z)==0 & z>-exp(-1) & z<0);
w(idx) = real(w(idx));

% If not converged...
if ~converged
    overwritten = (z==0) | isinf(z) | isnan(z);
    if any(~overwritten(:))
        warning('Iteration limit reached, result may be inaccurate');
    end
end
