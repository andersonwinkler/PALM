function w = lambertw(b, z)
%% Copyright (C) 1998 by Nicol N. Schraudolph
%%
%% This program is free software; you can redistribute and/or
%% modify it under the terms of the GNU General Public
%% License as published by the Free Software Foundation;
%% either version 3, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be
%% useful, but WITHOUT ANY WARRANTY; without even the implied
%% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
%% PURPOSE.  See the GNU General Public License for more
%% details.
%%
%% You should have received a copy of the GNU General Public
%% License along with this software; see the file COPYING.  If not,
%% see <https://www.gnu.org/licenses/>.

% ## -*- texinfo -*-
% ## @deftypefn {Function File} {@var{x} = } lambertw (@var{z})
% ## @deftypefnx {Function File} {@var{x} = } lambertw (@var{z}, @var{n})
% ## Compute the Lambert W function of @var{z}.
% ##
% ## This function satisfies W(z).*exp(W(z)) = z, and can thus be used to express
% ## solutions of transcendental equations involving exponentials or logarithms.
% ##
% ## @var{n} must be integer, and specifies the branch of W to be computed;
% ## W(z) is a shorthand for W(0,z), the principal branch.  Branches
% ## 0 and -1 are the only ones that can take on non-complex values.
% ##
% ## If either @var{n} or @var{z} are non-scalar, the function is mapped to each
% ## element; both may be non-scalar provided their dimensions agree.
% ##
% ## This implementation should return values within 2.5*eps of its
% ## counterpart in Maple V, release 3 or later.  Please report any
% ## discrepancies to the author, Nici Schraudolph <schraudo@@inf.ethz.ch>.
% ##
% ## For further details, see:
% ##
% ## Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), `On the Lambert
% ## W Function', Advances in Computational Mathematics 5(4):329-359.
% ## @end deftypefn

%% Author:   Nicol N. Schraudolph <schraudo@inf.ethz.ch>
%% Version:  1.0
%% Created:  07 Aug 1998
%% Keywords: Lambert W Omega special transcendental function
%% Notes:    Special-case handling aligned with MATLAB lambertw;
%%           minor edits for MATLAB/Octave portability.


  if nargin == 1
    z = b;
    b = 0;
  elseif nargin ~= 2
    error('lambertw: usage: lambertw(z) or lambertw(k, z)');
  else
    if any(round(real(b)) ~= b)
      error('lambertw: branch number must be integer');
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

  %% series expansion about -1/e
  %
  % p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  % w = (11/72)*p;
  % w = (w - 1/3).*p;
  % w = (w + 1).*p - 1
  %
  % first-order version suffices:
  %
  ee = exp(1);
  w = (1 - 2 * abs(b)) .* sqrt(2 * ee * z + 2) - 1;

  %% asymptotic expansion at 0 and Inf
  %
  % (~(z|b) and ~v) avoid log(0); written without ~ on complex values
  % so this runs under MATLAB as well as Octave.
  v = log(z + (z == 0 & b == 0)) + 2 * pi * 1i * b;
  v = v - log(v + (v == 0));

  %% choose strategy for initial guess
  %
  c = abs(z + 1 / ee);
  c = (c > 1.45 - 1.1 * abs(b));
  c = c | (b .* imag(z) > 0) | ((imag(z) == 0) & (b == 1));
  w = (1 - c) .* w + c .* v;

  %% Halley iteration
  %
  converged = false;
  for n = 1:10
    p = exp(w);
    t = w .* p - z;
    f = ones(size(w));
    f(w == -1) = 0;
    t = f .* t ./ (p .* (w + f) - 0.5 * (w + 2.0) .* t ./ (w + f));
    w = w - t;
    if all((abs(real(t)) < (2.48 * eps) * (1.0 + abs(real(w)))) ...
        &  (abs(imag(t)) < (2.48 * eps) * (1.0 + abs(imag(w)))))
      converged = true;
      break;
    end
  end

  %% special cases aligned with MATLAB lambertw
  w(isnan(w)) = NaN; % do not return NaN+NaN*1i
  w(isinf(z) & z > 0 & b == 0) = Inf;
  w(z == 0) = -Inf;
  w(z == 0 & b == 0) = 0;
  idx = (b == 0 & imag(z) == 0 & z > -exp(-1)) | ...
        (b == -1 & imag(z) == 0 & z > -exp(-1) & z < 0);
  w(idx) = real(w(idx));

  if ~converged
    overwritten = (z == 0) | isinf(z) | isnan(z);
    if any(~overwritten(:))
      warning('lambertw: iteration limit reached, result may be inaccurate');
    end
  end
end
