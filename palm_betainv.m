function x = palm_betainv(p,a,b,varargin)
% Inverse of the beta distribution.
%
% Usage:
% x = palm_betainv(p,a,b)
%
% Inputs:
% p     : Probabilities in [0,1].
% a, b  : Shape parameters.
%         p, a, and b must be the same size, or be scalars.
%
% Outputs:
% x     : Quantiles.
%
% A scalar input is treated as a constant array of the same size
% as the other inputs. Invalid values result in NaN.
%
% This is a MATLAB/Octave port of betainv from the GNU Octave
% statistics package, using only core functions (betainc, betaln)
% so that the statistics package/toolbox is not required.
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Sep/2026
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 1995-2016 Kurt Hornik
% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 2023 Andreas Bertsatos
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

% Check for valid number of input arguments
if nargin < 3
    error('Too few input arguments.');
elseif nargin > 3
    error('Too many input arguments');
end

% Check for common size
if ~isscalar(p) || ~isscalar(a) || ~isscalar(b)
    [retval,p,a,b] = common_size(p,a,b);
    if retval > 0
        error('Arguments must be of common size or scalars.');
    end
end

% Must be real
if ~isreal(p) || ~isreal(a) || ~isreal(b)
    error('Arguments must be real.');
end

% Check for class type
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x = zeros(size(p),'single');
else
    x = zeros(size(p));
end

k = (p<0) | (p>1) | ~(a>0) | ~(b>0) | isnan(p);
x(k) = NaN;

k = (p==1) & (a>0) & (b>0);
x(k) = 1;

k = find((p>0) & (p<1) & (a>0) & (b>0));
if ~isempty(k)
    if ~isscalar(a) || ~isscalar(b)
        a = a(k);
        b = b(k);
        y = a ./ (a+b);
    else
        y = a / (a+b) * ones(size(k));
    end
    p = p(k);

    if isa(y, 'single')
        myeps = eps('single');
    else
        myeps = eps;
    end

    l = find(y < myeps);
    if any(l)
        y(l) = sqrt(myeps) * ones(length(l),1);
    end
    l = find(y > 1 - myeps);
    if any(l)
        y(l) = 1 - sqrt(myeps) * ones(length(l),1);
    end

    % Newton iteration for the quantile
    y_new = y;
    loopcnt = 0;
    while true
        y_old = y_new;
        h     = (betacdf(y_old,a,b) - p) ./ betapdf(y_old,a,b);
        y_new = y_old - h;
        ind   = find(y_new <= myeps);
        if any(ind)
            y_new(ind) = y_old(ind) / 10;
        end
        ind = find(y_new >= 1-myeps);
        if any(ind)
            y_new(ind) = 1 - (1-y_old(ind)) / 10;
        end
        h = y_old - y_new;
        loopcnt = loopcnt + 1;
        if max(abs(h)) < sqrt(myeps) || loopcnt == 40
            break;
        end
    end
    if loopcnt == 40
        warning('Failed to converge for some values.');
    end
    x(k) = y_new;
end

% -------------------------------------------------------------------------
function [err, p, a, b] = common_size(p, a, b)
% Expand scalars to a common size, as in Octave's common_size.
args = {p, a, b};
is_array = cellfun(@numel, args) ~= 1;
aridx = find(is_array, 1);
if isempty(aridx)
    err = 0;
    return;
end
target = size(args{aridx});
for i = 1:3
    if is_array(i) && ~isequal(size(args{i}), target)
        err = 1;
        return;
    end
end
err = 0;
if ~is_array(1), p = repmat(p,target); end
if ~is_array(2), a = repmat(a,target); end
if ~is_array(3), b = repmat(b,target); end

% -------------------------------------------------------------------------
function cdf = betacdf(x,a,b)
% Regularized incomplete beta (core MATLAB/Octave).
cdf = betainc(x,a,b);

% -------------------------------------------------------------------------
function pdf = betapdf(x,a,b)
% Beta density on (0,1) via betaln (core MATLAB/Octave).
pdf = exp((a-1) .* log(x) + (b-1) .* log(1-x) - betaln(a,b));
