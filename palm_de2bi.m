% ## Copyright (C) 2001 Laurent Mazet
% ##
% ## This program is free software; you can redistribute it and/or modify it under
% ## the terms of the GNU General Public License as published by the Free Software
% ## Foundation; either version 3 of the License, or (at your option) any later
% ## version.
% ##
% ## This program is distributed in the hope that it will be useful, but WITHOUT
% ## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% ## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% ## details.
% ##
% ## You should have received a copy of the GNU General Public License along with
% ## this program; if not, see <http://www.gnu.org/licenses/>.
% 
% ## -*- texinfo -*-
% ## @deftypefn {Function File} {@var{b} = } de2bi (@var{d})
% ## @deftypefnx {Function File} {@var{b} = } de2bi (@var{d},@var{n})
% ## @deftypefnx {Function File} {@var{b} = } de2bi (@var{d},@var{n},@var{p})
% ## @deftypefnx {Function File} {@var{b} = } de2bi (@var{d},@var{n},@var{p},@var{f})
% ##
% ## Convert a non-negative integer to bit vector.
% ##
% ## The variable @var{d} must be a vector of non-negative integers. @dfn{de2bi}
% ## then returns a matrix where each row represents the binary representation
% ## of elements of @var{d}. If @var{n} is defined then the returned matrix
% ## will have @var{n} columns. This number of columns can be either larger
% ## than the minimum needed and zeros will be added to the msb of the
% ## binary representation or smaller than the minimum in which case the
% ## least-significant part of the element is returned.
% ##
% ## If @var{p} is defined then it is used as the base for the decomposition
% ## of the returned values. That is the elements of the returned value are
% ## between '0' and 'p-1'.
% ##
% ## The variable @var{f} defines whether the first or last element of @var{b}
% ## is considered to be the most-significant. Valid values of @var{f} are
% ## 'right-msb' or 'left-msb'. By default @var{f} is 'right-msb'.
% ##
% ## @seealso{bi2de}
% ## @end deftypefn

function b = palm_de2bi (d, n, p, f)

  if (nargin == 1)
    p = 2;
    n = floor ( log (max (max (d), 1)) ./ log (p) ) + 1;
    f = 'right-msb';
  elseif (nargin == 2)
    p = 2;
    f = 'right-msb';
  elseif (nargin == 3)
    if (ischar (p))
      f = p;
      p = 2;
    else
      f = 'right-msb';
    end
  elseif (nargin == 4)
    if (ischar (p))
      tmp = f;
      f = p;
      p = tmp;
    end
  else
    print_usage ();
  end

  d = d(:);
  if ( any (d < 0) || any (d ~= floor (d)) )
    error ('de2bi: d must only contain non-negative integers');
  end

  if (isempty (n))
    n = floor ( log (max (max (d), 1)) ./ log (p) ) + 1;
  end

  power = ones (length (d), 1) * (p .^ (0 : n-1) );
  d = d * ones (1, n);
  b = floor (rem (d, p*power) ./ power);

  if (strcmp (f, 'left-msb'))
    b = b(:,columns(b):-1:1);
  elseif (~strcmp (f, 'right-msb'))
    error ('de2bi: unrecognized flag');
  end

end
