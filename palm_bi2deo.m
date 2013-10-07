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
% ## @deftypefn {Function File} {@var{d} = } bi2de (@var{b})
% ## @deftypefnx {Function File} {@var{d} = } bi2de (@var{b},@var{f})
% ## @deftypefnx {Function File} {@var{d} = } bi2de (@var{b},@var{p})
% ## @deftypefnx {Function File} {@var{d} = } bi2de (@var{b},@var{p},@var{f})
% ##
% ## Convert bit matrix to a vector of integers
% ##
% ## Each row of the matrix @var{b} is treated as a single integer represented
% ## in binary form. The elements of @var{b}, must therefore be '0' or '1' 
% ##
% ## If @var{p} is defined then it is treated as the base of the decomposition
% ## and the elements of @var{b} must then lie between '0' and 'p-1'.
% ##
% ## The variable @var{f} defines whether the first or last element of @var{b}
% ## is considered to be the most-significant. Valid values of @var{f} are
% ## 'right-msb' or 'left-msb'. By default @var{f} is 'right-msb'.
% ##
% ## @seealso{de2bi}
% ## @end deftypefn

function d = palm_bi2deo (b, p, f)

  switch (nargin)
    case 1,
      p = 2;
      f = 'right-msb';
    case 2,
      if (ischar (p))
        f = p;
        p = 2;
      else
        f = 'right-msb';
      end
    case 3,
      if (ischar (p))
        tmp = f;
        f = p;
        p = tmp;
      end
    otherwise
      print_usage ();
  end

  if ( any (b(:) < 0) || any (b(:) ~= floor (b(:))) || any (b(:) > p - 1) )
    error ('bi2de: d must only contain integers in the range [0, p-1]');
  end

  if (strcmp (f, 'left-msb'))
    b = b(:,size(b,2):-1:1);
  elseif (~strcmp (f, 'right-msb'))
    error ('bi2de: unrecognized flag');
  end

  if (isempty (b))
    d = [];
  else
    d = b * ( p .^ ( 0 : (columns(b)-1) )' );
  end

end
