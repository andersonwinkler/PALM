function Z = palm_inormal(varargin)
% Applies a rank-based inverse normal transformation.
% 
% Usage: Z = inormal(X)
%            inormal(X,c)
%            inormal(X,method)
%
% X      : Original data. Can be a vector or an array.
% c      : Constant to be used in the transformation.
%          Default c=3/8 (Blom).
% method : Method to choose c. Accepted values are:
%              'Blom'   (c=3/8),
%              'Tukey'  (c=1/3),
%              'Bliss', (c=1/2)  and
%              'Waerden' or 'SOLAR' (c=0). 
% Z      : Transformed data.
% 
% References:
% * Van der Waerden BL. Order tests for the two-sample
%   problem and their power. Proc Koninklijke Nederlandse
%   Akademie van Wetenschappen. Ser A. 1952; 55:453�458
% * Blom G. Statistical estimates and transformed
%   beta-variables. Wiley, New York, 1958.
% * Tukey JW. The future of data analysis.
%   Ann Math Stat. 1962; 33:1�67.
% * Bliss CI. Statistics in biology. McGraw-Hill,
%   New York, 1967.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011
% http://brainder.org


% Accept inputs & defaults
if nargin == 1,
    c = 3/8;  % Default (Blom, 1958)
elseif nargin == 2 && ischar(varargin{2}),
    switch lower(varargin{2}),
        case 'blom'
            c = 3/8;
        case 'tukey'
            c = 1/2;
        case 'bliss'
            c = 1/2;
        case 'waerden'
            c = 0;
        case 'solar'
            c = 0; % SOLAR is the same as Van der Waerden
        otherwise
            error('Method %s unknown. Use either ''Blom'', ''Tukey'', ''Bliss'', ''Waerden'' or ''SOLAR''.',varargin{2});
    end
elseif nargin == 2 && isscalar(varargin{2}),
    c = varargin{2};  % For a user-specified value for c
else
    error('Invalid arguments.')
end
X = varargin{1};

% Get the rank for each value
[~,iX] = sort(X);
[~,ri] = sort(iX);

% Do the actual transformation
N = size(X,1);
p = ((ri-c)/(N-2*c+1));
Z = sqrt(2)*erfinv(2*p-1);

