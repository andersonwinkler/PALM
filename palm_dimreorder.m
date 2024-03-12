function Y = palm_dimreorder(X,sqz)
% Permute the dimensions of data read with palm_miscread such that
% the permutable dimension is the last. If the data had been permuted
% already, then reorder back to the original.
%
% Usage:
% Y = palm_dimreorder(X,sqz)
%
% X   : Data structure from palm_miscread.
% szq : Squeeze the permutable dimension? (true/false)
% Y   : Like X, but with the dimensions permuted.
%
% For reordering an array only, see palm_convNto2 and palm_conv2toN.
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Mar/2024
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2024 Anderson M. Winkler
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

switch X.readwith

    case 'cifti-matlab'

        % Operate on a copy
        Y = X;

        if isfield(X.extra,'oldorder')

            % If a field specifying an old order is present, either squeeze
            % the last dimension, or reorder back to original
            if sqz

                % Cannot drop to less than two dimensions
                if numel(Y.extra.diminfo) > 2
                    Y.extra.diminfo(end) = [];
                end
                idxdimdrop = Y.extra.oldorder == numel(Y.extra.oldorder);
                Y.extra.oldorder(idxdimdrop) = [];

                % Adjust the file extension
                switch lower(Y.extra.cifti_file_extension)
                    case {'pconnscalar','pconnseries'}
                        Y.extra.cifti_file_extension = 'pconn';
                    case {'ptseries'}
                        Y.extra.cifti_file_extension = 'pscalar';
                    case {'dtseries'}
                        Y.extra.cifti_file_extension = 'dscalar';
                end
            else
                % Reorder the dims back, keeping the file type
                Y.data = permute(X.data,X.extra.oldorder);
                Y.extra.diminfo = Y.extra.diminfo(X.extra.oldorder);
                Y.extra = rmfield(Y.extra,'oldorder');
            end
        else

            % Otherwise, we want to reorder such that the permutable
            % dimension is the last.
            % Locate the dimension type for each dimension
            nD = numel(X.extra.diminfo);
            Dtypes = cell(1,nD);
            for d = 1:nD
                Dtypes{d} = X.extra.diminfo{d}.type;
            end

            % Identify the permutable dimension
            pdim = find(strcmpi('scalars',Dtypes) | strcmpi('series',Dtypes));
            if numel(pdim) ~= 1
                error('Cannot determine the permutation dimension in CIFTI file: %s', X.filename);
            end

            % Reorder the data such that the permutable dimension is the last
            neworder = [setdiff(1:nD,pdim) pdim];
            Y.data = permute(X.data,neworder);
            Y.extra.diminfo = Y.extra.diminfo(neworder);
            [~,Y.extra.oldorder] = sort(neworder);
        end

    otherwise
        error('This function currently does not work with files read with %s.', X.readwith);
end
