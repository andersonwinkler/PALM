function Q = palm_semiortho(Z,Sel)
% Compute a semi-orthogonal matrix according to
% the Huh-Jhun or Theil methods.
% 
% Q = palm_semiortho(Z,Sel)
% 
% Z   : Matrix of nuisance variables
% Sel : Selection matrix
% Q   : Semi-orthogonal matrix
% 
% _______________________________________
% Anderson M. Winkler & Thomas E. Nichols
% NIH/NIMH & Univ. of Oxford
% Mar/2020 (first version)
% Apr/2026 (this version)
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
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

% Note that, due to a
% simplification of HJ, input here is Z, not Rz.
if isempty(Sel)
    % If Sel is empty, do Huh-Jhun
    % HJ here is simplified as in Winkler et al, 2020 (see the Appendix text of the paper)
    [Q,D,~] = svd(null(Z'),'econ');
    Q = Q*D;
else
    % Theil
    [N,R] = size(Z);
    if isvector(Sel)
        % If Sel is a vector of logical or integer indices
        if islogical(Sel)
            Sel = find(Sel);
        end
        if Sel(1) > 0
            % If Sel is a column of indices
            unSel = setdiff(1:N,Sel);
            if rank(Z(unSel,:)) < R
                error('Selected rows of nuisance not full rank')
            end
        else
            % If Sel is -1 or anything else but empty [].
            % First confirm that this is even possible
            rZ = rank(Z);
            if rZ < R
                error('Impossible to use the Theil method with this set of nuisance variables')
            end

            % Find unique rows; since unique sorts outputs, shuffle to avoid trends
            [~,iU,~] = unique(Z,'rows');
            pidx = randperm(numel(iU));
            iU = iU(pidx);

            % Go by trial and error
            unSel0 = [];
            rnk0   = 0;
            for u = iU'
                unSel = [unSel0 u];
                Zout  = Z(unSel,:);
                rnk   = rank(Zout);
                if rnk > rnk0
                    unSel0 = unSel;
                    rnk0   = rnk;
                end
                if rnk == R
                    break
                end
            end
            Sel = setdiff((1:N),unSel);
        end
        S = eye(N);
        S = S(:,Sel);
    else
        % Sel is a matrix proper
        S = Sel;
    end
    Rz = eye(N) - Z*pinv(Z);
    Q = Rz*S*sqrtm(inv(S'*Rz*S));
end