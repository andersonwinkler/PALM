function [X,Z,eCm,eCx] = palm_partition(M,C,meth,Y)
% Partition a design matrix into regressors of interest and
% nuisance according to a given contrast.
% 
% Usage
% [X,Z] = partition(M,C,meth,Y)
% 
% Inputs:
% M    : Design matrix, to be partitioned.
% C    : Contrast that will define the partitioning.
% meth : Method for the partitioning. It can be:
%        - 'Guttman'
%        - 'Beckmann'
%        - 'Winkler'
%        - 'Ridgway'
% Y    : (Optional) For the 'Winkler' method only.
% 
% Outputs:
% X    : Matrix with regressors of interest.
% Z    : Matrix with regressors of no interest.
% eCm  : Effective contrast, equivalent to the original,
%        for the partitioned model [X Z], and considering
%        all regressors.
% eCx  : Same as above, but considering only X.
%
% References:
% * Guttman I. Linear Models: An Introduction. Wiley,
%   New York, 1982.
% * Smith SM, Jenkinson M, Beckmann C, Miller K,
%   Woolrich M. Meaningful design and contrast estimability
%   in FMRI. Neuroimage 2007;34(1):127-36.
% * Ridgway GR. Statistical analysis for longitudinal MR
%   imaging of dementia. PhD thesis. 2009.
% * Winkler AM, Ridgway GR, Webster MG, Smith SM, Nichols TE.
%   Permutation inference for the general linear model.
%   Neuroimage. 2014 May 15;92:381-97.
% _____________________________________
% A. Winkler, G. Ridgway & T. Nichols
% FMRIB / University of Oxford
% Mar/2012 (1st version)
% Aug/2013 (this version)
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

switch lower(meth),
    case 'guttman'
        idx = any(C~=0,2);
        X   = M(:,idx);
        Z   = M(:,~idx);
        eCm = vertcat(C(idx,:),C(~idx,:));
        
    case 'beckmann'
        C2  = null(C');
        Q   = pinv(M'*M);
        F1  = pinv(C'*Q*C);
        Pc  = C*pinv(C'*Q*C)*C'*Q;
        C3  = C2 - Pc*C2;
        F3  = pinv(C3'*Q*C3);
        X   = M*Q*C*F1;
        Z   = M*Q*C3*F3;
        eCm = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'winkler'
        Q   = pinv(M'*M);
        X   = M*Q*C*pinv(C'*Q*C);
        Z   = (M*Q*M'-X*pinv(X))*Y;
        eCm = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'ridgway'
        X     = M*pinv(C');
        C0    = eye(size(M,2)) - C*pinv(C);
        [Z,~] = svd(M*C0);
        Z     = Z(:,1:rank(M)-rank(C));
        X     = X-Z*(pinv(Z)*X);
        eCm   = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'none'
        X     = M;
        Z     = [];
        eCm   = C;
        
    otherwise
        error('''%s'' - Unknown partitioning scheme',meth);
end
eCx = eCm(1:size(X,2),:);