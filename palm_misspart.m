function [X,Z,eCm,eCx,Y,Pidx,isdiscrete] = palm_misspart(M,C,meth,my,mm,mcar)
% Partition the model for missing data, generating all the sub-models
% that can later be subjected to NPC.
%
% Usage:
% [X,Z,eCm,eCx,Y,Pidx,isdiscrete] = palm_misspart(M,C,meth,my,mm)
%
% Inputs:
% M    : Design matrix, to be partitioned.
% C    : Contrast that will define the partitioning.
% meth : Method for the partitioning. It can be any of those available
%        in palm_partition.
% my   : Missing data indicators for the observations (Y).
% mm   : Missing data indicators for the design (M).
% mcar : Boolean, indicating whether the missing data process is completely
%        at random (true) or not (false).
%
% Outputs:
% X    : Cell array with sets of EVs of interest.
% Z    : Cell array with sets of nuisance EVs.
% eCm  : Cell array of effective contrasts.
% eCx  : Same as above, but considering only X.
% Y    : Cell array with indices (logical) or data for regression (double).
%        If empty, it's equivalent to a vector index full of ones.
% Pidx : Cell array of indices to select the rows of the permutation
%        matrices.
% isdiscrete : Vector indicating if the respective cell array contains both
%        discrete (binary) X and Y, such that the Chi^2 (Yates) test can be
%        be performed.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2016 Anderson M. Winkler
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

% Initial partitioning of the model:
[x,z,ecm,ecx] = palm_partition(M,C,'guttman');

% Partitioning to obtain the indices:
mx = zeros(size(mm,1),size(C,2));
for t = 1:size(C,2),
    tmp = palm_partition(mm,C(:,t),'guttman');
    mx(:,t) = any(tmp,2);
end
mx = all(mx,2);
[~,mz] = palm_partition(mm,C,'guttman');

% These are (logical) indices for the variables that are available:
iy = ~ my;
ix = ~ any(mx,2);
iz = ~ any(mz,2); % if mz is empty, iz is also empty

% These are the actual missing indicators (double) that go in the design:
my = double(my);
mx = double(mx);
mz = double(mz);

% It's simpler and faster to fork the code 16 times than have
% multiple loops and conditions for the up-to-eight equations.
isdiscrete = false;
if isempty(mz),
    if       all(iy) &   all(ix), % Case 1
        Pidx{1} = [];  Y{1} = [];  X{1} = x;  Z{1} = [];  eC{1} = ecm;
    elseif ~ all(iy) &   all(ix), % Case 2
        if mcar,
            Pidx{1} = iy;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
        else
            Pidx{1} = iy;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
            Pidx{2} = [];  Y{2} = my;       X{2} = x;               Z{2} = [];  eC{2} = ecm;
            isdiscrete = [false false];
        end
    elseif   all(iy) & ~ all(ix), % Case 3
        if mcar,
            Pidx{1} = ix;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
        else
            Pidx{1} = ix;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
            Pidx{2} = [];  Y{2} = [];       X{2} = mx;              Z{2} = [];  eC{2} = mkcon(X{2},Z{2});
            isdiscrete = [false false];
        end
    elseif ~ all(iy) & ~ all(ix), % Case 5
        if mcar,
            Pidx{1} = iy & ix;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
        else
            Pidx{1} = iy & ix;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = [];  eC{1} = ecm;
            Pidx{2} = ix;       Y{2} = my(Pidx{2});  X{2} = x(Pidx{2},:,:);  Z{2} = [];  eC{2} = ecm;
            Pidx{3} = iy;       Y{3} = Pidx{3};      X{3} = mx(Pidx{3});     Z{3} = [];  eC{3} = ecm;
            Pidx{4} = [];       Y{4} = my;           X{4} = mx;              Z{4} = [];  eC{4} = ecm; % discrete
            isdiscrete = [false false false true];
        end
    end
else
    if       all(iy) &   all(ix) &   all(iz), % Case 1
        Pidx{1} = [];  Y{1} = [];  X{1} = x;  Z{1} = z;  eC{1} = ecm;
    elseif ~ all(iy) &   all(ix) &   all(iz), % Case 2
        if mcar,
            Pidx{1} = iy;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = iy;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = [];  Y{2} = my;       X{2} = x;               Z{2} = z;               eC{2} = ecm;
            isdiscrete = [false false];
        end
    elseif   all(iy) & ~ all(ix) &   all(iz), % Case 3
        if mcar,
            Pidx{1} = ix;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = ix;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = [];  Y{2} = [];       X{2} = mx;              Z{2} = z;               eC{2} = mkcon(X{2},Z{2});
            isdiscrete = [false false];
        end
    elseif   all(iy) &   all(ix) & ~ all(iz), % Case 4
        if mcar,
            Pidx{1} = iz;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = iz;  Y{1} = Pidx{1};  X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = [];  Y{2} = [];       X{2} = x;               Z{2} = mz;              eC{2} = mkcon(ecx,Z{2});
            isdiscrete = [false false];
        end
    elseif ~ all(iy) & ~ all(ix) &   all(iz), % Case 5
        if mcar,
            Pidx{1} = iy & ix;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = iy & ix;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = ix;       Y{2} = my(Pidx{2});  X{2} = x(Pidx{2},:,:);  Z{2} = z(Pidx{2},:,:);  eC{2} = ecm;
            Pidx{3} = iy;       Y{3} = Pidx{3};      X{3} = mx(Pidx{3});     Z{3} = z(Pidx{3},:,:);  eC{3} = mkcon(X{3},Z{3});
            Pidx{4} = [];       Y{4} = my;           X{4} = mx;              Z{4} = z;               eC{4} = mkcon(X{4},Z{4}); % discrete
            isdiscrete = [false false false true];
        end
    elseif ~ all(iy) &   all(ix) & ~ all(iz), % Case 6
        if mcar,
            Pidx{1} = iy & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = iy & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = iz;       Y{2} = my(Pidx{2});  X{2} = x(Pidx{2},:,:);  Z{2} = z(Pidx{2},:,:);  eC{2} = ecm;
            Pidx{3} = iy;       Y{3} = Pidx{3};      X{3} = x(Pidx{3},:,:);  Z{3} = mz(Pidx{3});     eC{3} = mkcon(ecx,Z{3});
            Pidx{4} = [];       Y{4} = my;           X{4} = x;               Z{4} = mz;              eC{4} = mkcon(ecx,Z{4});
            isdiscrete = [false false false false];
        end
    elseif   all(iy) & ~ all(ix) & ~ all(iz), % Case 7
        if mcar,
            Pidx{1} = ix & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = ix & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = iz;       Y{2} = Pidx{2};      X{2} = mx(Pidx{2});     Z{2} = z(Pidx{2},:,:);  eC{2} = mkcon(X{2},Z{2});
            Pidx{3} = ix;       Y{3} = Pidx{3};      X{3} = x(Pidx{3},:,:);  Z{3} = mz(Pidx{3});     eC{3} = mkcon(ecx,Z{3});
            Pidx{4} = [];       Y{4} = [];           X{4} = mx;              Z{4} = mz;              eC{4} = mkcon(X{4},Z{4});
            isdiscrete = [false false false false];
        end
    elseif ~ all(iy) & ~ all(ix) & ~ all(iz), % Case 8
        if mcar,
            Pidx{1} = iy & ix & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
        else
            Pidx{1} = iy & ix & iz;  Y{1} = Pidx{1};      X{1} = x(Pidx{1},:,:);  Z{1} = z(Pidx{1},:,:);  eC{1} = ecm;
            Pidx{2} = ix & iz;       Y{2} = my(Pidx{2});  X{2} = x(Pidx{2},:,:);  Z{2} = z(Pidx{2},:,:);  eC{2} = ecm;
            Pidx{3} = iy & iz;       Y{3} = Pidx{3};      X{3} = mx(Pidx{3});     Z{3} = z(Pidx{3},:,:);  eC{3} = mkcon(X{3},Z{3});
            Pidx{4} = iy & ix;       Y{4} = Pidx{4};      X{4} = x(Pidx{4},:,:);  Z{4} = mz(Pidx{4});     eC{4} = mkcon(ecx,Z{4});
            Pidx{5} = iz;            Y{5} = my(Pidx{5});  X{5} = mx(Pidx{5});     Z{5} = z(Pidx{5},:,:);  eC{5} = mkcon(X{5},Z{5}); % discrete?
            Pidx{6} = ix;            Y{6} = my(Pidx{6});  X{6} = x(Pidx{6},:,:);  Z{6} = mz(Pidx{6});     eC{6} = mkcon(ecx,Z{6});
            Pidx{7} = iy;            Y{7} = Pidx{7};      X{7} = mx(Pidx{7});     Z{7} = mz(Pidx{7});     eC{7} = mkcon(X{7},Z{7});
            Pidx{8} = [];            Y{8} = my;           X{8} = mx;              Z{8} = mz;              eC{8} = mkcon(X{8},Z{8}); % discrete
            isdiscrete = [false false false false false false false true];
        end
    end
end

for o = 1:numel(Y),
    if ~isempty(Y{o}) && ~islogical(Y{o}),
        Y{o} = bsxfun(@minus,Y{o},mean(Y{o},1));
    end
    X{o} = bsxfun(@minus,X{o},mean(X{o},1));
    Z{o} = bsxfun(@minus,Z{o},mean(Z{o},1));
end

% PCA of Z:
for o = 1:numel(Z),
    sX = size(X{o},2);
    [u,s,~] = svd(Z{o},'econ');
    ds = diag(s);
    tol = 10 * max(size(Z{o})) * eps(max(ds));
    idx = ds > tol;
    Z{o} = u(:,idx)*s(idx,idx);
    eC{o} = [eC{o}(1:sX,:); zeros(size(Z{o},2),size(eC{o},2))];
end

% Partition each of these models using the method indicated by the user:
eCm = cell(size(X));
eCx = eCm;
for o = 1:numel(Y),
    [X{o},Z{o},eCm{o},eCx{o}] = palm_partition(horzcat(X{o},Z{o}),eC{o},meth);
end

% Remove bits that are all zeroes (this needs to be adapted for voxelwise):
for o = numel(Y):-1:1,
    if  (~isempty(Y{o}) && ~islogical(Y{o}) && all(Y{o} == 0,1)) || all(X{o}(:,1,1) == 0,1),
        Pidx(o) = []; Y(o) = []; X(o) = []; Z(o) = [];
        eCm(o) = []; eCx(o) = [];
        isdiscrete(o) = [];
    end
end

% ==============================================================
function eC = mkcon(X,Z);
% Shortcut to create the effective contrast.
if isempty(Z) || size(X,1) == size(Z,1),
    % Typical case, X is X and Z is Z.
    eC = vertcat(eye(size(X,2)),zeros(size(Z,2),size(X,2)));
else
    % Here X is C, and Z is Z.
    eC = vertcat(X,zeros(size(Z,2),size(X,2)));
end