function tfdestat = palm_tfde(X,y,opts,plm)
% Compute the TFDE statistic, for volume or surface
% data (vertexwise or facewise).
% 
% Usage:
% tfpestat = palm_tfde(X,y,opts,plm)
% 
% Inputs:
% - X    : Statistical map.
% - y    : Modality index (of those stored in the plm struct).
% - opts : Struct with PALM options.
% - plm  : Struct with PALM data.
% 
% Outputs:
% - tfpestat  : TFDE map.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Apr/2015
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

% Choose an appropriate mask struct.
if opts.NPC || opts.MV,
    S = plm.maskinter;
else
    if plm.nmasks == 1,
        S = plm.masks{1};
    else
        S = plm.masks{y};
    end
end

% Inject the data.
mask    = S.data;
D       = double(S.data);
D(mask) = X;

% "delta h"
dh = max(X(:))/100;

if plm.Yisvol(y),

    % Volume (voxelwise data)
    tfdestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        CC     = bwconncomp(D>=h,opts.tfce.conn);
        sizese = cellfun(@numel,CC.PixelIdxList);
        for c = 1:CC.NumObjects,
            tfdestat(CC.PixelIdxList{c}) = sum(D(CC.PixelIdxList{c})) / sizese(c);
        end
    end
    
elseif plm.Yisvtx(y),
    
    % Vertexwise surface data
    tfdestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        dpxl  = palm_vtxlabel(D>=h,plm.srf{y}.data.fac);
        U     = unique(dpxl(dpxl>0))';
        for u = 1:numel(U),
            idx = dpxl == U(u);
            tfdestat(idx) = sum(D(dpxl == U(u))) / sum(plm.Yarea{y}(idx));
        end
    end
    
elseif plm.Yisfac(y),
    
    % Facewise surface data
    tfdestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        dpxl  = palm_faclabel(D>=h,plm.srf{y}.data.fac);
        U     = unique(dpxl(dpxl>0))';
        for u = 1:numel(U),
            idx = dpxl == U(u);
            tfdestat(idx) = sum(D(dpxl == U(u))) / sum(plm.Yarea{y}(idx));
        end
    end
end

% Return as a vector with the same size as X.
tfdestat = tfdestat(mask);
tfdestat = tfdestat(:)';