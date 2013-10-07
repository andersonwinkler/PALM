function tfcestat = palm_tfce(X,y,opts,plm)

% Choose an appropriate mask struct.
if opts.NPC,
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
    tfcestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        CC    = bwconncomp(D>=h,opts.tfce.conn);
        integ = cellfun(@numel,CC.PixelIdxList).^opts.tfce.E * h^opts.tfce.H;
        for c = 1:CC.NumObjects,
            tfcestat(CC.PixelIdxList{c}) = ...
                tfcestat(CC.PixelIdxList{c}) + integ(c);
        end
    end
    
elseif plm.Yisvtx(y),
    
    % Vertexwise surface data
    tfcestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        dpxl  = vtxlabel(D>=h,plm.srf{y}.data.fac);
        U     = unique(dpxl(dpxl>0))';
        for u = 1:numel(U),
            idx = dpxl == U(u);
            tfcestat(idx) = tfcestat(idx) + ...
                sum(plm.Yarea{y}(idx)).^opts.tfce.E * h^opts.tfce.H;
        end
    end
    
elseif plm.Yisfac(y),
    
    % Facewise surface data
    tfcestat = zeros(size(D));
    for h = dh:dh:max(D(:));
        dpxl  = faclabel(D>=h,plm.srf{y}.data.fac);
        U     = unique(dpxl(dpxl>0))';
        for u = 1:numel(U),
            idx = dpxl == U(u);
            tfcestat(idx) = tfcestat(idx) + ...
                sum(plm.Yarea{y}(idx)).^opts.tfce.E * h^opts.tfce.H;
        end
    end
end

% Return as a vector with the same size as X.
tfcestat = tfcestat(mask)';