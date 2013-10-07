function [maxsize,clstat,sizes] = palm_clusterm(X,y,thr,opts,plm)

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

% Inject the data and threshold it.
mask    = S.data;
D       = double(S.data);
D(mask) = X;
Dt      = D >= thr;

% Do the labelling
if plm.Yisvol(y),
    
    % Volume (voxelwise data)
    % bwconncomp is slightly faster and
    % less memory intensive than bwlabel
    CC = bwconncomp(Dt);

elseif plm.Yisvtx(y),
    
    % Vertexwise surface data
    dpxl = vtxlabel(Dt,plm.srf{y}.data.fac);

elseif plm.Yisfac(y),
    
    % Facewise surface data
    dpxl = faclabel(Dt,plm.srf{y}.data.fac);
    
end

% Compute the sizes and the statistic
if plm.Yisvol(y),
    
    % A loop here is 4x faster than cellfun
    sizes = zeros(CC.NumObjects,1);
    for c = 1:CC.NumObjects,
        sizes(c) = sum(D(CC.PixelIdxList{c}));
    end
    
    % Compute the statistic image (this should be for the 1st perm only)
    if nargout > 1,
        clstat = zeros(size(D));
        for c = 1:CC.NumObjects,
            clstat(CC.PixelIdxList{c}) = sizes(c);
        end
        clstat = clstat(mask)';
    end
    
elseif plm.Yisvtx(y) || plm.Yisfac(y),
    
    % Now compute the cluster stats
    U     = unique(dpxl(dpxl>0))';
    sizes = zeros(size(U));
    for u = 1:numel(U),
        sizes(u) = sum(D(dpxl == U(u)));
    end
    
    % Compute the statistic image (this is normally for the 1st perm only)
    if nargout > 1,
        clstat = zeros(size(D));
        for u = 1:numel(U),
            clstat(dpxl == U(u)) = sizes(u);
        end
        clstat = clstat(mask)';
    end
end

% In fact, only the max matters, because the uncorrected cluster extent
% doesn't make much sense (probably TFCE doesn't make sense either for
% the same reason...). Make sure the output isn't empty.
if isempty(sizes),
    sizes = 0;
end
maxsize = max(sizes);