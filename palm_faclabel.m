function dpfl = palm_faclabel(dpf,fac)

% Ignore faces that are not in the mask
dpfl = zeros(size(dpf));
facm = bsxfun(@times,fac,dpf);

% For each face
k = 1; % label indices (to be incremented)
for f = find(dpf)',
    
    % Find other faces that share 2 vtx (1 edge)
    neifacidx = ...
        facm == fac(f,1) | ...
        facm == fac(f,2) | ...
        facm == fac(f,3);
    neifacidx = sum(neifacidx,2) >= 2; % only 2+ shared vertices, not 1
    numneigh  = sum(neifacidx); % number of neighbours
    if numneigh > 1, % if not isolated face
        flab = sort(dpfl(neifacidx));
        if flab(numneigh) == 0,
            dpfl(neifacidx) = k;
        else
            for n = numneigh:-1:1,
                if flab(n) ~= 0,
                    dpfl(dpfl == flab(n)) = flab(numneigh);
                end
            end
            dpfl(f) = flab(numneigh);
        end
    else % isolated face
        dpfl(f) = k;
    end
    k = k + 1;
end
