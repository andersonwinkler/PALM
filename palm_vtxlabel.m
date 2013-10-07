function dpvl = palm_vtxlabel(dpv,fac)

% For the labelled data
dpvl = zeros(size(dpv));

% Identify faces that are entirely (3 vertices)
% within the supra-trreshold vertices. These are
% the only that will be labelled.
facm = reshape(dpv(fac(:)),size(fac));
faci = all(facm,2);

% For each face
k = 1; % label indices (to be incremented)
for f = find(faci)',
    flab = sort(dpvl(fac(f,:)));
    if flab(3) == 0,
        dpvl(fac(f,:)) = k;
    else
        if flab(2) ~= 0;
            dpvl(dpvl == flab(2)) = flab(3);
        end
        if flab(1) ~= 0;
            dpvl(dpvl == flab(1)) = flab(3);
        end
        dpvl(fac(f,:)) = flab(3);
    end
    k = k + 1;
end
