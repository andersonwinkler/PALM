function A = palm_calcarea(S,isvtx)

nV = size(S.vtx,1);
nF = size(S.fac,1);

facvtx = [S.vtx(S.fac(:,1),:) S.vtx(S.fac(:,2),:) S.vtx(S.fac(:,3),:)];
facvtx0(:,1:6) = facvtx(:,1:6) - [facvtx(:,7:9) facvtx(:,7:9)];
cp  = cross(facvtx0(:,1:3),facvtx0(:,4:6),2);
A   = sqrt(sum(cp.^2,2))./2;

if isvtx,
    dpv = zeros(nV,1);
    A   = A/3;
    for f = 1:nF,
        dpv(S.fac(f,:)) = dpv(S.fac(f,:)) + A(f);
    end
    A = dpv;
end