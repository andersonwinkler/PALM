function palm_hemisplit(varargin)
% Split back bh into lh and rh for files that have previously
% been merged with palm_hemimerge.
% Output file names have prefix "lh" and "rh".
%
% palm_hemisplit <files>
%
% Wildcards are accepted.
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jun/2016
% http://brainder.org

% List of files (with wildcards)
j = 1;
for a = 1:nargin,
    F = dir(varargin{a});
    for f = numel(F):-1:1,
        if F(f).name(1) == '.',
            F(f) = [];
        else
            Flist{j} = F(f).name;
            j = j + 1;
        end
    end
end
Flist = flipud(Flist');

% For each input file.
for f = 1:numel(Flist),
    fprintf('Working on: %s\n',Flist{f});
    
    B = palm_miscread(Flist{f});
    L = B; R = B;
    switch B.readwith,
        case {'load','csvread','fs_load_mgh'},
            nV = size(B.data,1);
            nV2 = nV/2;
            L.data = B.data(1:nV2,:,:,:);
            R.data = B.data(nV2+1:end,:,:,:);
            
        case 'dpxread',
            nV  = size(B.data,1);
            nV2 = nV/2;
            L.data      = B.data(1:nV2,:,:,:);
            R.data      = B.data(nV2+1:end,:,:,:);
            L.extra.crd = B.extra.crd(1:nV2,:,:,:);
            R.extra.crd = B.extra.crd(nV2+1:end,:,:,:);
            L.extra.idx = (1:size(L.data,1))';
            R.extra.idx = L.extra.idx;
            
        case {'srfread','fs_read_surf'},
            
            nV = size(B.data.vtx,1);
            nF = size(B.data.fac,1);
            nV2 = nV/2;
            nF2 = nF/2;
            L.data.vtx = B.data.vtx(1:nV2,:);
            R.data.vtx = B.data.vtx(nV2+1:end,:);
            L.data.fac = B.data.fac(1:nF2,:);
            R.data.fac = B.data.fac(nF2+1:end,:)-nV2;
            
        case 'fs_read_curv',
            L.data      = B.data(1:nV2,:,:,:);
            R.data      = B.data(nV2+1:end,:,:,:);
            L.extra.fnum = B.extra.fnum/2;
            R.extra.fnum = L.extra.fnum;
            
        otherwise
            warning('Cannot deal with files read with %s. Skipping: %s\n',B.readwith,Flist{f});
    end
    
    % Adjust filenames and save
    filename = Flist{f};
    if any(strcmpi(filename(end-3:end),{'.mgh','.mgz'})),
        filename = filename(1:end-4);
    end
    if strcmpi(Flist{f}(1:2),'bh'),
        L.filename = filename;
        L.filename(1:2) = 'lh';
        R.filename = filename;
        R.filename(1:2) = 'rh';
    else
        L.filename = strcat('lh_',filename);
        R.filename = strcat('rh_',filename);
    end
    palm_miscwrite(L);
    palm_miscwrite(R);
end