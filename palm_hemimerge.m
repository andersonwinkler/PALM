function palm_hemimerge(varargin)
% Merge data for lh and rh (left and right hemispheres) into a
% single file. Various input files can be given, starting with
% 'lh' and/or 'rh'. The script finds the correct pair, removes
% repeated files, and merges them into files beginning with 'bh'.
%
% palm_hemimerge <files>
%
% Wildcards are accepted.
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jun/2016
% http://brainder.org

% List of files (with wildcards)
j = 1;
for fl = 1:nargin
    F = dir(varargin{fl});
    for f = numel(F):-1:1
        if F(f).name(1) == '.'
            F(f) = [];
        else
            Flist{j} = F(f).name;
            j = j + 1;
        end
    end
end
Flist = flipud(Flist');

% Prepare pairs to merge and output names
F = cell(0,3);
f = 1;
for fl = numel(Flist):-1:1
    if strcmpi(Flist{fl}(1:2),'lh')
        if ~ any(strcmpi(Flist{fl},F(:,1)))
            F{f,1} = Flist{fl};
            F{f,2} = strcat('rh',Flist{fl}(3:end));
            F{f,3} = strcat('bh',Flist{fl}(3:end));
            f = f + 1;
        end
    elseif strcmpi(Flist{fl}(1:2),'rh')
        if ~ any(strcmpi(Flist{fl},F(:,2)))
            F{f,1} = strcat('lh',Flist{fl}(3:end));
            F{f,2} = Flist{fl};
            F{f,3} = strcat('bh',Flist{fl}(3:end));
            f = f + 1;
        end
    elseif strfind(Flist{fl},'_hemi-L')
        if ~ any(strcmpi(Flist{fl},F(:,1)))
            F{f,1} = Flist{fl};
            F{f,2} = strrep(Flist{fl},'_hemi-L','_hemi-R');
            F{f,3} = strrep(Flist{fl},'_hemi-L','_hemi-B');
            f = f + 1;
        end
    elseif strfind(Flist{fl},'_hemi-R')
        if ~ any(strcmpi(Flist{fl},F(:,2)))
            F{f,1} = strrep(Flist{fl},'_hemi-R','_hemi-L');
            F{f,2} = Flist{fl};
            F{f,3} = strrep(Flist{fl},'_hemi-R','_hemi-B');
            f = f + 1;
        end
    else
        warning([...
            'Files for merger must start with "lh" or "rh", \n' ...
            'or follow the BIDS convention and contain "hemi-L" or "hemi-R".\n' ...
            'Skipping: %s\n'], ...
            Flist{fl});
    end
end

% For each valid pair
for f = 1:size(F,1)
    
    fprintf('Working on: %s and %s\n',F{f,1},F{f,2});
    
    % Load L and R
    L = palm_miscread(F{f,1});
    R = palm_miscread(F{f,2});
    if ~ strcmpi(L.readwith,R.readwith)
        error('Left and right must be of the same type and format');
    end
    
    % Prepare data to save
    B = L;
    switch R.readwith
        case {'load','csvread','fs_load_mgh'}
            
            dim = 1;
            checkdim(L.data,R.data,dim);
            B.data = cat(dim,L.data,R.data);
            
        case 'dpxread'

            dim = 1;
            checkdim(L.data,R.data,dim);
            B.data = cat(dim,L.data,R.data);
            B.extra.crd = cat(dim,L.extra.crd,R.extra.crd);
            B.extra.idx = (0:size(B.data,1)-1)';
            
        case {'srfread','fs_read_surf'}

            dim = 1;
            checkdim(L.data,R.data,dim);
            B.data.vtx = cat(dim,L.data.vtx,R.data.vtx);
            B.data.fac = cat(dim,L.data.fac,R.data.fac + size(L.data.vtx,1));
            
        case 'fs_read_curv'

            dim = 1;
            checkdim(L.data,R.data,dim);
            B.data = cat(dim,L.data,R.data);
            B.extra.fnum = L.extra.fnum + R.extra.fnum;

        case 'gifti'

            if isfield(B.data,'vtx') && isfield(B.data,'fac')
                dim = 1;
                checkdim(L.data.vtx,R.data.vtx,dim);
                checkdim(L.data.fac,R.data.fac,dim);
                B.data.vtx = cat(dim,L.data.vtx,R.data.vtx);
                B.data.fac = cat(dim,L.data.fac,R.data.fac + size(L.data.vtx,1));
            else
                dim = 2;
                checkdim(L,R,dim);
                B.data = cat(dim,L.data,R.data);
            end
            for d = 1:numel(B.extra.data)
                B.extra.data{d}.attributes.Dim(1) = L.extra.data{d}.attributes.Dim(1) + R.extra.data{d}.attributes.Dim(1);
            end
            
        otherwise
            warning('Cannot deal with files read with %s. Skipping: ?%s\n',L.readwith,Flist{fl}(2:end));
    end
    
    % Save
    B.filename = F{f,3};
    if any(strcmpi(B.filename(end-3:end),{'.mgh','.mgz'}))
        B.filename = B.filename(1:end-4);
    end
    palm_miscwrite(B);
end

function checkdim(L,R,dim)
if size(L,dim) ~= size(R,dim)
    error('Left and Right sides must have the same number of elements (e.g., same number of vertices)')
end