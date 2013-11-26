function extras = read_extras(fname)
% Read extra bits of information
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% Id: read_extras.m 2237 2008-09-29 17:39:53Z guillaume 

%
% niftilib $Id: read_extras.m,v 1.3 2012/03/22 18:36:33 fissell Exp $
%



extras = struct;
[pth,nam,ext] = fileparts(fname);
switch ext
case {'.hdr','.img','.nii'}
    mname = fullfile(pth,[nam '.mat']);
case {'.HDR','.IMG','.NII'}
    mname = fullfile(pth,[nam '.MAT']);
otherwise
    mname = fullfile(pth,[nam '.mat']);
end

if spm_existfile(mname),
    try,
        extras = load(mname);
    catch,
        warning('Can not load "%s" as a binary MAT file.\n', mname);
    end;
end;
