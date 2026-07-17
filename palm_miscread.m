function X = palm_miscread(filespec,varargin)
% Read various scalar data formats based on the file extension.
%
% X = palm_miscread(filespec,useniiclass,precision,mz3surf);
%
% filespec    : Path to file to be read or specification of an HDF5
%               datablock (see palm_hdf5spec.m).
% useniiclass : True/False. For NIFTI files, use the NIFTI class
%               when reading them. It requires less memory.
% precision   : Ensure the output data is 'single' or 'double'
%               precision.
% mz3surf     : True/False. For MZ3 files, keep as data the surface
%               geometry (if present in the file) as opposed to the
%               scalars).
%
% X is a struct that contains the fields:
% X.filename  : Contains the name of the file.
% X.readwith  : This tells which program or function was used
%               to read the data. This is useful when saving the
%               data back, to use a compatible function.
% X.data      : Array with the actual data. The size can vary
%               according to what was read.
% X.affine    : Affine matrix, to be used only for information hence
%               here in a consistent place for different formats.
%               The affine matrix that matters when saving the data is
%               the one inside extras.
% X.extra     : Contain extra information, depending on the kind
%               of data that was read and the function or
%               program used for reading.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013 (first version)
% Jul/2026 (this version)
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

% Defaults (from palm_defaults.m, but can be overridden)
optsx       = palm_defaults;
useniiclass = optsx.useniiclass;
precision   = optsx.precision;
mz3surf     = optsx.mz3surf;

% Parse arguments
narginchk(1,4);
nA = numel(varargin);
if nA >= 1, useniiclass = varargin{1}; end
if nA >= 2, precision   = varargin{2}; end
if nA >= 3, mz3surf     = varargin{3}; end
if ~ischar(filespec) && ~isstring(filespec)
    error('Input must be a string');
end

% If the filename has wildcards, verify that it resolves to a unique name
if contains(filespec,'*') || contains(filespec,'?')
    filelist = dir(filespec);
    if isscalar(filelist)
        filespec = filelist(1).name;
    elseif numel(filelist) == 0
        error('File not found: %s',filespec);
    else
        error('More than one file match: %s',filespec);
    end
end

% Figure out the file extension and check if file exists
[filename,datapath,permdim] = palm_filespec(filespec);
[~,fnam,fext] = fileparts(filename);
fext = palm_tokenize(strcat(fnam,fext),'.');
if ~ exist(filename,'file')
    error('File not found: %s',filename);
end

% A small exception since .mat can be a MATLAB workspace or an
% FSL VEST file. The logic breaks slightly here
if strcmpi(fext{end},'mat') && ischar(datapath)
    if contains(datapath,'/')
        error('Variable names for MATLAB files must not contain the symbol "/".')
    else
        fext = {'matlab'};
    end
end

% Store the filename, in case there is a need to overwrite this later
X.filename = filename;

% Some formats use external i/o functions that use random numbers. Save
% current state of the random number generator, then restore at the end.
if palm_isoctave
    state = rand('state'); %#ok
else
    state = rng;
end

% Check for external programs
ext = palm_checkprogs;

% For each file type, act accordingly
switch lower(fext{end})

    case 'txt'

        % Read a generic text file
        X.readwith = 'textscan';
        fid        = fopen(X.filename);
        X.data     = textscan(fid,'%s');
        X.data     = X.data{1};
        fclose(fid);
        X.extra    = [];
        X.affine   = NaN;
        X.size     = size(X.data);

    case 'csv'

        % Read a CSV file. It has to contain numeric values only.
        % The command 'csvwrite' is a frontend to 'dlmwrite', which calls
        % 'textscan', which on its turn has a limitation of 100k columns.
        % Using 'load' bypass this issue.
        X.readwith = 'load';
        X.data     = load(X.filename);
        X.extra    = [];
        X.affine   = NaN;
        X.size     = size(X.data);

    case {'mat','con','fts','grp'} % note the 'mat' exception above

        % Read an FSL "VEST" file.
        X.readwith = 'vestread';
        [X.data,X.extra.PPH] = palm_vestread(X.filename);
        X.affine   = NaN;
        X.size     = size(X.data);

    case 'mset'

        % Set of matrices
        X.readwith = 'mset';
        X.data     = palm_msetread(X.filename);
        X.affine   = NaN;
        X.size     = size(X.data);

    case optsx.hdf5

        % HDF5 files
        X.readwith       = 'hdf5read';
        X.data           = palm_hdf5read(filename,datapath);
        X.extra.datapath = datapath;
        X.extra.permdim  = permdim;
        if ~isnan(permdim) && ndims(X.data) < permdim
            error('Dataset "%s" in file "%s" has %d dimension(s), but you asked to permute dimension %d.', ...
                X.extra.datapath,filename,ndims(X.data),X.extra.permdim);
        end
        X.affine = NaN;
        X.size   = size(X.data);

    case 'matlab'

        % MATLAB workspace files
        X.readwith      = 'matlab';
        X.data          = load(filename,datapath);
        X.data          = X.data.(datapath);
        X.extra.varname = datapath;
        X.extra.permdim = permdim;
        X.extra.version = matversion(filename);
        if ~isnan(permdim) && ndims(X.data) < permdim
            error('Variable "%s" in file "%s" has %d dimension(s), but you asked to permute dimension %d.', ...
                X.extra.varname,filename,ndims(X.data),X.extra.permdim);
        end
        X.affine = NaN;
        X.size   = size(X.data);

    case 'parquet'

        % Apache Parquet files
        X.readwith = 'parquet';
        [X.data,X.extra.VariableNames] = palm_parquetread(X.filename);
        X.affine   = NaN;
        X.size     = size(X.data);

    case 'gz'

        % Handle (or not) a gzipped NIFTI or CIFTI file.
        if strcmpi(fext{end-1},'nii')

            if any(strcmpi(fext{end-2},optsx.ciftitypes))

                % Until CIFTI migrates to HDF5, users will have to uncompress manually.
                error('CIFTI files must be uncompressed before they can be read. Use gunzip and try again.');

            else
                % Read as NIFTI proper (not CIFTI)
                if useniiclass
                    error([
                        'Reading of gzipped NIFTI files (.nii.gz) is currently disabled\n' ...
                        'If you are sure that your gzipped files, once uncompressed, are not\n' ...
                        'too large to exceed memory limits, you can include the option ''-noniiclass''\n' ...
                        'in the command line. Otherwise, uncompress manually and try again using\n' ...
                        'as input the .nii files instead.\n' ...
                        'File: %s'],X.filename);
                else
                    if ext.matlab_ipt
                        X.readwith  = 'ipt';
                        X.extra.hdr = niftiinfo(X.filename);
                        X.data      = niftiread(X.filename);
                        X.affine    = X.extra.hdr.Transform.T';
                        X.size      = size(X.data);
                    else
                        X.readwith  = 'fs_load_nifti';
                        X.extra.hdr = load_nifti(X.filename);
                        X.data      = X.extra.hdr.vol;
                        X.extra.hdr.vol = [];
                        if X.extra.hdr.qform_code > 0 % qform visited first
                            X.affine = X.extra.hdr.qform;
                        end
                        if X.extra.hdr.sform_code > 0 % but sform will prevail
                            X.affine = X.extra.hdr.sform;
                        end
                        X.size = size(X.data);
                    end
                end
            end
        else
            error('Unrecognised format with extension %s%s',fext0,fext);
        end

    case {'nii','hdr','img'}

        % Handle NIFTI and CIFTI files.
        if strcmpi(fext{end},'nii') && any(strcmpi(fext{end-1},optsx.ciftitypes))

            % Read a CIFTI file
            X.readwith = 'cifti-matlab';
            tmp = cifti_read(X.filename);
            X.data = tmp.cdata;
            X.extra.diminfo  = tmp.diminfo;
            X.extra.metadata = tmp.metadata;
            X.extra.cifti_file_extension = fext{end-1};
            X.affine = NaN;
            X.size = size(X.data);

        else
            % Read a NIFTI file. Note that this will should not
            % be used for ANALYZE.
            if useniiclass
                X.readwith = 'nifticlass';
                X.extra    = nifti(X.filename);
                X.data     = X.extra.dat;
                X.affine   = X.extra.mat;
                X.size     = size(X.data);
            else
                if ext.matlab_ipt
                    X.readwith  = 'ipt';
                    X.extra.hdr = niftiinfo(X.filename);
                    X.data      = niftiread(X.filename);
                    X.affine    = X.extra.hdr.Transform.T';
                    X.size      = size(X.data);
                else
                    X.readwith  = 'fs_load_nifti';
                    X.extra.hdr = load_nifti(X.filename);
                    X.data      = X.extra.hdr.vol;
                    X.extra.hdr.vol = [];
                    if X.extra.hdr.qform_code > 0 % qform visited first
                        X.affine = X.extra.hdr.qform;
                    end
                    if X.extra.hdr.sform_code > 0 % but sform will prevail
                        X.affine = X.extra.hdr.sform;
                    end
                    X.size = size(X.data);
                end
            end
        end

    case {'dpv','dpf','dpx','asc'}

        % Read a DPV/DPF file, in ASCII
        X.readwith = 'dpxread';
        [X.data,X.extra.crd,X.extra.idx] = palm_dpxread(X.filename);
        X.affine   = NaN;
        X.size     = size(X.data);

    case 'srf'

        % Read a SRF file, in ASCII
        X.readwith = 'srfread';
        [X.data.vtx,X.data.fac] = palm_srfread(X.filename);
        X.affine   = NaN;
        X.size     = NaN;

    case 'obj'

        % Read a Wavefront file
        X.readwith = 'wavefront';
        [X.data.vtx,X.data.fac,X.extra] = palm_objread(X.filename);
        X.affine   = NaN;
        X.size     = NaN;

    case 'mz3'

        % Read a MZ3 file
        X.readwith = 'mz3';
        [vtx,fac,colour] = readMz3(X.filename);
        if mz3surf
            X.data.vtx = vtx;
            X.data.fac = fac;
            X.extra.colour = colour;
            X.size     = NaN;
        else
            X.data      = colour;
            X.extra.vtx = vtx;
            X.extra.fac = fac;
            X.size      = size(X.data);
        end
        X.affine   = NaN;

    case optsx.fscurv

        % Read a FreeSurfer curvature file
        X.readwith = 'fs_read_curv';
        [X.data,X.extra.fnum] = read_curv(X.filename);
        X.affine   = NaN;
        X.size     = size(X.data);

    case optsx.fssurf

        % Read a FreeSurfer surface file
        X.readwith = 'fs_read_surf';
        [X.data.vtx,X.data.fac] = read_surf(X.filename);
        X.data.fac = X.data.fac + 1;
        X.affine   = NaN;
        X.size     = NaN;

    case {'mgh','mgz'}

        % Read a FreeSurfer MGH/MGZ file
        X.readwith = 'fs_load_mgh';
        [X.data,X.extra.M,X.extra.mr_parms,X.extra.volsz] = load_mgh(X.filename);
        X.affine   = X.extra.M;
        X.size     = size(X.data);

    case 'annot'

        % Read a FreeSurfer annotation file
        X.readwith = 'fs_load_annot';
        [X.extra.vertices,X.extra.label,X.extra.colourtab] = read_annotation(X.filename);

        % For each structure, replace its label by its index, which
        % is the actual label
        X.data = X.extra.label;
        for s = 1:X.extra.colourtab.numEntries
            X.data(X.extra.label == X.extra.colourtab.table(s,5)) = s;
        end
        X.data(X.data == 0) = 1;
        X.affine = NaN;
        X.size   = size(X.data);

        % Create a Matlab colourmap, useful to make figures
        X.extra.colourmap = X.extra.colourtab.table(:,1:3)/255;

    case 'gii'

        % Read a GIFTI file (no mapped file arrays)
        X.readwith = 'gifti';
        gii = gifti(X.filename);
        if isfield(gii,'vertices') && isfield(gii,'faces')
            X.data.vtx = gii.vertices;
            X.data.fac = gii.faces;
            if isfield(gii,'mat')
                vtx = [X.data.vtx ones(size(X.data.vtx,1),1)];
                vtx = vtx * gii.mat;
                X.data.vtx = vtx(:,1:3);
                X.extra.mat = gii.mat';
                X.affine = X.extra.mat;
            else
                X.affine = NaN;
            end
            X.size = NaN;
        elseif isfield(gii,'cdata')
            X.data = gii.cdata';
            if isfield(gii,'mat')
                X.extra.mat = gii.mat';
                X.affine = X.extra.mat;
            else
                X.affine = NaN;
            end
            X.size = size(X.data);
        else
            error('Invalid GIFTI file: %s',X.filename);
        end
        X.extra.gifti = gii;

    otherwise
        error('File extension %s not known. Data cannot be loaded\n',fext{end});
end

% Restore the state of the random number generator.
if palm_isoctave
    rand('state',state); %#ok
else
    rng(state);
end

% Enforce a certain precision defined by the user:
if ~ (isstruct(X.data) || iscell(X.data))
    if strcmpi(precision,'double')
        X.data = double(X.data);
    elseif strcmpi(precision,'single')
        X.data = single(X.data);
    end
end

% ==============================================================
function result = contains(str,ch)
% Test is a character exists in a string
result = any(str == ch);

% ==============================================================
function ver = matversion(filename)
% Determines the version of a MATLAB .mat file

% Read first 128 bytes (header size for Level 5)
fid    = fopen(filename,'r');
header = fread(fid,128,'*char')';
fclose(fid);

% v7.3 (HDF5)
sig = 'MATLAB 7.3 MAT-file';
if strcmp(header(1:length(sig)),sig)
    ver = '-v7.3';
    return;
end

% Level 5 files (v6 / v7)
sig = 'MATLAB 5.0 MAT-file';
if strcmp(header(1:length(sig)),sig)
    % Distinguishing v6 vs v7 isn't always possible from header alone
    % (both use Level 5 format; v7 adds compression + Unicode).
    % Will use v7 as its the current format
    ver = '-v7';  % change from -v7 to -v6 to drop compression and use in older MATLABs
    return;
end

% v4 v4 files do not have the "MATLAB 5.0" header.
% Check first 4 bytes (MOPT) for valid v4 signature.
fid  = fopen(filename,'r');
mopt = fread(fid,4,'uint8');
fclose(fid);

% Valid MOPT for v4: first byte 0-4, second=0, etc.
if mopt(1) <= 4 && mopt(2) == 0 && mopt(3) <= 5 && mopt(4) <= 2
    % Additional heuristic: v4 often has printable chars early on
    fid   = fopen(filename,'r');
    early = fread(fid,32,'*char')';
    fclose(fid);
    if any(early >= 32 & early <= 126) || sum(mopt == 0) > 0
        ver = 'v4';
        return;
    end
end

% Not a recognized .mat file
ver = NaN;
