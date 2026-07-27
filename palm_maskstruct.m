function S = palm_maskstruct(mask,readwith,extra,affine,size)
% Create a struct for a mask, as if it had been read from a file.
% This is useful to save the data later.
%
% Usage:
% M = palm_maskstruct(mask,readwith,extra)
%
% Inputs:
% mask     : A (1 by m) real array.
% readwith : A string telling which function was used to read
%            original data. See 'palm_miscread.m' for help.
% extra    : A struct that varies according to which function
%            was used to read the data.
% affine   : Affine matrix, used for testing sizes.
% size     : Size of the original data
%
% Usage:
% S        : A struct derived from the 'extra' argument along
%            The mask itself will then be in M.data.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013 (first version)
% Mar/2024 (this version)
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

% This is common to all cases below.
S.filename = '';
S.readwith = readwith;

% If more information is given, use to create the structs
% that will be used later to mask and save files. Note that
% the field 'filename' should remain empty.
switch lower(readwith)
    
    case {'load','csvread','vestread'}
        
        % If the original data is a CSV or VEST file.
        S.data  = mask;
        S.extra = extra;

    case 'nifticlass'

        % If the original data is NIFTI and was read with the NIFTI class
        S.data          = palm_conv2toN(mask,extra.dat.dim(1:3),4);
        S.extra.mat     = extra.mat;
        S.extra.dat.dim = extra.dat.dim;

    case 'fs_load_nifti'

        % If the original data is NIFTI and was read with FreeSurfer.
        S.data                 = palm_conv2toN(mask,extra.hdr.dim(2:4),4);
        S.extra                = extra;
        S.extra.hdr.scl_slope  = 1;
        S.extra.hdr.dim([1 5]) = [3 1];
        S.extra.hdr.pixdim(5)  = 0;
        S.extra.hdr.datatype   = 64;
        S.extra.hdr.bitpix     = 64;

    case 'ipt'

        % If the original data is NIFTI and was read with the MATLAB's
        % Image Processing Toolbox or Octave's equivalent commands
        S.data      = palm_conv2toN(mask,extra.hdr.ImageSize(1:3),4);
        S.extra.hdr = extra.hdr;

    case 'parquet'

        % If the original data is an Apache Parquet file
        S.data  = mask;
        S.extra = extra;

    case 'hdf5read'

        % If the original data is an HDF5 file
        S.data  = mask;
        S.extra = extra;

    case 'matlab'

        % If the original data is in a MATLAB .mat file
        S.data  = mask;
        S.extra = extra;

    case {'fs_read_curv','dpxread'}
        
        % If the original data is an FS curvature.
        S.data  = mask;
        S.extra = extra;
        
    case 'fs_load_mgh'
        
        % If the original data is an FS MGH/MGZ file.
        S.data  = palm_conv2toN(mask,extra.volsz(1:3));
        S.extra = extra;
 
    case 'wb_command'

        % If the original data is CIFTI and was read after convering via
        % via wb_command
        S.data  = mask;
        S.extra = extra;

    case 'cifti-matlab'

        % If the original data is CIFTI and was read with
        % the cifti-matlab toolbox
        nD = numel(extra.diminfo);
        Dlength = ones(1,nD);
        for d = 1:nD
            Dlength(d) = extra.diminfo{d}.length;
        end
        S.data  = palm_conv2toN(mask,Dlength(1:end-1));
        S.extra = extra;
        S.extra.diminfo{end}.length = 1;
        S.extra.diminfo{end}.maps = extra.diminfo{end}.maps(1);
        S = palm_dimreorder(S,true);

    case 'gifti'
        
        % If the original data is a GIFTI file.
        S.data  = mask;
        S.extra = extra;
        S.extra.data = S.extra.data(1);

    otherwise
        error('Unknown format: %s',readwith);
end
S.affine = affine;
S.size   = size;