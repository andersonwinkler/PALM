function S = palm_maskstruct(mask,readwith,extra)
% Create a struct for a mask, as if it had been read from a file.
% This is useful to save later the data.
%
% Usage:
% M = maskstruct(mask,readwith,extra)
%
% Inputs:
% Y        : A (n by m) real array.
% readwith : A string telling which function was used to read
%            original data. See 'miscread.m' for help.
% extra    : A struct that varies according to which function
%            was used to read the data.
%
% Usage:
% S        : A struct derived from the 'extra' argument along
%            The mask itself will then be in M.data.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org


% This is common to all cases below.
S.filename = '';
S.readwith = readwith;

% If more information is given, use to create the structs
% that will be used later to mask and save files. Note that
% the field 'filename' should remain empty.
switch lower(readwith),
    
    case {'load','csvread','vestread'},
        
        % If the original data is a CSV or VEST file.
        S.data = mask;
        
    case 'spm_spm_vol'
        
        % If the original data is NIFTI and was read with SPM.
        S.data           = palm_conv2to4(mask,extra(1).dim(1:3));
        S.extra          = extra(1);
        S.extra.dt(1)    = spm_type('float64');
        S.extra.pinfo(1) = 1;
        
    case 'fs_load_nifti',
        
        % If the original data is NIFTI and was read with FreeSurfer.
        S.data                 = palm_conv2to4(mask,extra.hdr.dim(2:4));
        S.extra                = extra;
        S.extra.hdr.scl_slope  = 1;
        S.extra.hdr.dim([1 5]) = [3 1];
        S.extra.hdr.pixdim(5)  = 0;
        S.extra.hdr.datatype   = 64;
        S.extra.hdr.bitpix     = 64;
        
    case 'fsl_read_avw',
        
        % If the original data is NIFTI and was read with FSL.
        S.data  = palm_conv2to4(mask,extra.dims(1:3));
        S.extra = extra;
        if ~ isfield(S.extra,'vtype'),
            S.extra.vtype = 'd';
        end
        
    case 'nii_load_nii',
        
        % If the original data is NIFTI and was read with the NIFTI toolbox.
        S.data                      = palm_conv2to4(mask,extra.hdr.dime.dim(2:4));
        S.extra                     = extra;
        S.extra.hdr.dime.dim([1 5]) = [3 1];
        S.extra.hdr.dime.pixdim(5)  = 0;
        S.extra.hdr.dime.datatype   = 64;
        S.extra.hdr.dime.bitpix     = 64;
        
    case 'fs_read_curv',
        
        % If the original data is an FS curvature.
        S.data  = mask(:);
        S.extra = extra;
        
    case 'fs_load_mgh',
        
        % If the original data is an FS MGH/MGZ file.
        S.data  = palm_conv2to4(mask,extra.volsz(1:3));
        S.extra = extra;
end
