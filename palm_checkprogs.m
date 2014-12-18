function ext = palm_checkprogs
% Test whether some external programs or toolboxes
% are available, namely, FSL, FreeSurfer, SPM and
% Jimmy Shen's NIFTI toolbox.
%
% ext = checkprogs
% 
% 'ext' is a struct containing one field for each
% of these applications, each being containing
% a 0 (false) or 1 (true) depending on whether
% these programs are available or not.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

persistent palm_extern;
if isempty(palm_extern),
    
    % Check the path of PALM and add the path for the NIFTI class.
    palm_extern.palmpath = fileparts(mfilename('fullpath'));
    addpath(fullfile(palm_extern.palmpath,'niftimatlib','matlab'));
    addpath(fullfile(palm_extern.palmpath,'niftimatlib','matlab'));
        
    % Check for FSL
    palm_extern.fsl = false;
    fsldir = getenv('FSLDIR');
    if ~isempty(fsldir),
        palm_extern.fsl = true;
        addpath(fullfile(fsldir,'etc','matlab'));
        fprintf('Found FSL in %s\n',fsldir);
    end
    
    % Check for FreeSurfer
    palm_extern.fs  = false;
    fshome = getenv('FREESURFER_HOME');
    if ~isempty(fshome),
        palm_extern.fs = true;
        addpath(fullfile(fshome,'matlab'));
        fprintf('Found FreeSurfer in %s\n',fshome);
    end
    
    % Check for the NIFTI toolbox
    palm_extern.nii = false;
    if exist('load_nii') == 2 && exist('save_nii') == 2, %#ok
        palm_extern.nii = true;
    end
    
    % Check for SPM
    palm_extern.spm = false;
    try %#ok
        spm_check_installation('basic');
        palm_extern.spm = true;
        spmpath = which('spm');
        spmpath = fileparts(which('spm'));
        fprintf('Found SPM in %s\n',spmpath);
    end
end
ext = palm_extern;
