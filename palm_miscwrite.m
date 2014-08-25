function palm_miscwrite(X)
% Write various scalar data formats based on the struct X.
% 
% palm_miscwrite(X);
% 
% See 'palm_miscread.m' for a description of the contents of X.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

switch lower(X.readwith),
    
    case 'textscan',
        
        % Write a generic text file. Note that this doesn't recover
        % the original file (spaces become newlines).
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.txt');
        end
        fid = fopen(X.filename,'w');
        fprintf(fid,'%s\n',X.data{:});
        fclose(fid);
    
    case {'load','csvread'},
        
        % Write a CSV file.
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.csv');
        end
        csvwrite(X.filename,X.data);
        
    case 'vestread',
        
        % Write an FSL "VEST" file
        vestwrite(X.filename,X.data);
        
    case 'nifticlass',
        
        % Write using the NIFTI class.
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.nii');
        end
        dat = file_array(  ...
            X.filename,    ...
            size(X.data),  ...
            'FLOAT32-LE',  ...
            ceil(348/8)*8);
        nii     = nifti;
        nii.dat = dat;
        nii.mat = X.extra.mat;
        create(nii);
        nii.dat(:,:,:) = X.data(:,:,:);
        
    case 'spm_spm_vol',
        
        % Write NIFTI with SPM
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.nii');
        end
        X.extra.fname = X.filename;
        X.extra.dt(1) = spm_type('float32'); % for now, save everything as double
        spm_write_vol(X.extra,X.data);
        
    case 'fs_load_nifti',
        
        % Write NIFTI with FreeSurfer
        X.filename = horzcat(X.filename,'.nii.gz');
        X.extra.hdr.vol = X.data;
        X.extra.hdr.datatype = 64; % for now, save everything as double
        X.extra.hdr.bitpix = 64; % for now, save everything as double
        save_nifti(X.extra.hdr,X.filename);
        
    case 'fsl_read_avw',
        
        % Write NIFTI with FSL
        if ~ isfield(X.extra,'vtype'),
            X.extra.vtype = 'd'; % for now, save everything as double
        end
        save_avw(X.data,X.filename,X.extra.vtype,X.extra.voxsize);
        
    case 'nii_load_nii',
        
        % Write NIFTI with the NIFTI toolbox
        X.filename = horzcat(X.filename,'.nii');
        X.extra.img = X.data;
        save_nii(X.extra,X.filename);
        
    case 'dpxread',
        
        % Write a DPX (DPV or DPF) file
        palm_dpxwrite(X.filename,X.data,X.extra.crd,X.extra.idx);
        
    case 'srfread',
        
        % Write a SRF file
        srfwrite(X.data.vtx,X.data.fac,X.filename);
        
    case 'fs_read_curv',
        
        % Write a FreeSurfer curvature file
        write_curv(X.filename,X.data,X.extra.fnum);
        
    case 'fs_read_surf',
        
        % Write a FreeSurfer surface file
        write_surf(X.filename,X.data.vtx,X.data.fac);
        
    case 'fs_load_mgh',
        
        % Write a FreeSurfer MGH file
        save_mgh(X.data,X.filename,X.extra.M,X.extra.mr_parms);
        
end
