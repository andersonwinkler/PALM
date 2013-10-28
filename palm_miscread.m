function X = palm_miscread(filename)
% Read various scalar data formats based on the file extension.
%
% X = miscread(filename);
%
% X is a struct that contains the fields:
% X.filename : Contains the name of the file.
% X.readwith : This tells which program or function was used
%              to read the data. This is useful when saving the
%              data back, to use a compatible function.
% X.data     : Array with the actual data. The size can vary
%              according to what was read.
% X.extra    : Contain extra information, depending on the kind
%              of data that was read and the function or
%              program used for reading.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Check if the file actually exists before doing anything else
if ~ exist(filename,'file')
    error('File not found: %s',filename);
end

% Store the filename, in case there is a need to overwrite this later
X.filename = filename;

% Take the file extension and try to load accordingly
[~,~,fext] = fileparts(X.filename);
switch lower(fext),
    
    case '.txt',
        
        % Read a generic text file
        X.readwith = 'textscan';
        fid = fopen(X.filename);
        X.data = textscan(fid,'%s');
        X.data = X.data{1};
        fclose(fid);
        
    case '.csv',
        
        % Read a CSV file. It has to contain numeric values only.
        % The command 'csvwrite' is a frontend to 'dlmwrite', which calls
        % 'textscan', which on its turn has a limitation of 100k columns.
        % Using 'load' bypass this issue.
        X.readwith = 'load';
        X.data = load(X.filename,'-ascii');
        
    case {'.mat','.con','.fts','.grp'},
        
        % Read an FSL "VEST" file.
        X.readwith = 'vestread';
        [X.data,X.extra.PPH] = palm_vestread(X.filename);
        
    case '.gz',
        
        % Read a gzipped NIFTI file.
        extern = palm_checkprogs;
        if extern.fs,       % Read with FreeSurfer
            X.readwith = 'fs_load_nifti';
            X.extra.hdr = load_nifti(X.filename);
            X.data = X.extra.hdr.vol;
            X.extra.hdr.vol = [];
        elseif extern.fsl,  % Read with FSL
            X.readwith = 'fsl_read_avw';
            [X.data,X.extra.dims,X.extra.voxsize, ...
                X.extra.bpp,X.extra.endian] = read_avw(X.filename);
        else
            error([
                'Neither FSL or FreeSurfer were found.\n' ...
                'To use this data you must do one of:\n' ...
                '- Make sure FSL is correctly installed and configured,\n' ...
                '  and that your ''FSLDIR'' environmental variable is\n' ...
                '  properly set;\n' ...
                '- Make sure FreeSurfer is correctly installed and configured,\n' ...
                '  and that your ''FREESURFER_HOME'' environmental variable is\n' ...
                '  properly set;\n' ...
                '- If you do not have FSL or FS, but you have SPM, then\n' ...
                '  uncompress the .gz file and try again.\n' ...
                'File: %s'],X.filename);
        end
        
    case {'.nii','.hdr','.img'},
        
        % Read a non-gzipped NIFTI file.
        extern = palm_checkprogs;
        if extern.spm,      % Read with SPM
            X.readwith = 'spm_spm_vol';
            X.extra = spm_vol(X.filename);
            X.data = spm_read_vols(X.extra);
        elseif extern.fs,   % Read with FreeSurfer
            X.readwith = 'fs_load_nifti';
            X.extra.hdr = load_nifti(X.filename);
            X.data = X.extra.hdr.vol;
            X.extra.hdr.vol = [];
        elseif extern.fsl,  % Read with FSL
            X.readwith = 'fsl_read_avw';
            [X.data,X.extra.dims,X.extra.voxsize, ...
                X.extra.bpp,X.extra.endian] = read_avw(X.filename);
        elseif extern.nii,  % Read with the NIFTI toolbox
            X.readwith = 'nii_load_nii';
            X.extra = load_nii(X.filename);
            X.data = X.extra.img;
            X.extra.img = [];
        else
            error([
                'No FSL, FreeSurfer, SPM or the NIFTI toolbox were\n' ...
                'found. To use this data you must have one of these\n' ...
                'installed and correctly configured.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case {'.dpv','.dpf','.dpx'},
        
        % Read a DPV/DPF file, in ASCII
        X.readwith = 'dpxread';
        [X.data,X.extra.crd,X.extra.idx] = dpxread(X.filename);
        
    case '.srf',
        
        % Read a SRF file, in ASCII
        X.readwith = 'srfread';
        [X.data.vtx,X.data.fac] = srfread(X.filename);
        
    case {'.area','.avg_curv','.crv','.curv', ...
            '.h','.k','jacobian_white','.mid', ...
            '.pial','.sulc','.thickness','.volume'},
        
        % Read a FreeSurfer curvature file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_read_curv';
            [X.data,X.extra.fnum] = read_curv(X.filename);
        else
            error([
                'FreeSurfer was not found. To use this data, makes sure\n' ...
                'that FreeSurfer is correctly installed and configured and,\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case {'inflated','inflated.nofix','orig', ...
            'orig.nofix','pial','qsphere.nofix','smoothwm', ...
            'smoothwm.nofix','sphere','sphere.reg','white'},
        
        % Read a FreeSurfer surface file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_read_surf';
            [X.data.vtx,X.data.fac] = read_surf(X.filename);
        else
            error([
                'FreeSurfer was not found. To use this data, makes sure\n' ...
                'that FreeSurfer is correctly installed and configured and,\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case {'.mgh','.mgz'},
        
        % Read a FreeSurfer MGH/MGZ file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_load_mgh';
            [X.data,X.extra.M,X.extra.mr_parms,X.extra.volsz] = load_mgh(X.filename);
        else
            error([
                'FreeSurfer was not found. To use this data, makes sure\n' ...
                'that FreeSurfer is correctly installed and configured and,\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    otherwise
        error('File extension %s not known. Data cannot be loaded\n',fext);
end
