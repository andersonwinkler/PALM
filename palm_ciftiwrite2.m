function palm_ciftiwrite(fileprefix,C)

% For surfaces:
if C.MapsToSurface,
    
    % Convert the data to a GIFTI object
    G = gifti([]);
    F = fieldnames(C.gifti);
    for f = 1:numel(F),
        G.private.(F{f}) = C.gifti.extra.(F{f});
    end
    G.cdata = C.gifti.data;
    
    % Edit the XML-tree for the new file size
    xt = xmltree(G.private.data{1}.metadata.value);
    uid = find(xt,'/CIFTI/Matrix/MatrixIndicesMap');
    for u = 1:numel(uid),
        attr = attributes(xt,'get',uid(u));
        for a = 1:numel(attr),
            if strcmp(attr{a}.key,'NumberOfSeriesPoints'),
                xt = attributes(xt,'set',uid(u),a,'NumberOfSeriesPoints',num2str(size(G.cdata,2)));
                break;
            end
        end
    end
    
    % Save the GIFTI. Note that "save" is overloaded.
    % This will create a .gii and a .gii.dat.
    save(G,strcat(fileprefix,'.gii'),'ExternalFileBinary');
    
    % ---
    % Use the Workbench to convert the GIFTI to NIFTI.
    [~] = system(sprintf('wb_command -cifti-convert -from-gifti-ext %s %s',...
        strcat(fileprefix,'.gii'),strcat(fileprefix,'.dscalar.nii')));
    
    % Delete the temporary GIFTI file.
    delete(strcat(fileprefix,'.gii'));
    delete(strcat(fileprefix,'.dat'));
end

% For volumes:
if C.MapsToVolume,
    
    % Write as a NIFTI object.
    dat = file_array(            ...
        strcat(filename,'.nii'), ...
        size(C.nifti.data),      ...
        'FLOAT32-LE',            ...
        ceil(348/8)*8);
    nii      = nifti;
    nii.dat  = dat;
    nii.mat  = C.nifti.extra.mat;
    if isfield(C.nifti.extra,'mat0')
        nii.mat0 = C.nifti.extra.mat0;
    else
        nii.mat0 = C.nifti.extra.mat;
    end
    create(nii);
    nii.dat(:,:,:) = C.nifti.data(:,:,:);
    clear nii;
    
    % Use the Workbench to convert the GIFTI to NIFTI.
    [~] = system(sprintf('wb_command -cifti-convert -from-nifti %s %s',...
        strcat(fileprefix,'.nii'),strcat(fileprefix,'.dscalar.nii')));
    
    % Delete the temporary NIFTI file.
    delete(strcat(fileprefix,'.nii'));
end

