function C = palm_ciftiread(varargin)

% - filename    : Name of the CIFTI file to be loaded.
% - user_prefix : (optional) Prefix (possibly with a path) for temporary files.
% - wb_command  : (optional) Full path to the wb_command executable.

% Handle input arguments
narginchk(1,3);
filename        = varargin{1};
if nargin >= 2,
    user_prefix = varargin{2};
else
    user_prefix = '';
end
if nargin == 3,
    wb_command  = varargin{3};
else
    wb_command  = 'wb_command';
end

% Deal with the location of the temporary files
alphabet    = ['a':'z' 'A':'Z' '0':'9'];
temp_prefix = alphabet(randi(numel(alphabet),[1 6]));

% If no temporary prefix has been supplied, the file will be
% saved to the current directory, with temp_prefix only.
if ~ isempty(user_prefix),
    temp_prefix = strcat(user_prefix,'_',temp_prefix);
end
temp_txt = strcat(temp_prefix,'.txt');
temp_gii = strcat(temp_prefix,'.gii');
temp_nii = strcat(temp_prefix,'.nii');

% Check the CIFTI file contents:
% - If surface, convert to GIFTI.
% - If volume, convert to NIFTI.
% - If both, convert to both.
[~] = system(sprintf('%s -file-information %s | head -n 16 > %s',...
    wb_command,filename,temp_txt));
fileinfo = palm_strcsvread(temp_txt,':');
delete(temp_txt);
idx = strcmpi('Maps to Surface',fileinfo(:,1));
C.MapsToSurface = eval(fileinfo(idx,2));
idx = strcmpi('Maps to Volume', fileinfo(:,1));
C.MapsToVolume  = eval(fileinfo(idx,2));

% For surfaces:
if C.MapsToSurface,

    % Convert to GIFTI and load it
    [~] = system(sprintf('%s -cifti-convert -to-gifti-ext %s %s',...
        wb_command,filename,temp_gii));
    G = gifti(temp_gii);
    C.gifti.extra = G.private;
    C.gifti.data = C.gifti.extra.data{1}.data';
    C.gifti.extra.data{1}.data = [];
                
    % Delete the temporary GIFTI file.
    delete(temp_gii);
    delete(strcat(temp_gii,'.data'));
end

% For volumes:
if C.MapsToVolume,
    
    % Convert to NIFTI and load it
    [~] = system(sprintf('%s -cifti-convert -to-nifti %s %s',...
        wb_command,filename,temp_nii));
    
    % Note this unlinks the mapped file array, and uses more memory.
    % For large files, consider loading directly as NIFTI, not as CIFTI.
    C.nifti.extra = nifti(temp_nii);
    C.nifti.data = double(N.dat);
    C.nifti.extra.dat = [];
    
    % Delete the temporary NIFTI file.
    delete(temp_nii);
end


