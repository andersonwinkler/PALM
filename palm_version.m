function varargout = palm_version


% Read the file with the version
fid = fopen(fullfile(fileparts(mfilename('fullpath')),'palm_version.txt'),'r');
vstr = textscan(fid,'%s');
fclose(fid);

% Assemble back as a string
vstr = sprintf('%s ',vstr{1}{:});
vstr = sprintf('%s\n',vstr(1:end-1));

% Print in the screen if no other output
if nargout == 0,
    fprintf(vstr);
else
    varargout{1} = vstr;
end