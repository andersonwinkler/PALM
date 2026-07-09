function palm_hdf5_pkg_add(octdir)
% Register hdf5oct autoloads for PALM without pkg install.
%
% palm_hdf5_pkg_add(octdir)
%
% octdir : Directory containing hdf5oct.oct (lib/hdf5/src).
%
% Mirrors the PKG_ADD directives in hdf5oct.cc.

octfile = fullfile(octdir,'hdf5oct.oct');
if exist(octfile,'file') == 0
    error('hdf5oct.oct not found in %s',octdir);
end

autoload('__h5read__',    octfile);
autoload('__h5readatt__', octfile);
autoload('__h5write__',   octfile);
autoload('__h5writeatt__',octfile);
autoload('__h5create__',  octfile);
autoload('__h5delete__',  octfile);
autoload('h5info',        octfile);
