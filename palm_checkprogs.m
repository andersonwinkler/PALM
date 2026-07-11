function ext = palm_checkprogs
% Test whether some external programs or toolboxes
% are available.
%
% ext = palm_checkprogs
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

persistent palm_extern;
if isempty(palm_extern)

    % Check the path of PALM and add the paths for file I/O and other pieces
    palm_extern.palmpath = fileparts(mfilename('fullpath'));
    fprintf('PALM is located at %s\n',palm_extern.palmpath);
    addpath(fullfile(palm_extern.palmpath,'lib'));
    addpath(fullfile(palm_extern.palmpath,'lib','extras'));
    addpath(fullfile(palm_extern.palmpath,'lib','freesurfer'));
    addpath(fullfile(palm_extern.palmpath,'lib','cifti-matlab'));
    addpath(fullfile(palm_extern.palmpath,'lib','arrow3'));
    addpath(fullfile(palm_extern.palmpath,'lib','colourmaps'));

    % External programs - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Check FSL
    palm_extern.fsl = false;
    fsldir = getenv('FSLDIR');
    if ~isempty(fsldir)
        palm_extern.fsl = true;
        addpath(fullfile(fsldir,'etc','matlab'));
        fprintf('Found FSL in %s\n',fsldir);
    end

    % Check FreeSurfer
    palm_extern.fs  = false;
    fshome = getenv('FREESURFER_HOME');
    if ~isempty(fshome)
        palm_extern.fs = true;
        addpath(fullfile(fshome,'matlab'));
        fprintf('Found FreeSurfer in %s\n',fshome);
    end

    % Check SPM
    palm_extern.spm = false;
    try %#ok
        spm_check_installation('basic');
        palm_extern.spm = true;
        spmpath = fileparts(which('spm'));
        fprintf('Found SPM in %s\n',spmpath);
    end

    % Check HCP Workbench
    palm_extern.wb_command = false;
    [status,wb_command] = system('which wb_command');
    if status == 0
        palm_extern.wb_command = true;
        fprintf('Found HCP Workbench executable in %s',wb_command);
    end

    % Check DuckDB (to read Parquet files in Octave)
    palm_extern.duckdb = false;
    [status,duckdb] = system('which duckdb');
    if status == 0
        palm_extern.duckdb = true;
        fprintf('Found DuckDB executable in %s',duckdb);
    end

    % Octave packages - - - - - - - - - - - - - - - - - - - - - - - - - - -
    palm_extern.octave_image      = false;
    palm_extern.octave_specfun    = false;
    palm_extern.octave_statistics = false;
    if palm_isoctave
        pkg_installed = pkg('list');
        pkg_names     = cellfun(@(p)p.name,pkg_installed,'UniformOutput',false);
        palm_extern.octave_image      = any(strcmp(pkg_names,'image'));
        palm_extern.octave_specfun    = any(strcmp(pkg_names,'specfun'));
        palm_extern.octave_statistics = any(strcmp(pkg_names,'statistics'));
        if palm_extern.octave_image
            fprintf('Octave package "image" is available.\n');
        else
            fprintf('Octave package "image" is not available.\n');
        end
        if palm_extern.octave_specfun
            fprintf('Octave package "specfun" is available.\n');
        else
            fprintf('Octave package "specfun" is not available.\n');
        end
        if palm_extern.octave_statistics
            fprintf('Octave package "statistics" is available.\n');
        else
            fprintf('Octave package "statistics" is not available.\n');
        end

        % Internalized hdf5oct package (will shadow an existing hdf5oct)
        hdf5dir = fullfile(palm_extern.palmpath,'lib','hdf5','src');
        addpath(hdf5dir);
        octfile = fullfile(hdf5dir,'hdf5oct.oct');
        if exist(octfile,'file') == 0
            fprintf('Internal HDF5 library is not compiled for your platform. If you plan to use HDF5 files, consult the documentation.\n');
        else
            fprintf('Internal HDF5 library is available.\n');
            autoload('__h5read__',    octfile);
            autoload('__h5readatt__', octfile);
            autoload('__h5write__',   octfile);
            autoload('__h5writeatt__',octfile);
            autoload('__h5create__',  octfile);
            autoload('__h5delete__',  octfile);
            autoload('h5info',        octfile);
        end
    end

    % MATLAB toolboxes  - - - - - - - - - - - - - - - - - - - - - - - - - -
    palm_extern.matlab_ipt      = false;
    palm_extern.matlab_symbolic = false;
    if ~ palm_isoctave
        % Image Processing Toolbox
        if license('test','Image_Toolbox')
            palm_extern.matlab_ipt = true;
            fprintf('Image Processing Toolbox is available.\n');
        elseif  (exist('niftiread', 'builtin') == 5 || exist('niftiread', 'file') == 2) && ...
                (exist('niftiwrite','builtin') == 5 || exist('niftiwrite','file') == 2) && ...
                (exist('niftiinfo', 'builtin') == 5 || exist('niftiinfo', 'file') == 2)
            palm_extern.matlab_ipt = true;
            fprintf('Internal NIFTI read/write functions are available.');
        end
        % Symbolic Math Toolbox
        if license('test','Symbolic_Toolbox')
            palm_extern.matlab_symbolic = true;
            fprintf('Symbolic Math Toolbox is available.\n');
        end
    end
end
ext = palm_extern;
