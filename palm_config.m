function varargout = palm_config(args,action)
% Read and write PALM configuration files.
%
% Usage:
% cfg = palm_config(args,'read')
% palm_config(args,'write')
% prefix = palm_config(args,'prefix')
%
% cfg   : Configurations (cell array).
% fname : Text-file with the configurations.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2013
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

switch lower(action)

    case 'read'

        % Read files
        otmp = args{1};
        fid = fopen(otmp,'r');
        args = textscan(fid,'%s','CommentStyle','#');
        fclose(fid);
        if nargout == 1
            varargout = args;
        end

    case 'prefix'

        % Establish a prefix for the _config.txt
        % and _elapsed.csv file. This needs separate
        % treatment because the user may have specified
        % as output an HDF5 or MATLAB .mat file
        optsx = palm_defaults;
        idxa = find(strcmpi(args,'-o'));
        if isempty(idxa)
            otmp = optsx.o;
        else
            otmp = args{idxa+1};
        end
        [filename,~,~] = palm_filespec(otmp);
        fext = palm_tokenize(strcat(filename),'.');
        if any(strcmpi(fext{end},optsx.hdf5)) || ...
               strcmpi(fext{end},'mat')
            [fpth,fnam,~] = fileparts(filename);
            otmp = fullfile(fpth,fnam);
        end
        if ~ strcmp(otmp(end),'_')
            otmp = horzcat(otmp,'_');
        end
        varargout = {otmp};

    case 'write'

        % Establish the prefix for the file
        otmp = palm_config(args,'prefix');
        cfgname = sprintf('%sconfig.txt',otmp);
        [opth,~,~] = fileparts(otmp);
        if ~isempty(opth) && ~exist(opth,'dir')
            mkdir(opth);
        end

        % Version & environment
        [~,ver] = strtok(palm_help('version'),'(');
        if palm_isoctave
            envrun = 'Octave';
        else
            envrun = 'MATLAB';
        end

        % Write files
        fid = fopen(cfgname,'w');
        fprintf(fid,'# Configuration file for PALM.\n');
        fprintf(fid,'# Version %s, running on %s %s.\n',ver,envrun,version);
        fprintf(fid,'# %s\n',datestr(now)); %#ok<TNOW1,DATST>
        fprintf('Running PALM %s on %s %s with the following options:',ver,envrun,version);
        for c = 1:numel(args)
            s2d = str2double(args{c});
            if strcmp(args{c}(1),'-') && (isnan(s2d) || ~isreal(s2d))
                fprintf(    '\n%s',args{c});
                fprintf(fid,'\n%s',args{c});
            else
                if ischar(args{c})
                    fprintf(    ' %s',args{c});
                    fprintf(fid,' %s',args{c});
                else
                    fprintf(    ' %s',num2str(args{c}));
                    fprintf(fid,' %s',num2str(args{c}));
                end
            end
        end
        fprintf(    '\n');
        fprintf(fid,'\n');
        fclose(fid);
        varargout = {otmp};

    otherwise
        error('Incorrect number of input arguments.');
end
