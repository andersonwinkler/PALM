function varargout = palm_configrw(varargin)
% Read and write PALM configuration files.
% 
% Usage:
% cfg = palm_configrw(fname)
% palm_configrw(cfg,fname)
% 
% cfg   : Configurations (cell array).
% fname : Text-file with the configurations.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2013
% http://brainder.org

if nargin == 1,
    
    % Read files
    plmfile = varargin{1};
    fid = fopen(plmfile,'r');
    cfg = textscan(fid,'%s','CommentStyle','#');
    fclose(fid);
    if nargout == 1,
        varargout = cfg;
    end
    
elseif nargin == 2,
    
    % Write files
    cfg = varargin{1};
    plmfile = varargin{2};
    fid = fopen(plmfile,'w');
    fprintf(fid,'# Configuration file for PALM.\n');
    fprintf(fid,'# %s\n',datestr(now));
    fprintf('Running PALM with the following user supplied options:');
    for c = 1:numel(cfg),
        if strcmp(cfg{c}(1),'-'),
            fprintf('\n%s',cfg{c});
            fprintf(fid,'\n%s',cfg{c});
        else
            if ischar(cfg{c}),
                fprintf(' %s',cfg{c});
                fprintf(fid,' %s',cfg{c});
            else
                fprintf(' %s',num2str(cfg{c}));
                fprintf(fid,' %s',num2str(cfg{c}));
            end
        end
    end
    fclose(fid);
    fprintf('\n');
else
    error('Incorrect number of input arguments.');
end
