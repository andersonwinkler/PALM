function palm(varargin)
% Type 'palm' without arguments for help.

% This is redundant in when running as a function as all files should
% be together, but it helps when running from the shell
addpath(fileparts(mfilename('fullpath')));

% If Octave
if palm_isoctave,
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % If running as a script, take the input arguments
    cmdname = program_invocation_name();
    if ~ strcmpi(cmdname(end-5:end),'octave'),
        varargin = argv();
    end
    
    % Be sure to print to the screen immediately
    page_screen_output(0);
    page_output_immediately(1);
else
    % This line marks the place up to nothing will be printed. It's long as
    % this because if it fails, at least it's not ugly and looks
    % purposeful from the outside :-)
    fprintf('.......................................................................\n');
end

% This is probably redundant but fix a bug in an old Matlab version
nargin = numel(varargin);

% Print usage if no inputs are given
if nargin == 0 || strcmp(varargin{1},'-q'),
    palm_help;
    return;
end

% Now run what matters
palm_backend(varargin{:});
