function palm(varargin)
% Type 'palm -help' for help.

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
