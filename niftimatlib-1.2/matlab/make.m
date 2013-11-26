%
% Makefile equivalent for use in MatLab
% requires availability of MEX
%
% calling 'make' from the matlab prompt:
% 1. compiles the c-sources in @file_array/private/src
% 2. moves the compiled files to @file_array/private
%
% Alle Meije Wink 27/03/2012
% a.wink@vumc.nl
%

%
% niftilib $Id: make.m,v 1.1 2012/03/30 15:25:41 fissell Exp $
%

% check for existence of MEX compiler
if (exist('mex') ~=2)
    error('MEX not installed on your system. Exiting.');
end

% go to the directory where the c sources are
cd(['@file_array' filesep 'private' filesep 'src'])

% compile the 2 sources
if(mex('file2mat.c') | mex('mat2file.c'))
    error('something went wrong during compilation');
else
    % move them to the parent directory
    if(~movefile('*.mex*','..'))
        error('compiled files could not be moved')
    end
end

% go bake to the directory of make.m
cd(['..' filesep '..' filesep '..']);

% report the good news
disp ('files successfully compiled and moved')


