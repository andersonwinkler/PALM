function y = palm_isoctave
% Test whether the function is running in OCTAVE or not.
% Retuns 1 (true) if yes, 0 (false) otherwise.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2012
% http://brainder.org

persistent palm_isoct;
if isempty(palm_isoct),
    palm_isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end;
y = palm_isoct;
