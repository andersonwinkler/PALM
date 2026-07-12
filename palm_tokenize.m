function tok = palm_tokenize(str,sep)
% Split a string at the separator "sep" (e.g., '.', ':')
%
% Usage:
% T = palm_tokenize(str,sep)
% 
% Inputs:
% str : String to be separated.
% sep : Separator.
% 
% Output:
% tok : Cell array of separated strings.
% 
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Nov/2015 (first version)
% Mar/2026 (this version)
% http://brainder.org

idx  = find(str == sep);
idxb = [1 idx+1];
idxe = [idx-1 numel(str)];
tok    = cell(numel(idxb),1);
for s = 1:numel(idxb)
    tok{s} = str(idxb(s):idxe(s));
end