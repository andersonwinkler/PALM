function b = palm_d2b(d,n)
% Converts decimal integers to binary.
% 
% b = palm_d2b(d,n)
%
% d : Decimal value(s).
% n : Word size (number of columns of b).
% b : Binary value(s). The least significant
%     is the rightmost column.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2013
% http://brainder.org

[~,e] = log2(max(d));
b = rem(floor(d(:)*pow2(1-max(n,e):0)),2);