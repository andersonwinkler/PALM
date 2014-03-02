function lfac = palm_factorial(N)
% Computes the log(factorial(0:N)), so dealing with
% precision issues.
%
% lfac = palm_factorial(N)
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2012
% http://brainder.org

persistent lf;
if isempty(lf) || length(lf) < N+1,
    lf = zeros(N+1,1);
    for n = 1:N,
        lf(n+1) = log(n) + lf(n);
    end
end
lfac = lf;