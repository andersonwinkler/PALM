function [g, mu_opt, norm_g] = kurt_gradient_optstep(w, X, s, P, wreal)

% Computes optimal step size in the gradient-based optimization of the normalized kurtosis contrast
% (single iteration).
%
% Data-based version.
%
% See references below for details.
%
%
% SYNTAX: [g, mu_opt, norm_g] = kurt_gradient_optstep(w, X, s, P, wreal); 
%
%
% OUTPUTS:
%         g      : search direction (normalized gradient vector)
%
%         mu_opt : optimal step size globally optimizing the normalized kurtosis contrast function
%                  along direction g from f
%
%         norm_g : non-normalized gradient vector norm. 
%
%
% INPUTS:
%         w      : current extracting vector coefficients
%
%         X      : sensor-output data matrix (one signal per row, one sample per column)
%
%         s      : source kurtosis sign; if zero, the maximum absolute value of the contrast is sought
%
%         P      : projection matrix (used in deflationary orthogonalization; identity matrix otherwise)
%
%         wreal  : if different from zero, keep extracting vector real valued by retaining only the
%                  real part of the gradient (useful, for instance, in separating real-valued mixtures
%                  in the frequency domain, as in the RobustICA-f algorithm).
%
%
% REFERENCES:
%
% - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust independent component analysis by iterative maximization</a>
%   <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the kurtosis contrast with algebraic optimal step size"</a>, 
%   IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
%
% - V. Zarzoso, P. Comon and M. Kallel,  <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco06.pdf">"How fast is FastICA?"</a>, 
%   in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference, 
%       Florence, Italy, September 4-8, 2006. 
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% HISTORY:
% 
%    <Please add modification date here>: - <please add modification details here>
%    
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
%    2014/06/25: - added 'wreal' input parameter to allow the separation of real-valued mixtures
%                  in complex (e.g., frequency) domain
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
%    2010/02/09: - include gradient norm as output parameter (for use as an additional termination criterion)
%
%    2009/03/04: - removed test for constant contrast; sometimes the algorithm stopped too early, 
%                  because the contrast was not actually constant, leading to suboptimal extraction results
%                - if best candidate root is complex valued, its real part can be retained as optimal
%                  step size, but contrast is not guaranteed to increase monotonically in that case; 
%                  to avoid this problem, only the real parts of the roots are considered
%
%    2009/03/02: - simplified expressions of gradient and optimal step size, as in TNN paper
% 
%    2008/04/01: - problem encountered on 2008/03/25 using orthogonalization: due to nearly-zero
%                  gradient appearing when just one source is left, since the contrast
%                  function then becomes constant; normalization after projection in such a case
%                  destroys co-linearity between gradient and extracting vector
%                  (checking for zero gradient should probably be used as additional termination test
%                   in the next version of the algorithm; see modification on 2010/02/09)
% 
% -- 2008/03/31: Version 1 release ----------------------------------------------------------------
%
%    2008/03/26: - added this help
%
%    2008/03/25: - projecting the gradient after normalization seems to improve conditioning
%                  and accelerate convergence in the extraction of the last sources
%
%    2008/03/24: - created by Vicente Zarzoso (University of Nice - Sophia Antipolis, France).


verbose = 0;

[L, T] = size(X);

mu_opt = 0;     % default optimal step-size value
norm_g = 0;     % initialize gradient norm


%%% Compute search direction (gradient vector)

% compute necessary interim values

y = w'*X;

ya2 = y.*conj(y);
y2 = y.*y;
ya4 = ya2.*ya2;

Eya2 = mean(ya2);
Ey2 = mean(y2);
Eya4 = mean(ya4);

if abs(Eya2) < eps      % check for zero denominator
    if verbose, disp('>>> OPT STEP SIZE: zero power'), end
    g = zeros(L, 1);
    norm_g = 0;
else
    
% compute gradient if contrast denominator is not null
    Eycx = X*y'/T;
    Eyx = X*y.'/T;
    Ey3x = X*(ya2.*y)'/T;

    % contrast numerator and denominator at current point
    p1 = Eya4 - abs(Ey2)^2;
    p2 = Eya2;

    g = 4*( (Ey3x - Eyx*Ey2')*p2 - p1*Eycx )/p2^3;
 
    g = P*g;            % project if required (normalize later)
 
    norm_g = norm(g);
    
 if norm_g < eps
    if verbose,  disp('>>> OPT STEP SIZE: zero gradient'); end
 else       
 
 if wreal, g = real(g); end;    % keep only real part if real-valued extracting vectors are required    
     
 g = g/norm_g;          % normalize the gradient -> parameter of interest: direction 
                        % improves conditioning of opt step-size polynomial  
  
                         
%%% Compute optimal step size

gg = g'*X;


% calculate interim values for contrast rational function

ya2 = y.*conj(y);
ga2 = gg.*conj(gg);
reygc = real(y.*conj(gg));

g2 = gg.*gg; yg = y.*gg;

Eya2reygc = mean(ya2.*reygc);
Ereygc2 = mean(reygc.^2);
Ega2reygc = mean(ga2.*reygc);
Ega4 = mean(ga2.^2);
Eya2ga2 = mean(ya2.*ga2);

Ega2 = mean(ga2);
Ereygc = mean(reygc);

Eg2 = mean(g2); Eyg = mean(yg);

h0 = Eya4 - abs(Ey2)^2;
h1 = 4*Eya2reygc - 4*real(Ey2*Eyg');
h2 = 4*Ereygc2 + 2*Eya2ga2 - 4*abs(Eyg)^2 - 2*real(Ey2*Eg2');
h3 = 4*Ega2reygc - 4*real(Eg2*Eyg');
h4 = Ega4 - abs(Eg2)^2;

P = [h4, h3, h2, h1, h0];

i0 = Eya2; i1 = 2*Ereygc; i2 = Ega2;

Q = [i2, i1, i0];

% normalized kurtosis contrast = P/Q^2 - 2

a0 = -2*h0*i1 + h1*i0;
a1 = -4*h0*i2 - h1*i1 + 2*h2*i0;
a2 = -3*h1*i2 + 3*h3*i0;
a3 = -2*h2*i2 + h3*i1 + 4*h4*i0;
a4 = -h3*i2 + 2*h4*i1;

p = [a4, a3, a2, a1, a0];

% normalized kurtosis contrast derivative = p/Q^3


% %%% ALTERNATIVE METHOD to compute optimal step-size polynomial oefficients
% 
% % obtain contrast-function polynomials
% 
% p11 = [Ega4, 4*Ega2reygc, 4*Ereygc2+2*Eya2ga2, 4*Eya2reygc, Eya4];
% p13 = [Eg2, 2*Eyg, Ey2];
% P = p11 - conv(p13, conj(p13));     % numerator
% Q = [Ega2, 2*Ereygc, Eya2];         % square-root of denominator
% 
% % compute derivatives
% Pd = [4, 3, 2, 1].*P(1:4);
% Qd = [2, 1].*Q(1:2);
% 
% % contrast derivative numerator
% p = conv(Pd, Q) - 2*conv(Qd, P);   


rr = real(roots(p));        % keep real parts only

Pval = polyval(P, rr);
Q2val = polyval(Q, rr).^2;

nonzero_Q2val = find(Q2val > eps);     % check roots not shared by denominator
                            % NOTE: in theory, the denominator can never
                            % cancel out if the gradient is used as a search direction, due to the orthogonality
                            % between the extracting vector and the corresponding gradient 
                            % (only exception: if it is the last source to be extracted; 
                            %  but this scenario is detected by the gradient norm)
                            
if isempty(nonzero_Q2val)  
    if verbose
       disp('>>> OPT STEP SIZE: all roots shared by denominator')
       Pval'
       Q2val'
       p, P, Q
       Q2 = conv(Q, Q); P_Q2 = P./Q2
       pause
   end
else
   Pval = Pval(nonzero_Q2val);
   Q2val = Q2val(nonzero_Q2val);
   rr = rr(nonzero_Q2val);
    
   Jkm_val = Pval./Q2val - 2;       % normalized kurtosis

   if s
    Jkm_val = real(s*Jkm_val);      % maximize or minimize kurtosis value, depending on kurtosis sign
   else
    Jkm_val = abs(Jkm_val);         % maximize absolute kurtosis value, if no sign is given
   end

   [Jmax, im] = max(Jkm_val);   
   mu_opt = rr(im);                 % optimal step size
   
end % if isempty(nonzero_Q2val)

end % if norm(g) < eps

end % if abs(Eya2) < eps

