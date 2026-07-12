% RobustICA algorithm for independent component analysis (Release 3 - November 21, 2014)
% --------------------------------------------------------------------------------------
%
% RobustICA is based on the normalized kurtosis contrast function, which is optimized by a
% computationally efficient iterative technique. This technique computes algebraically
% the step size (adaption coefficient) globally optimizing the contrast in the search direction
% at each iteration. Any independent component with non-zero kurtosis can be extracted in this
% manner.
%
% The present implementation performs the deflationary separation of statistically independent
% sources under the instantaneous linear mixture model. Full separation is achieved if at most
% one source has zero kurtosis.
%
%
% ADVANTAGES OF RobustICA:
% 
% - Real- and complex-valued signals are treated by exactly the same algorithm. Both type of
%   source signals can be present simultaneously in a given mixture. Complex sources need not be
%   circular. The mixing matrix coefficients may be real or complex, regardless of the source type.
%
% - Prewhitening is not required, so that the performance limitations it imposes can be avoided.
%   In particular, the absence of prewhitening improves asymptotic performance. In that case, 
%   sequential extraction (deflation) can be performed via linear regression.
%
% - The algorithm can target sub-Gaussian or super-Gaussian sources in the order defined by a
%   by a kurtosis-sign vector provided by the user. Full separation, as well as the consequent
%   increase in complexity and estimation error, can be spared if the Gaussianity character of
%   the source(s) of interest is known in advance.
%
% - The optimal step-size technique provides some robustness to the presence of saddle points and
%   spurious local extrema in the contrast function, which tend to appear when processing short
%   data sizes.
%
% - The method shows a very high convergence speed measured in terms of source extraction quality
%   versus number of operations. In the real-valued two-signal case, the algorithm theoretically 
%   converges in a single iteration, even without prewhitening.
%
% More details about the RobustICA algorithm and its comparative performance analysis can be found
% in references [1]-[3] below.
%
%
% WHAT'S NEW ON RELEASE 3:
% 
% - If required, the user can force the separating matrix to be real-valued. This is useful, e.g., 
%   when the mixing matrix takes real values but the separation is carried out in a complex domain
%   (e.g., after Fourier transform).
%
% - The calling syntax has been simplified by using a cell-array input argument.
%
% - The code is now licensed under a Creative Commons Non-Commercial License
%   (see 'Conditions of use' for details).
%
%
% WHAT WAS NEW ON RELEASE 2 (Feb. 16, 2010):
% 
% - The optimal step-size computation has been refined
%   ('kurt_gradient_optstep.m' function).
%  
% - Circular and noncircular complex sources are now included in the demonstration
%   ('robustica_demo.m' function).
% 
% - For the sake of reproducible research, the code used to generate some of the results
%   reported in reference [1] has been added to the package
%   ('robustica_tnn_fig3_sim.m' and related functions).
%
%
% M-FILES:
%
% The package is composed of the following M-files:
%                                          
%  - <a href = "matlab:doc robustica">robustica.m</a>:             implements the algorithm itself;
%                             the algorithm calls the three functions below.
%
%       o <a href = "matlab:doc kurt_gradient_optstep">kurt_gradient_optstep.m</a>: computes the optimal step-size of the normalized kurtosis
%                                  contrast using the gradient vector as search direction. 
%
%       o <a href = "matlab:doc deflation_regression">deflation_regression.m</a>:  performs deflation via linear regression.
%
%       o <a href = "matlab:doc get_arguments">get_arguments.m</a>:        gets algorithm inputs parameters from RobustICA method call.
%
%  - <a href = "matlab:doc robustica_demo">robustica_demo.m</a>:        a simple demonstration illustrating the performance of RobustICA
%                             on synthetic mixtures.
% 
%  - <a href = "matlab:doc robustica_tnn_fig3_sim">robustica_tnn_fig3_sim.m</a>:  reproduces the simulation of Fig. 3 (Section IV.B) of reference [1].
%
%       o <a href = "matlab:doc compute_smse">compute_smse.m</a>:   computes performance parameter from estimated sources.
%
%       o <a href =  "matlab:doc completion_bar">completion_bar.m</a>: shows the elapsed and remaining time of a running simulation.
%   
%       o <a href =  "matlab:doc fastica">fastica.m</a>:        FastICA algorithm with cubic nonlinearity.
%   
%
% CONDITIONS OF USE:
% 
%    This is open-source code. Please, feel free to edit the M-files and modify them as you wish.
%    It is only requested that original authorship be acknowledged and modifications be clearly
%    indicated in the code (history section and elsewhere as appropriate).
%     
%    If you use this package to generate results in your publications, please cite at least reference [1].
%
%    This code is licensed under a <a href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
%
%
% COMPATIBILITY:
% 
%    The package has been tested on Matlab 7.14.
%
%
% UPDATES:
% 
%    To download the latest version of the package, please visit the RobustICA webpage:
%
%       <a href = "http://www.i3s.unice.fr/~zarzoso/robustica.html">http://www.i3s.unice.fr/~zarzoso/robustica.html</a>
%
%
% FEEDBACK:
% 
%    Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% REFERENCES:
%       
% [1] V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust Independent Component Analysis by Iterative Maximization</a>
%     <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the Kurtosis Contrast with Algebraic Optimal Step Size"</a>, 
%     IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
%
% [2] V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/ica07.pdf">"Comparative Speed Analysis of FastICA"</a>, 
%     in: Proceedings ICA-2007, 7th International Conference on Independent Component Analysis
%         and Signal Separation, London, UK, September 9-12, 2007, pp. 293-300.
% 
% [3] V. Zarzoso, P. Comon and M. Kallel, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco06.pdf">"How Fast is FastICA?"</a>, 
%     in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference, 
%         Florence, Italy, September 4-8, 2006. 
%