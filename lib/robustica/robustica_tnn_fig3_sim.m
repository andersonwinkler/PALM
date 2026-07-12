% RobustICA's comparative quality-cost performance for different mixture sizes
% in the real-valued case, fixed sample size, with and without prewhitening.
%
% Reproduces the simulation of Fig. 3 (Section IV.B) reported in
%
% V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust Independent Component Analysis by Iterative Maximization</a>
% <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the Kurtosis Contrast with Algebraic Optimal Step Size"</a>, 
% IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
%
%
% HISTORY:
%
%    <Please add modification date here>: - <please add modification details here>
%
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
% 2014/11/21: - complies with new calling structure of 'robustica.m'.
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
% 2010/02/10: - created and added to RobustICA package by Vicente Zarzoso
%               (I3S Laboratory, University of Nice Sophia Antipolis, CNRS, France).


clear all; clc

filename = 'robustica_tnn_fig3_sim_results';    % file where simulation results will be stored

disp(' ');
disp('=====================================================================================================');
disp(' ');


%%%%% Simulation parameters

T = 150;                % sample size
T_label = num2str(T);
read_me = ['Separation quality-cost trade-off, real case, noiseless orthogonal mixtures, T = ', T_label, ' samples'];
disp(read_me);

nruns = 5               % number of Monte Carlo runs ('nruns = 1000' in Fig. 3 of [1])

K = [5, 10, 20]         % number of sources used in the simulation
nK = length(K);

verbose = false         % verbose operation

npoints = 11;           % number of points per curve in quality-cost plot

ITER = round(logspace(log10(1), log10(200), npoints))    % number of iterations -> determines complexity axis
all_iter = sum(ITER);

tol = -1                % negative termination-test threshold -> fix overall complexity

methods = {'fastica', 'fastica', 'robustica', 'robustica'}
nmet = length(methods);
par_prewhi = [false, true, false, true]           % prewhitening
par_deftype = {'orthogonalization', 'orthogonalization', 'orthogonalization', 'orthogonalization'}  % deflation: orthogonalization
par_dimred = [false, false, false, false]         % do not perform dimensionality reduction

flops_iter = zeros(nK, nmet);       % flops per iteration per source for each method and mixture size
flops_prewhi = flops_iter;          % cost of prewhitening per source for each method and mixture size

SMSE = zeros(nmet, npoints, nK);    % average performance index for each method, at each independent parameter value
                                    % and mixture size
                                    
figname = ['Quality-cost trade-off, T = ', T_label];
linplotK = {'-', '--', ':'};        % line types for different mixture sizes
linplotmet = {'x', '+', '^', 'o'};  % line markers for different methods 
fig_axis = [10^1, 10^5, -30, 0];
fig_xtick = 10.^[1:5];
fig_xlabel = 'complexity per source per sample (flops)'; 
fig_ylabel = 'SMSE (dB)';

disp(' '); disp('Please, press any key to start the simulations...'); disp(' ');
pause

%%%%%% Perform a simulation for each mixture size (number of sources and observations)

for kpos = 1:nK     % loop over number of sources
    
k = K(kpos);
k_label = num2str(k);

disp(' ');
disp('----------------------------------------------------------------------------');
disp(' ');
read_me_k = [read_me, ', K = ', k_label, ' sources'];
disp(read_me_k);


%%% generate sources and mixtures

S_total = (rand(k, T*nruns) > 0.5);     % pseudo-random binary sequences
S_total = S_total - mean(S_total')'*ones(1, T*nruns);               % remove mean
S_total = diag(1./sqrt(diag(S_total*S_total')/(T*nruns)))*S_total;  % unit-power normalization

H_total = zeros(k, k*nruns);
for r = 1:nruns
    Hr = randn(k);                      % random Gaussian mixing matrix
    [Uh, Sh, Vh] = svd(Hr);
    Hr = Uh*Vh';                        % unitary mixture
    H_total(:, (r-1)*k+1:r*k) = Hr;     % store
end % for r

% cost per iteration for current mixture size (real-valued mixtures)
ffica_iter = 2*(k +  1)*T;   % flops per iteration for FastICA with cubic nonlinearity
frica_iter = (5*k+12)*T;     % flops per iteration for RobustICA

flops_iter(kpos, :) = [ffica_iter, ffica_iter, frica_iter, frica_iter];     % complexity per source per iteration
flops_prewhi(kpos, :) = (2*k^2*T)/k*par_prewhi;                             % cost of prewhitening per source, if employed

Wini = eye(k);       % canonical basis initialization; same for all methods compared

[t0_bar, old_done, last_text] = completion_bar(0, nruns*all_iter);   % initialize completion bar


%%% loop over independent parameter values

for p = 1:npoints       % loop over independent parameter values (number of iterations)

    max_it = ITER(p);       % all methods will perform this number of iterations per source extraction

    smse = zeros(nruns, nmet);  % extraction results for each method-parameter (+MMSE) and each run


    % Monte Carlo runs

    for r = 1:nruns         % loop over signal realizations
  
        % obtain mixture realization 
        S = S_total(:, (r-1)*T+1:r*T);          % retrieve current source realization
        S = S - mean(S')'*ones(1, T);           % remove (residual) mean from source realization
        S = diag(1./sqrt(diag(S*S')/T))*S;      % unit-power sources

        H = H_total(:, (r-1)*k+1:r*k);          % retrieve current mixing matrix

        X = H*S;    % generate current mixture realization
          

        for m = 1:nmet % loop over methods

            % obtain source extraction for given method and operation parameters
            arguments = {'tol', tol,   'maxiter', max_it, 'prewhi', par_prewhi(m), 'deftype', par_deftype{m}, 'dimred', par_dimred(m), ...
                         'Wini', Wini, 'verbose', verbose};
                      
%             [Se, He, it] = feval(methods{m}, X, [], tol, max_it, par_prewhi(m), par_deftype(m), par_dimred(m), Wini, verbose);
            [Se, He, it] = feval(methods{m}, X, arguments);


            % measure performance
            smse(r, m) = compute_smse(S, Se);   % compute source mean square error for current signal realization and method,
                                                % after appropriate ordering, scaling and phase correction

        end % for m (loop over methods)


        % show elapsed and remaning time for completion after each MC run
        [t0_bar, old_done, last_text] = completion_bar(nruns*sum(ITER(1:p-1))+r*ITER(p), nruns*all_iter, t0_bar, old_done, last_text);


    end % for r (loop over MC runs)


    SMSE(:, p, kpos) = mean(smse)';     % average performance index for each method over mixture realizations


end % for p (loop over iterations)



%%% plot results for current mixture size

figure; set(gcf, 'Name', [figname, ', K = ', k_label]);
for m = 1:nmet
    lin = [linplotK{kpos}, linplotmet{m}];  % line features
    semilogx((ITER*flops_iter(kpos, m) + flops_prewhi(kpos, m))/T, 10*log10(SMSE(m, :, kpos)), lin); hold on
end
grid on
title(read_me_k); xlabel(fig_xlabel); ylabel(fig_ylabel); 
axis(fig_axis); set(gca, 'XTick', fig_xtick);
legend('FastICA', 'pw+FastICA', 'RobustICA', 'pw+RobustICA', 4);

drawnow;

disp(' ');
disp(['FIG. ', num2str(gcf), ': Average extraction quality as a function of computational cost for K = ', num2str(k), ' sources']);
disp(['with signal blocks composed of T = ', num2str(T), ' samples and ', num2str(nruns), ' mixture realizations.']);


% just in case, save current results before starting simulation for new K

clear H_total S_total   % first delete mixing matrix and source realizations (they take too much space)
save(filename)

end % for kpos (loop over number of sources)



%%%%% Generate overall summary figure, with average quality-cost performance for all values of K
           
disp(' ');
disp('----------------------------------------------------------------------------');
disp(' ');

figure; set(gcf, 'Name', figname);
semilogx(1, 2, linplotmet{1}); hold on;     % these are just dummy plots to get the legend right
semilogx(1, 2, linplotmet{2});
semilogx(1, 2, linplotmet{3});
semilogx(1, 2, linplotmet{4});


for kpos = 1:nK

k = K(kpos);

for m = 1:nmet
   lin = [linplotK{kpos}, linplotmet{m}];
   semilogx((ITER*flops_iter(kpos, m) + flops_prewhi(kpos, m))/T, 10*log10(SMSE(m, :, kpos)), lin); hold on
end

end % for kpos

grid on
title(read_me); xlabel(fig_xlabel); ylabel(fig_ylabel); 
axis(fig_axis); set(gca, 'XTick', fig_xtick);
legend('FastICA', 'pw+FastICA', 'RobustICA', 'pw+RobustICA', 4)

disp(['FIG. ', num2str(gcf), ': Average extraction quality as a function of computational cost for different mixture sizes K']);
disp(['with signal blocks composed of T = ', num2str(T), ' samples and ', num2str(nruns), ' mixture realizations.']);
disp(['Solid lines: K = ', num2str(K(1)), '. Dashed lines: K = ', num2str(K(2)), '. Dotted lines: K = ', num2str(K(3)), '.']);
disp(' '); disp(' ');
