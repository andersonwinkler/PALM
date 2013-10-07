function palm_backend(varargin)

% Start by taking what matters from the arguments
global plm opts; % uncomment for debugging
[opts,plm] = palm_takeargs(varargin);

% To store the statistic name for each contrast, to be used later when
% saving the statistic image to a file
plm.Gname = cell(plm.nC,1);

% Variables to store stuff for later.
plm.tmp.rC       = zeros(plm.nC,1);
G                = cell(plm.nY,plm.nC);   % to store G at each permutation (volatile)
df2              = cell(plm.nY,plm.nC);   % to store df2 at each permutation (volatile)
Gpperm           = cell(plm.nY,plm.nC);   % counter, for the permutation p-value (volatile)
plm.G            = cell(plm.nY,plm.nC);   % for the unpermutted G (and to be saved)
plm.df2          = cell(plm.nY,plm.nC);   % for the unpermutted df2 (and to be saved)
plm.Gmax         = cell(plm.nC,1);        % to store the max statistic
plm.nP           = zeros(plm.nC,1);       % number of permutations for each contrast
T                = cell(1,plm.nC);        % to store T at each permutation (volatile)
Tpperm           = cell(1,plm.nC);        % counter, for the combined p-value (volatile)
Tzstat           = cell(1,plm.nC);        % for the combined parametric p-value (volatile)
Tppara           = cell(1,plm.nC);        % for the combined parametric z-score (volatile)
plm.Tmax         = cell(plm.nC,1);        % to store the max combined statistic
if opts.clustere_t.do || opts.clustere_F.do,
    plm.checkcle = zeros(plm.nC,1);       % to store whether do cluster extent for each contrast
    plm.Gcle     = cell(plm.nY,plm.nC);   % to store cluster extent statistic
    plm.Gclemax  = cell(plm.nC,1);        % for the max cluster extent
end
if opts.clusterm_t.do || opts.clusterm_F.do,
    plm.checkclm = zeros(plm.nC,1);       % to store whether do cluster mass for each contrast
    plm.Gclm     = cell(plm.nY,plm.nC);   % to store cluster mass statistic
    plm.Gclmmax  = cell(plm.nC,1);        % for the max cluster mass
end
if opts.tfce.do,
    Gtfce        = cell(plm.nY,plm.nC);   % to store TFCE at each permutation (volatile)
    Gtfcepperm   = cell(plm.nY,plm.nC);   % counter, for the TFCE p-value (volatile)
    plm.Gtfce    = cell(plm.nY,plm.nC);   % for the unpermuted TFCE
    plm.Gtfcemax = cell(plm.nC,1);        % to store the max TFCE statistic
end
if opts.NPC && opts.clustere_npc.do,
    plm.Tcle     = cell(plm.nC,1);        % to store cluster extent NPC statistic
    plm.Tclemax  = cell(plm.nC,1);        % for the max cluster extent NPC
end
if opts.NPC && opts.clusterm_npc.do,
    plm.Tclm     = cell(plm.nC,1);        % to store cluster mass NPC statistic
    plm.Tclmmax  = cell(plm.nC,1);        % for the max cluster mass NPC
end
if opts.NPC && opts.tfce_npc.do,
    Ttfce        = cell(plm.nC,1);        % to store TFCE at each permutation (volatile)
    Ttfcepperm   = cell(plm.nC,1);        % counter, for the TFCE p-value (volatile)
    plm.Ttfce    = cell(plm.nC,1);        % for the unpermuted TFCE
    plm.Ttfcemax = cell(plm.nC,1);        % to store the max TFCE statistic
end

% For each contrast.
for c = 1:plm.nC,
    
    % Partition the model.
    [plm.tmp.X,plm.tmp.Z,plm.tmp.eC] = ...
        palm_partition(plm.M,plm.Cset{c},opts.pmethod);
    
    % Remove from the data and design those observations that are zeroed
    % out in X and belongs to variance groups that are not considered in
    % this contrast.
    idx1 = any(logical(plm.tmp.X),2);
    idx2 = any(bsxfun(@eq,plm.VG,unique(plm.VG(idx1))'),2);
    plm.idx{c} = idx1 | idx2;
    plm.tmp.Yset = cell(plm.nY,1);
    for y = 1:plm.nY,
        plm.tmp.Yset{y} = plm.Yset{y}(plm.idx{c},:);
    end
    plm.tmp.N   = sum(plm.idx{c});
    plm.tmp.X   = plm.tmp.X(plm.idx{c},:);
    plm.tmp.Z   = plm.tmp.Z(plm.idx{c},:);
    plm.tmp.EB  = plm.EB(plm.idx{c},:);
    plm.tmp.VG  = plm.VG(plm.idx{c},:);
    plm.tmp.nVG = numel(unique(plm.tmp.VG));
    
    % Some other variables to be used in the function handles below.
    plm.tmp.Mp = [plm.tmp.X plm.tmp.Z]; % partitioned design matrix, joined
    plm.tmp.rC(c) = rank(plm.tmp.eC);   % rank(C), also df1 for all methods

    % Residual-forming matrix. This is used by the ter Braak method and
    % also to compute some of the stats later. Note that, even though the
    % residual-forming matrix does change at every permutation, the trace
    % for each VG remains unchanged, hence it's not necessary to recompute
    % it for every permutation, and just one works for all.
    plm.tmp.Hm  = plm.tmp.Mp*pinv(plm.tmp.Mp);
    plm.tmp.Rm  = eye(plm.tmp.N) - plm.tmp.Hm;
    plm.tmp.dRm = diag(plm.tmp.Rm); % this is used for the pivotal statistic
    plm.tmp.rM  = plm.tmp.N - sum(plm.tmp.dRm); % this is faster than rank(M)
    
    % Decide which method is going to be used for the regression and
    % permutations, compute some useful matrices for later and create
    % the appropriate function handle to prepare for the model fit.
    % Each of these little functions is a replacement for the generic
    % function 'permglm.m', which is much slower.
    % Note that this swich needs to remain inside the for-loop over
    % contrasts, because the variable plm.tmp, which depends on the
    % partitioning for each contrast, is a constant for these
    % function handles.
    checkterbraak = false;
    switch lower(opts.rmethod),
        case 'exact',
            prepglm         = @(P,Y0)palm_exact(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'draper-stoneman',
            prepglm         = @(P,Y0)palm_draperstoneman(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'still-white',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0)palm_stillwhite(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'freedman-lane',
            plm.tmp.Hz      = plm.tmp.Z*pinv(plm.tmp.Z);
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Hz;
            prepglm         = @(P,Y0)palm_freedmanlane(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'terbraak',
            prepglm         = @(P,Y0)palm_terbraak(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.Mp);
            checkterbraak   = true;
        case 'kennedy',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0,ysel)palm_kennedy(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'manly',
            prepglm         = @(P,Y0)palm_manly(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.Mp);
        case 'huh-jhun',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            [plm.tmp.Q,D]   = schur(plm.tmp.Rz);
            D               = diag(D);
            D               = abs(D) > 10*eps;
            plm.tmp.Q(:,~D) = [];
            prepglm         = @(P,Y0)palm_huhjhun(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
        case 'smith',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0)palm_smith(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.tmp.X);
    end
    
    % To gain in speed, choose an appropriate faster replacement for the
    % function 'pivotal.m', depending on the rank of the contrast and on
    % the presence or not of variance groups. Also define the name of the
    % statistic to save as a file later.
    if     plm.tmp.rC(c) == 1 && plm.tmp.nVG == 1,
        plm.Gname{c} = 'tstat';
        fastpiv = @(M,psi,res)palm_fastt(M,psi,res,plm);
    elseif plm.tmp.rC(c) >  1 && plm.tmp.nVG == 1,
        plm.Gname{c} = 'fstat';
        fastpiv = @(M,psi,res)palm_fastf(M,psi,res,plm,c);
    elseif plm.tmp.rC(c) == 1 && plm.tmp.nVG >  1,
        plm.Gname{c} = 'vstat';
        fastpiv = @(M,psi,res)palm_fastv(M,psi,res,plm);
    elseif plm.tmp.rC(c) >  1 && plm.tmp.nVG >  1,
        plm.Gname{c} = 'gstat';
        fastpiv = @(M,psi,res)palm_fastg(M,psi,res,plm,c);
    end
    
    % Choose the thresholds for the cluster-level.
    plm.checkcle(c) = false;
    plm.checkclm(c) = false;
    if     opts.clustere_t.do && plm.tmp.rC(c) == 1,
        plm.checkcle(c) = true;
        clethr = opts.clustere_t.thr;
    elseif opts.clustere_F.do && plm.tmp.rC(c) >  1,
        plm.checkcle(c) = true;
        clethr = opts.clustere_F.thr;
    end
    if     opts.clusterm_t.do && plm.tmp.rC(c) == 1,
        plm.checkclm(c) = true;
        clmthr = opts.clusterm_t.thr;
    elseif opts.clusterm_F.do && plm.tmp.rC(c) >  1,
        plm.checkclm(c) = true;
        clmthr = opts.clusterm_F.thr;
    end
    
    % Create an appropriate function handle for the NPC. As above, the
    % definition of these handles must stay inside the loop because of the
    % partitioning, which changes for every contrast.
    if opts.NPC,
        checknichols = false;
        plm.Tname = lower(opts.cmethod);
        switch plm.Tname,
            case 'tippett',
                fastnpc      = @(G,df2)palm_tippett(G,df2,plm,c);
                pparanpc     = @(T)palm_tippettp(T,plm);
            case 'fisher',
                fastnpc      = @(G,df2)palm_fisher(G,df2,plm,c);
                pparanpc     = @(T)palm_fisherp(T,plm);
            case 'pearson-david',
                fastnpc      = @(G,df2)palm_pearsondavid(G,df2,plm,c);
                pparanpc     = @(T)palm_pearsondavidp(T,plm);
            case 'stouffer',
                fastnpc      = @(G,df2)palm_stouffer(G,df2,plm,c);
                pparanpc     = @(T)palm_stoufferp(T);
            case 'wilkinson',
                fastnpc      = @(G,df2)palm_wilkinson(G,df2,plm,c);
                pparanpc     = @(T)palm_wilkinsonp(T,plm);
            case 'winer'
                fastnpc      = @(G,df2)palm_winer(G,df2,plm,c);
                pparanpc     = @(T)palm_winerp(T);
            case 'edgington',
                fastnpc      = @(G,df2)palm_edgington(G,df2,plm,c);
                pparanpc     = @(T)palm_edgingtonp(T,plm);
            case 'mudholkar-george',
                fastnpc      = @(G,df2)palm_mudholkargeorge(G,df2,plm,c);
                pparanpc     = @(T)palm_mudholkargeorgep(T,plm);
            case 'friston',
                fastnpc      = @(G,df2)palm_fristonnichols(G,df2,plm,c);
                pparanpc     = @(T)palm_fristonp(T,plm);
            case 'darlington-hayes',
                fastnpc      = @(G,df2)palm_darlingtonhayes(G,df2,plm,c);
            case 'zaykin',
                fastnpc      = @(G,df2)palm_zaykin(G,df2,plm,c);
                pparanpc     = @(T)palm_zaykinp(T,plm);
            case 'dudbridge-koeleman',
                fastnpc      = @(G,df2)palm_dudbridgekoeleman(G,df2,plm,c);
                pparanpc     = @(T)palm_dudbridgekoelemanp(T,plm);
            case 'dudbridge-koeleman2',
                fastnpc      = @(G,df2)palm_dudbridgekoeleman2(G,df2,plm,c);
                pparanpc     = @(T)palm_dudbridgekoeleman2p(T,plm);
            case 'nichols',
                fastnpc      = @(G,df2)palm_fristonnichols(G,df2,plm,c);
                pparanpc     = @(T)palm_nicholsp(T);
                checknichols = true;
            case 'taylor-tibshirani',
                fastnpc      = @(G,df2)palm_taylortibshirani(G,df2,plm,c);
                pparanpc     = @(T)palm_taylortibshiranip(T,plm);
            case 'jiang',
                fastnpc      = @(G,df2)palm_jiang(G,df2,plm,c);
        end
    end
    
    % Define the set of permutations
    if c == 1 || ~any(strcmpi(opts.rmethod,{'terbraak','manly'})),
        [plm.tmp.Pset,plm.nP(c)] = palm_makeperms(opts,plm);
    end
    
    % Some vars for later
    if checkterbraak, psi0 = cell(plm.nY,1); end
    if opts.draft,    ysel = cell(plm.nY,1); end
    plm.Gmax{c} = zeros(plm.nP(c),plm.nY);
    if plm.checkcle(c),  plm.Gclemax{c}  = zeros(plm.nP(c),plm.nY); end
    if plm.checkclm(c),  plm.Gclmmax{c}  = zeros(plm.nP(c),plm.nY); end
    if opts.tfce.do,     plm.Gtfcemax{c} = zeros(plm.nP(c),plm.nY); end
    if opts.clustere_npc.do, plm.Tclemax{c}  = zeros(plm.nP(c),1); end
    if opts.clusterm_npc.do, plm.Tclmmax{c}  = zeros(plm.nP(c),1); end
    if opts.tfce_npc.do,     plm.Ttfcemax{c} = zeros(plm.nP(c),1); end
    if opts.NPC,
        if checknichols,
            plm.Tmax{c} = zeros(plm.nP(c),plm.nY);
        else
            plm.Tmax{c} = zeros(plm.nP(c),1);
        end
    end
    
    % For each permutation
    for p = 1:plm.nP(c),
        
        % Some feedback
        fprintf('Contrast %d/%d, Permutation %d/%d: [ ', ...
            c,plm.nC,p,plm.nP(c));
        
        % For each input dataset
        for y = 1:plm.nY,
            fprintf('%d ',y);
            
            % Shuffle the data and/or design.
            if opts.draft,
                if p == 1,
                    ysel{y} = true(1,size(plm.tmp.Yset{y},2));
                end
                [M,Y] = prepglm(plm.tmp.Pset{p},plm.tmp.Yset{y}(:,ysel{y}));
            else
                [M,Y] = prepglm(plm.tmp.Pset{p},plm.tmp.Yset{y});
            end
            
            % Do the GLM fit.
            psi = pinv(M)*Y; % faster than psi = M\Y
            res = Y - M*psi;
            
            % ter Braak permutes under alternative.
            if checkterbraak,
                if p == 1,
                    psi0{y} = psi;
                else
                    psi = psi - psi0{y};
                end
            end
            
            % Compute the pivotal statistic.
            [G{y,c},df2{y,c}] = fastpiv(M,psi,res);
            
            % In the "draft" mode, the Gpperm variable isn't a counter,
            % but the number of permutations until a statistic larger than
            % the unpermuted was found.
            if opts.draft,
                if p == 1,
                    % In the first permutation, keep G and df2,
                    % and start the counter at 0.
                    plm.G{y,c}   = G{y,c};
                    plm.df2{y,c} = df2{y,c};
                    Gpperm{y,c}  = zeros(size(G{y,c}));
                else
                    % Otherwise, store the permutation in which a larger
                    % statistic happened, and remove this voxel/vertex/face
                    % from further runs.
                    ysel{y}(ysel{y}) = G{y,c} >= plm.G{y,c}(ysel{y});
                    Gpperm{y,c}(ysel{y}) = p;
                    ysel{y} = ~logical(Gpperm{y,c});
                end
            else
                
                % In the first permutation, keep G and df2,
                % and start the counter at 1.
                if p == 1,
                    plm.G{y,c}   = G{y,c};
                    plm.df2{y,c} = df2{y,c};
                    Gpperm{y,c}  = zeros(size(G{y,c}));
                end
                Gpperm{y,c}      = Gpperm{y,c} + (G{y,c} >= plm.G{y,c});
                plm.Gmax{c}(p,y) = max(G{y,c},[],2);
                
                % Cluster extent is here
                if plm.checkcle(c),
                    if p == 1,
                        [plm.Gclemax{c}(p,y),plm.Gcle{y,c}] = ...
                            palm_clustere(G{y,c},y,clethr,opts,plm);
                    else
                        plm.Gclemax{c}(p,y) = ...
                            palm_clustere(G{y,c},y,clethr,opts,plm);
                    end
                end
                
                % Cluster mass is here
                if plm.checkclm(c),
                    if p == 1,
                        [plm.Gclmmax{c}(p,y),plm.Gclm{y,c}] = ...
                            palm_clusterm(G{y,c},y,clmthr,opts,plm);
                    else
                        plm.Gclmmax{c}(p,y) = ...
                            palm_clusterm(G{y,c},y,clmthr,opts,plm);
                    end
                end
                
                % TFCE is here
                if opts.tfce.do,
                    Gtfce{y,c} = palm_tfce(G{y,c},y,opts,plm);
                    if p == 1,
                        plm.Gtfce{y,c} = Gtfce{y,c};
                        Gtfcepperm{y,c} = zeros(size(G{y,c}));
                    end
                    Gtfcepperm{y,c} = Gtfcepperm{y,c} + ...
                        (Gtfce{y,c} >= plm.Gtfce{y,c});
                    plm.Gtfcemax{c}(p,y) = max(Gtfce{y,c},[],2);
                end
            end
        end
        
        % NPC is here
        if opts.NPC,
            fprintf('C ');
            
            % Compute the combined statistic
            Gnpc = cat(1,G{:,c});
            T{c} = fastnpc(Gnpc,cat(1,df2{:,c}));
            
            % Increment counters
            if p == 1,
                plm.T{c}  = T{c};
                Tpperm{c} = zeros(size(T{c}));
            end
            if checknichols,
                Tpperm{c} = Tpperm{c} + ...
                    sum(bsxfun(@ge,Gnpc,plm.T{c}),1);
                plm.Tmax{c}(p,:) = max(T{c},[],2)';
            else
                Tpperm{c} = Tpperm{c} + ...
                    bsxfun(@ge,T{c},plm.T{c});
                plm.Tmax{c}(p) = max(T{c},[],2);
            end
            
            % Just a feedback message for some situations.
            if p == 1                     &&   ...
                    opts.savepara         &&   ...
                    ~ plm.nonpcppara      &&   ...
                    ~ any([                    ...
                    opts.clustere_npc.do       ...
                    opts.clusterm_npc.do       ...
                    opts.tfce_npc.do]')   &&   ...
                    any(strcmpi(opts.cmethod,{ ...
                    'dudbridge-koeleman',      ...
                    'dudbridge-koeleman2'})),
                fprintf('(1st perm is slower) ');
            end
            
            % Since computing the parametric p-value for some methods
            % can be quite slow, it's faster to run all these tests
            % to ensure that 'pparanpc' runs just once.
            if any([ ...
                    opts.clustere_npc.do   ...
                    opts.clusterm_npc.do   ...
                    opts.tfce_npc.do]') || ...
                    (p == 1             && ...
                    opts.savepara       && ...
                    ~ plm.nonpcppara),
                Tppara{c} = pparanpc(T{c});
                
                % Reserve the p-parametric to save later.
                if p == 1,
                    plm.Tppara{c} = Tppara{c};
                end
            end
            
            % Now compute the NPC spatial statistics.
            if any([ ...
                    opts.clustere_npc.do   ...
                    opts.clusterm_npc.do   ...
                    opts.tfce_npc.do]'),
            
                % Convert to z-score.
                Tzstat{c} = norminv(Tppara{c});
                
                % Cluster extent NPC is here
                if opts.clustere_npc.do,
                    if p == 1,
                        [plm.Tclemax{c}(p),plm.Tcle{c}] = ...
                            palm_clustere(Tzstat{c},1, ...
                            opts.clustere_npc.thr,opts,plm);
                    else
                        plm.Tclemax{c}(p) = ...
                            palm_clustere(Tzstat{c},1, ...
                            opts.clustere_npc.thr,opts,plm);
                    end
                end
                
                % Cluster mass NPC is here
                if opts.clusterm_npc.do,
                    if p == 1,
                        [plm.Tclmmax{c}(p),plm.Tclm{c}] = ...
                            palm_clusterm(Tzstat{c},1, ...
                            opts.clusterm_npc.thr,opts,plm);
                    else
                        plm.Tclmmax{c}(p) = ...
                            palm_clusterm(Tzstat{c},1, ...
                            opts.clusterm_npc.thr,opts,plm);
                    end
                end
                
                % TFCE NPC is here
                if opts.tfce_npc.do,
                    Ttfce{c} = palm_tfce(Tzstat{c},1,opts,plm);
                    if p == 1,
                        plm.Ttfce{c} = Ttfce{c};
                        Ttfcepperm{c} = zeros(size(Tzstat{c}));
                    end
                    Ttfcepperm{c} = Ttfcepperm{c} + ...
                        (Ttfce{c} >= plm.Ttfce{c});
                    plm.Ttfcemax{c}(p) = max(Ttfce{c},[],2);
                end
            end 
        end
        fprintf(']\n');
    end
    
    % Save the statistic and the uncorrected & FWER within modality p-values.
    fprintf('Saving p-values (uncorrected and corrected within modality & contrast).\n')
    for y = 1:plm.nY,
        
        % Statistic
        palm_quicksave(plm.G{y,c},opts,plm,y, ...
            sprintf('%s%s_%s_%s_mod%d_con%d', ...
            opts.o,plm.margstr,plm.Ykindstr{y},plm.Gname{c},y,c));

        % Only permutation p-value and its FDR ajustment are saved in the
        % draft mode.
        if opts.draft,
            
            % Permutation p-value, uncorrected
            P = 1./Gpperm{y,c};
            palm_quicksave(1-P,opts,plm,y, ...
                sprintf('%s%s_%s_uncp_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            
            % Permutation p-value, FDR adjusted
            if opts.FDR,
                palm_quicksave(1-palm_fastfdr(P),opts,plm,y, ...
                    sprintf('%s%s_%s_fdrp_mod%d_con%d', ...
                    opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            end
        else
            
            % Permutation p-value
            P = Gpperm{y,c}/plm.nP(c);
            palm_quicksave(1-P,opts,plm,y, ...
                sprintf('%s%s_%s_uncp_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            
            % FWER-corrected within modality and contrast.
            palm_quicksave(1-palm_datacdf(plm.G{y,c},plm.Gmax{c}(:,y)), ...
                opts,plm,y,sprintf('%s%s_%s_fwep_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            
            % Permutation p-value, FDR adjusted
            if opts.FDR,
                palm_quicksave(1-palm_fastfdr(P),opts,plm,y, ...
                    sprintf('%s%s_%s_fdrp_mod%d_con%d', ...
                    opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            end
            
            % Cluster extent results.
            if plm.checkcle(c),
                
                % Cluster extent statistic.
                palm_quicksave(plm.Gcle{y,c},opts,plm,y, ...
                    sprintf('%s%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.margstr,'clustere',plm.Gname{c},y,c));
                
                % Cluster extent FWER p-value
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},plm.Gclemax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwep_mod%d_con%d', ...
                    opts.o,plm.margstr,'clustere',y,c));
            end
            
            % Cluster mass results.
            if plm.checkclm(c),
                
                % Cluster mass statistic.
                palm_quicksave(plm.Gclm{y,c},opts,plm,y, ...
                    sprintf('%s%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.margstr,'clusterm',plm.Gname{c},y,c));
                
                % Cluster mass FWER p-value.
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},plm.Gclmmax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwep_mod%d_con%d', ...
                    opts.o,plm.margstr,'clusterm',y,c));
            end
            
            % TFCE results
            if opts.tfce.do,
                
                % TFCE statistic
                palm_quicksave(plm.Gtfce{y,c},opts,plm,y, ...
                    sprintf('%s%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',plm.Gname{c},y,c));
                
                % TFCE p-value
                P = Gtfcepperm{y,c}/plm.nP(c);
                palm_quicksave(1-P,opts,plm,y, ...
                    sprintf('%s%s_%s_uncp_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',y,c));
                
                % TFCE FWER-corrected within modality and contrast.
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},plm.Gtfcemax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwep_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',y,c));
                
                % TFCE p-value, FDR adjusted.
                if opts.FDR,
                    palm_quicksave(1-palm_fastfdr(P),opts,plm,y, ...
                        sprintf('%s%s_%s_fdrp_mod%d_con%d', ...
                        opts.o,plm.margstr,'tfce',y,c));
                end
            end
        end
        
        % Parametric p-value and its FDR adjustment
        if opts.savepara,
            P = palm_gcdf(plm.G{y,c},plm.tmp.rC(c),plm.df2{y,c});
            palm_quicksave(P, ...
                opts,plm,y,sprintf('%s%s_%s_uncparap_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            if opts.FDR,
                palm_quicksave(1-palm_fastfdr(1-P),opts,plm,y, ...
                    sprintf('%s%s_%s_fdrparap_mod%d_con%d', ...
                    opts.o,plm.margstr,plm.Ykindstr{y},y,c));
            end
        end
    end
    
    % Save the NPC results for this contrast
    if opts.NPC,
        fprintf('Saving NPC p-values (uncorrected and corrected within contrast).\n')
        
        % NPC Statistic
        palm_quicksave(plm.T{c},opts,plm,[], ...
            sprintf('%s%s_%s_%s_con%d', ...
            opts.o,plm.npcstr,plm.Ykindstr{1},plm.Tname,c));
        
        % For the Nichols method, the maxima for all modalities are pooled
        if checknichols,
            plm.Tmax{c} = plm.Tmax{c}(:);
        end
        
        % NPC p-value
        P = Tpperm{c}/numel(plm.Tmax{c});
        palm_quicksave(1-P,opts,plm,[], ...
            sprintf('%s%s_%s_uncp_con%d', ...
            opts.o,plm.npcstr,plm.Ykindstr{1},c));
        
        % NPC FWER-corrected within modality and contrast.
        palm_quicksave(1-palm_datacdf(plm.T{c},plm.Tmax{c}), ...
            opts,plm,[],sprintf('%s%s_%s_fwep_con%d', ...
            opts.o,plm.npcstr,plm.Ykindstr{1},c));
        
        % NPC FDR
        if opts.FDR,
            palm_quicksave(1-palm_fastfdr(P),opts,plm,[], ...
                sprintf('%s%s_%s_fdrp_con%d', ...
                opts.o,plm.npcstr,plm.Ykindstr{1},c));
        end
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara,
            palm_quicksave(plm.Tppara{c},opts,plm,[], ...
                sprintf('%s%s_%s_uncppara_con%d', ...
                opts.o,plm.npcstr,plm.Ykindstr{1},c));
        end
        
        % Cluster extent NPC results.
        if opts.clustere_npc.do,
            
            % Cluster extent statistic.
            palm_quicksave(plm.Tcle{c},opts,plm,[], ...
                sprintf('%s%s_%s_%s_con%d', ...
                opts.o,plm.npcstr,'clustere',plm.Tname,c));
            
            % Cluster extent FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Tcle{c},plm.Tclemax{c}), ...
                opts,plm,y,sprintf('%s%s_%s_fwep_con%d', ...
                opts.o,plm.npcstr,'clustere',c));
        end
        
        % Cluster mass NPC results.
        if opts.clusterm_npc.do,
            
            % Cluster mass statistic.
            palm_quicksave(plm.Tclm{c},opts,plm,[], ...
                sprintf('%s%s_%s_%s_con%d', ...
                opts.o,plm.npcstr,'clusterm',plm.Tname,c));
            
            % Cluster mass FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Tclm{c},plm.Tclmmax{c}), ...
                opts,plm,y,sprintf('%s%s_%s_fwep_con%d', ...
                opts.o,plm.npcstr,'clusterm',c));
        end
        
        % TFCE NPC results.
        if opts.tfce_npc.do,
            
            % TFCE statistic.
            palm_quicksave(plm.Ttfce{c},opts,plm,[], ...
                sprintf('%s%s_%s_%s_con%d', ...
                opts.o,plm.npcstr,'tfce',plm.Tname,c));
            
            % TFCE FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Ttfce{c},plm.Ttfcemax{c}), ...
                opts,plm,y,sprintf('%s%s_%s_fwep_con%d', ...
                opts.o,plm.npcstr,'tfce',c));
        end
    end
end

% Free up some memory after the loop.
clear M Y psi res G df2 Gpperm T Tpperm Tppara Tzstat Ttfce;

% Save FWER corrected across modalities.
if opts.corrmod,
    fprintf('Saving p-values (corrected across modalities).\n')
    for c = 1:plm.nC,
        distmax = max(plm.Gmax{c},[],2);
        for y = 1:plm.nY,
            palm_quicksave(1-palm_datacdf(plm.G{y,c},distmax), ...
                opts,plm,y,sprintf('%s%s_%s_fwepm_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
        end
    end
    
    % Cluster extent
    if all(plm.Yisvol) || all(plm.Yissrf),
        for c = 1:plm.nC,
            if plm.checkcle(c),
                distmax = max(plm.Gclemax{c},[],2);
                for y = 1:plm.nY,
                    palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax), ...
                        opts,plm,y,sprintf('%s%s_%s_fwepm_mod%d_con%d', ...
                        opts.o,plm.margstr,'clustere',y,c));
                end
            end
        end
    end
    
    % Cluster mass
    if all(plm.Yisvol) || all(plm.Yissrf),
        for c = 1:plm.nC,
            if plm.checkclm(c),
                distmax = max(plm.Gclmmax{c},[],2);
                for y = 1:plm.nY,
                    palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},distmax), ...
                        opts,plm,y,sprintf('%s%s_%s_fwepm_mod%d_con%d', ...
                        opts.o,plm.margstr,'clusterm',y,c));
                end
            end
        end
    end
    
    % TFCE
    if opts.tfce.do && (all(plm.Yisvol) || all(plm.Yissrf)),
        for c = 1:plm.nC,
            distmax = max(plm.Gtfcemax{c},[],2);
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},distmax), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepm_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',y,c));
            end
        end
    end
end

% Save FWER corrected across contrasts.
if opts.corrcon,
    fprintf('Saving p-values (corrected across contrasts).\n')
    plm.Gmax = cat(3,plm.Gmax{:});
    distmax = max(plm.Gmax,[],3);
    for c = 1:plm.nC,
        for y = 1:plm.nY,
            palm_quicksave(1-palm_datacdf(plm.G{y,c},distmax(:,y)), ...
                opts,plm,y,sprintf('%s%s_%s_fwepc_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
        end
    end
    
    % Cluster extent
    if all(plm.checkcle) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        plm.Gclemax = cat(3,plm.Gclemax{:});
        distmax = max(plm.Gclemax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepc_mod%d_con%d', ...
                    opts.o,plm.margstr,'clustere',y,c));
            end
        end
    end
    
    % Cluster mass
    if all(plm.checkclm) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        plm.Gclmmax = cat(3,plm.Gclmmax{:});
        distmax = max(plm.Gclmmax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepc_mod%d_con%d', ...
                    opts.o,plm.margstr,'clusterm',y,c));
            end
        end
    end
    
    % TFCE
    if opts.tfce.do && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        plm.Gtfcemax = cat(3,plm.Gtfcemax{:});
        distmax = max(plm.Gtfcemax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepc_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',y,c));
            end
        end
    end
end

% Save FWER corrected across modalities and contrasts.
if opts.corrmod && opts.corrcon,
    fprintf('Saving p-values (corrected across modalities and contrasts).\n')
    distmax = max(max(plm.Gmax,[],3),[],2);
    for c = 1:plm.nC,
        for y = 1:plm.nY,
            palm_quicksave(1-palm_datacdf(plm.G{y,c},distmax), ...
                opts,plm,y,sprintf('%s%s_%s_fwepmc_mod%d_con%d', ...
                opts.o,plm.margstr,plm.Ykindstr{y},y,c));
        end
    end
    
    % Cluster extent
    if all(plm.checkcle) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        distmax = max(max(plm.Gclemax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepmc_mod%d_con%d', ...
                    opts.o,plm.margstr,'clustere',y,c));
            end
        end
    end
    
    % Cluster mass
    if all(plm.checkclm) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        distmax = max(max(plm.Gclmmax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},distmax), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepmc_mod%d_con%d', ...
                    opts.o,plm.margstr,'clusterm',y,c));
            end
        end
    end
    
    % TFCE
    if opts.tfce.do && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.tmp.rC)) == 1,
        distmax = max(max(plm.Gtfcemax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},distmax), ...
                    opts,plm,y,sprintf('%s%s_%s_fwepmc_mod%d_con%d', ...
                    opts.o,plm.margstr,'tfce',y,c));
            end
        end
    end
end

% Save FWER corrected across contrasts for NPC.
if opts.NPC && opts.corrcon,
    fprintf('Saving NPC p-values (corrected across contrasts).\n')
    plm.Tmax = cat(2,plm.Tmax{:});
    distmax = max(plm.Tmax,[],2);
    for c = 1:plm.nC,
        palm_quicksave(1-palm_datacdf(plm.T{c},distmax), ...
            opts,plm,[],sprintf('%s%s_%s_fwepc_con%d', ...
            opts.o,plm.npcstr,plm.Ykindstr{1},c));
    end
    
    % Cluster extent NPC
    if opts.clustere_npc.do,
        plm.Tclemax = cat(3,plm.Tclemax{:});
        distmax = max(plm.Tclemax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Tcle{c},distmax), ...
                opts,plm,[],sprintf('%s%s_%s_fwepc_con%d', ...
                opts.o,plm.npcstr,'clustere',c));
        end
    end
    
    % Cluster mass NPC
    if opts.clusterm_npc.do,
        plm.Tclmmax = cat(3,plm.Tclmmax{:});
        distmax = max(plm.Tclmmax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Tclm{c},distmax), ...
                opts,plm,[],sprintf('%s%s_%s_fwepc_con%d', ...
                opts.o,plm.npcstr,'clusterm',c));
        end
    end
    
    % TFCE NPC
    if opts.tfce_npc.do,
        plm.Ttfcemax = cat(3,plm.Ttfcemax{:});
        distmax = max(plm.Ttfcemax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Ttfce{c},distmax), ...
                opts,plm,[],sprintf('%s%s_%s_fwepc_con%d', ...
                opts.o,plm.npcstr,'tfce',c));
        end
    end
end

% Finished! :-)