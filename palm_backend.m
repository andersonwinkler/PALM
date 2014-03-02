function palm_backend(varargin)

% Start by taking what matters from the arguments
global plm opts; % uncomment for debugging
[opts,plm] = palm_takeargs(varargin{1});

% To store the statistic name for each contrast, to be used later when
% saving the statistic image to a file
plm.Gname = cell(plm.nC,1);

% Variables to store stuff for later.
plm.rC           = zeros(plm.nC,1);       % to store the rank of each contrast
G                = cell(plm.nY,plm.nC);   % to store G at each permutation (volatile)
df2              = cell(plm.nY,plm.nC);   % to store df2 at each permutation (volatile)
Gpperm           = cell(plm.nY,plm.nC);   % counter, for the permutation p-value (volatile)
if opts.draft,
    Gppermp      = Gpperm;                % number of perms done, for the draft mode
end
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
    
    % Partition the model, now using the method chosen by the user
    [plm.tmp.X,plm.tmp.Z,plm.tmp.eCm,plm.tmp.eCx] = ...
        palm_partition(plm.M,plm.Cset{c},opts.pmethod);
    
    % Some methods don't work well if Z is empty
    if isempty(plm.tmp.Z),
        plm.tmp.rmethod = 'noz';
    else
        plm.tmp.rmethod = opts.rmethod;
    end
    
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
    plm.rC(c)  = rank(plm.tmp.eCm);     % rank(C), also df1 for all methods

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
    switch lower(plm.tmp.rmethod),
        case 'noz',
            prepglm         = @(P,Y0)noz(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCx;
        case 'exact',
            prepglm         = @(P,Y0)exact(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCx;
        case 'draper-stoneman',
            prepglm         = @(P,Y0)draperstoneman(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCm;
        case 'still-white',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0)stillwhite(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCx;
        case 'freedman-lane',
            plm.tmp.Hz      = plm.tmp.Z*pinv(plm.tmp.Z);
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Hz;
            prepglm         = @(P,Y0)freedmanlane(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCm;
        case 'terbraak',
            checkterbraak   = true;
            prepglm         = @(P,Y0)terbraak(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCm;
        case 'kennedy',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0,ysel)kennedy(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCx;
        case 'manly',
            prepglm         = @(P,Y0)manly(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCm;
        case 'huh-jhun',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            [plm.tmp.Q,D]   = schur(plm.tmp.Rz);
            D               = diag(D);
            D               = abs(D) > 10*eps;
            plm.tmp.Q(:,~D) = [];
            prepglm         = @(P,Y0)huhjhun(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCx;
        case 'smith',
            plm.tmp.Rz      = eye(plm.tmp.N) - plm.tmp.Z*pinv(plm.tmp.Z);
            prepglm         = @(P,Y0)smith(P,Y0,plm);
            plm.tmp.seq     = palm_mat2seq(plm.Xset{c});
            plm.tmp.eC      = plm.tmp.eCm;
    end
    
    % To gain in speed, choose an appropriate faster replacement for the
    % function 'pivotal.m', depending on the rank of the contrast and on
    % the presence or not of variance groups. Also define the name of the
    % statistic to save as a file later.
    if     plm.rC(c) == 1 && plm.tmp.nVG == 1,
        plm.Gname{c} = 'tstat';
        fastpiv = @(M,psi,res)fastt(M,psi,res,plm);
    elseif plm.rC(c) >  1 && plm.tmp.nVG == 1,
        plm.Gname{c} = 'fstat';
        fastpiv = @(M,psi,res)fastf(M,psi,res,plm,c);
    elseif plm.rC(c) == 1 && plm.tmp.nVG >  1,
        plm.Gname{c} = 'vstat';
        fastpiv = @(M,psi,res)fastv(M,psi,res,plm);
    elseif plm.rC(c) >  1 && plm.tmp.nVG >  1,
        plm.Gname{c} = 'gstat';
        fastpiv = @(M,psi,res)fastg(M,psi,res,plm,c);
    end
    
    % Choose the thresholds for the cluster-level.
    plm.checkcle(c) = false;
    plm.checkclm(c) = false;
    if     opts.clustere_t.do && plm.rC(c) == 1,
        plm.checkcle(c) = true;
        clethr = opts.clustere_t.thr;
    elseif opts.clustere_F.do && plm.rC(c) >  1,
        plm.checkcle(c) = true;
        clethr = opts.clustere_F.thr;
    end
    if     opts.clusterm_t.do && plm.rC(c) == 1,
        plm.checkclm(c) = true;
        clmthr = opts.clusterm_t.thr;
    elseif opts.clusterm_F.do && plm.rC(c) >  1,
        plm.checkclm(c) = true;
        clmthr = opts.clusterm_F.thr;
    end
    
    % Create an appropriate function handle for the NPC. As above, the
    % definition of these handles must stay inside the loop because of the
    % partitioning, which changes for every contrast.
    if opts.NPC,
        chknichols = false;
        plm.Tname = lower(opts.cmethod);
        switch plm.Tname,
            case 'tippett',
                fastnpc    = @(G,df2)tippett(G,df2,plm,c);
                pparanpc   = @(T)tippettp(T,plm);
            case 'fisher',
                fastnpc    = @(G,df2)fisher(G,df2,plm,c);
                pparanpc   = @(T)fisherp(T,plm);
            case 'pearson-david',
                fastnpc    = @(G,df2)pearsondavid(G,df2,plm,c);
                pparanpc   = @(T)pearsondavidp(T,plm);
            case 'stouffer',
                fastnpc    = @(G,df2)stouffer(G,df2,plm,c);
                pparanpc   = @(T)stoufferp(T);
            case 'wilkinson',
                fastnpc    = @(G,df2)wilkinson(G,df2,plm,c);
                pparanpc   = @(T)wilkinsonp(T,plm);
            case 'winer'
                fastnpc    = @(G,df2)winer(G,df2,plm,c);
                pparanpc   = @(T)winerp(T);
            case 'edgington',
                fastnpc    = @(G,df2)edgington(G,df2,plm,c);
                pparanpc   = @(T)edgingtonp(T,plm);
            case 'mudholkar-george',
                fastnpc    = @(G,df2)mudholkargeorge(G,df2,plm,c);
                pparanpc   = @(T)mudholkargeorgep(T,plm);
            case 'friston',
                fastnpc    = @(G,df2)fristonnichols(G,df2,plm,c);
                pparanpc   = @(T)fristonp(T,plm);
            case 'darlington-hayes',
                fastnpc    = @(G,df2)darlingtonhayes(G,df2,plm,c);
            case 'zaykin',
                fastnpc    = @(G,df2)zaykin(G,df2,plm,c);
                pparanpc   = @(T)zaykinp(T,plm);
            case 'dudbridge-koeleman',
                fastnpc    = @(G,df2)dudbridgekoeleman(G,df2,plm,c);
                pparanpc   = @(T)dudbridgekoelemanp(T,plm);
            case 'dudbridge-koeleman2',
                fastnpc    = @(G,df2)dudbridgekoeleman2(G,df2,plm,c);
                pparanpc   = @(T)dudbridgekoeleman2p(T,plm);
            case 'nichols',
                fastnpc    = @(G,df2)fristonnichols(G,df2,plm,c);
                pparanpc   = @(T)nicholsp(T);
                chknichols = true;
            case 'taylor-tibshirani',
                fastnpc    = @(G,df2)taylortibshirani(G,df2,plm,c);
                pparanpc   = @(T)taylortibshiranip(T,plm);
            case 'jiang',
                fastnpc    = @(G,df2)jiang(G,df2,plm,c);
        end
    end
    
    % Define the set of permutations
    if c == 1 || ~any(strcmpi(opts.rmethod,{'terbraak','manly'})),
        [plm.tmp.Pset,plm.nP(c)] = palm_shuftree(opts,plm);
    end

    % If the user wants to save the permutations, save the vectors now.
    % This has 3 benefits: (1) the if-test below will run just once, rather
    % than many times inside the loop, (2) if the user only wants the
    % vectors, not the images, he/she can cancel immediately after the
    % text file has been created and (3) having all just as a single big
    % file is more convenient than hundreds of small ones.
    if opts.saveperms,
        fid = fopen(sprintf('%s_%s_%s_mod%d_con%d_permidx.csv', ...
            opts.o,plm.Ykindstr{y},plm.Gname{c},y,c),'w');
        for p = 1:plm.nP(c),
            fprintf(fid,'%d,',palm_perm2idx(plm.tmp.Pset{p})');
            fseek(fid,-1,'eof');
            fprintf(fid,'\n');
        end
        fclose(fid);
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
        if chknichols,
            plm.Tmax{c} = zeros(plm.nP(c),plm.nY);
        else
            plm.Tmax{c} = zeros(plm.nP(c),1);
        end
    end
    
    % For each permutation
    for p = 1:plm.nP(c),
        
        % Some feedback
        fprintf('Contrast %d/%d, Shuffling %d/%d: [ ', ...
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
                    Gppermp{y,c} = zeros(size(G{y,c}));
                    
                    % Save the degrees of freedom
                    csvwrite(sprintf('%s_%s_%s_mod%d_con%d.dof', ...
                        opts.o,plm.Ykindstr{y},plm.Gname{c},y,c), ...
                        [plm.rC(c) plm.df2{y,c}]);
                else
                    % Otherwise, store the permutation in which a larger
                    % statistic happened, and remove this voxel/vertex/face
                    % from further runs.
                    Gpperm{y,c}(ysel{y}) = Gpperm{y,c}(ysel{y}) + ...
                        (G{y,c} >= plm.G{y,c}(ysel{y}));
                    Gppermp{y,c}(ysel{y}) = p;
                    ysel{y} = Gpperm{y,c} < opts.draft;
                end
            else
                if p == 1,
                    % In the first permutation, keep G and df2,
                    % and start the counter at 1.
                    plm.G{y,c}   = G{y,c};
                    plm.df2{y,c} = df2{y,c};
                    Gpperm{y,c}  = zeros(size(G{y,c}));
                    
                    % Save the degrees of freedom
                    csvwrite(sprintf('%s_%s_%s_mod%d_con%d.dof', ...
                        opts.o,plm.Ykindstr{y},plm.Gname{c},y,c), ...
                        [plm.rC(c) plm.df2{y,c}]);
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
                
                % If the user wants to save the statistic for each
                % permutation, save it now. This isn't obviously allowed
                % in draft mode, as the images are not complete.
                if opts.saveperms,
                    palm_quicksave(G{y,c},opts,plm,y, ...
                        sprintf('%s_%s_%s_mod%d_con%d_perm%06d', ...
                        opts.o,plm.Ykindstr{y},plm.Gname{c},y,c,p));
                end
            end
        end
        
        % NPC is here
        if opts.NPC,
            fprintf('C ');
            
            % Compute the combined statistic
            Gnpc = cat(1,G{:,c});
            T{c} = fastnpc(Gnpc,cat(1,df2{:,c}));
            
            % If the user wants to save the NPC statistic for each
            % permutation, save it now.
            if opts.saveperms,
                palm_quicksave(T{c},opts,plm,[], ...
                    sprintf('%s_%s_%s_con%d_perm%06d', ...
                    opts.o,plm.Ykindstr{1},'npc',c,p));
            end
            
            % Increment counters
            if p == 1,
                plm.T{c}  = T{c};
                Tpperm{c} = zeros(size(T{c}));
            end
            if chknichols,
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
            sprintf('%s_%s_%s_mod%d_con%d', ...
            opts.o,plm.Ykindstr{y},plm.Gname{c},y,c));

        % Only permutation p-value and its FDR ajustment are saved in the
        % draft mode.
        if opts.draft,
            
            % Permutation p-value, uncorrected
            P = (Gpperm{y,c}+1)./Gppermp{y,c}; % it could be (Gpperm{y,c}+1)./Gppermp{y,c};, slightly more conservative
            palm_quicksave(1-P,opts,plm,y, ...
                sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'uncp',y,c));
            
            % Permutation p-value, FDR adjusted
            if opts.FDR,
                palm_quicksave(1-fastfdr(P),opts,plm,y, ...
                    sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.Ykindstr{y},plm.Gname{c},'fdrp',y,c));
            end
        else
            
            % Permutation p-value
            P = Gpperm{y,c}/plm.nP(c);
            palm_quicksave(1-P,opts,plm,y, ...
                sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'uncp',y,c));
            
            % FWER-corrected within modality and contrast.
            palm_quicksave(1-palm_datacdf(plm.G{y,c},plm.Gmax{c}(:,y)), ...
                opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'fwep',y,c));
            
            % Permutation p-value, FDR adjusted
            if opts.FDR,
                palm_quicksave(1-fastfdr(P),opts,plm,y, ...
                    sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.Ykindstr{y},plm.Gname{c},'fdrp',y,c));
            end
            
            % Cluster extent results.
            if plm.checkcle(c),
                
                % Cluster extent statistic.
                palm_quicksave(plm.Gcle{y,c},opts,plm,y, ...
                    sprintf('%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clustere',plm.Gname{c},y,c));
                
                % Cluster extent FWER p-value
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},plm.Gclemax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clustere',plm.Gname{c},'fwep',y,c));
            end
            
            % Cluster mass results.
            if plm.checkclm(c),
                
                % Cluster mass statistic.
                palm_quicksave(plm.Gclm{y,c},opts,plm,y, ...
                    sprintf('%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clusterm',plm.Gname{c},y,c));
                
                % Cluster mass FWER p-value.
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},plm.Gclmmax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clusterm',plm.Gname{c},'fwep',y,c));
            end
            
            % TFCE results
            if opts.tfce.do,
                
                % TFCE statistic
                palm_quicksave(plm.Gtfce{y,c},opts,plm,y, ...
                    sprintf('%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},y,c));
                
                % TFCE p-value
                P = Gtfcepperm{y,c}/plm.nP(c);
                palm_quicksave(1-P,opts,plm,y, ...
                    sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},'uncp',y,c));
                
                % TFCE FWER-corrected within modality and contrast.
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},plm.Gtfcemax{c}(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},'fwep',y,c));
                
                % TFCE p-value, FDR adjusted.
                if opts.FDR,
                    palm_quicksave(1-fastfdr(P),opts,plm,y, ...
                        sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                        opts.o,'tfce',plm.Gname{c},'fdrp',y,c));
                end
            end
        end
        
        % Parametric p-value and its FDR adjustment
        if opts.savepara,
            P = palm_gcdf(plm.G{y,c},plm.rC(c),plm.df2{y,c});
            palm_quicksave(P, ...
                opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'uncparap',y,c));
            if opts.FDR,
                palm_quicksave(1-fastfdr(1-P),opts,plm,y, ...
                    sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,plm.Ykindstr{y},plm.Gname{c},'fdrparap',y,c));
            end
        end
    end
    
    % Save the NPC results for this contrast
    if opts.NPC,
        fprintf('Saving NPC p-values (uncorrected and corrected within contrast).\n')
        
        % NPC Statistic
        palm_quicksave(plm.T{c},opts,plm,[], ...
            sprintf('%s_%s_%s_con%d', ...
            opts.o,plm.Ykindstr{1},'npc',c));
        
        % For the Nichols method, the maxima for all modalities are pooled
        if chknichols,
            plm.Tmax{c} = plm.Tmax{c}(:);
        end
        
        % NPC p-value
        P = Tpperm{c}/numel(plm.Tmax{c});
        palm_quicksave(1-P,opts,plm,[], ...
            sprintf('%s_%s_%s_%s_con%d', ...
            opts.o,plm.Ykindstr{1},'npc','uncp',c));
        
        % NPC FWER-corrected within modality and contrast.
        palm_quicksave(1-palm_datacdf(plm.T{c},plm.Tmax{c}), ...
            opts,plm,[],sprintf('%s_%s_%s_%s_con%d', ...
            opts.o,plm.Ykindstr{1},'npc','fwep',c));
        
        % NPC FDR
        if opts.FDR,
            palm_quicksave(1-fastfdr(P),opts,plm,[], ...
                sprintf('%s_%s_%s_%s_con%d', ...
                opts.o,plm.Ykindstr{1},'npc','fdrp',c));
        end
        
        % Parametric combined pvalue
        if opts.savepara && ~ plm.nonpcppara,
            palm_quicksave(plm.Tppara{c},opts,plm,[], ...
                sprintf('%s_%s_%s_%s_con%d', ...
                opts.o,plm.Ykindstr{1},'npc','uncparap',c));
        end
        
        % Cluster extent NPC results.
        if opts.clustere_npc.do,
            
            % Cluster extent statistic.
            palm_quicksave(plm.Tcle{c},opts,plm,[], ...
                sprintf('%s_%s_%s_con%d', ...
                opts.o,'clustere','npc',c));
            
            % Cluster extent FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Tcle{c},plm.Tclemax{c}), ...
                opts,plm,y,sprintf('%s_%s_%s_%s_con%d', ...
                opts.o,'clustere','npc','fwep',c));
        end
        
        % Cluster mass NPC results.
        if opts.clusterm_npc.do,
            
            % Cluster mass statistic.
            palm_quicksave(plm.Tclm{c},opts,plm,[], ...
                sprintf('%s_%s_%s_con%d', ...
                opts.o,'clusterm','npc',c));
            
            % Cluster mass FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Tclm{c},plm.Tclmmax{c}), ...
                opts,plm,y,sprintf('%s_%s_%s_%s_con%d', ...
                opts.o,'clusterm','npc','fwep',c));
        end
        
        % TFCE NPC results.
        if opts.tfce_npc.do,
            
            % TFCE statistic.
            palm_quicksave(plm.Ttfce{c},opts,plm,[], ...
                sprintf('%s_%s_%s_con%d', ...
                opts.o,'tfce','npc',c));
            
            % TFCE FWER p-value
            palm_quicksave(1-palm_datacdf(plm.Ttfce{c},plm.Ttfcemax{c}), ...
                opts,plm,y,sprintf('%s_%s_%s_%s_con%d', ...
                opts.o,'tfce','npc','fwep',c));
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
                opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'fwemp',y,c));
        end
    end
    
    % Cluster extent
    if all(plm.Yisvol) || all(plm.Yissrf),
        for c = 1:plm.nC,
            if plm.checkcle(c),
                distmax = max(plm.Gclemax{c},[],2);
                for y = 1:plm.nY,
                    palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax), ...
                        opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                        opts.o,'clustere',plm.Gname{c},'fwemp',y,c));
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
                        opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                        opts.o,'clusterm',plm.Gname{c},'fwemp',y,c));
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
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},'fwemp',y,c));
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
                opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'fwecp',y,c));
        end
    end
    
    % Cluster extent
    if all(plm.checkcle) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        plm.Gclemax = cat(3,plm.Gclemax{:});
        distmax = max(plm.Gclemax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clustere',plm.Gname{c},'fwecp',y,c));
            end
        end
    end
    
    % Cluster mass
    if all(plm.checkclm) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        plm.Gclmmax = cat(3,plm.Gclmmax{:});
        distmax = max(plm.Gclmmax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clusterm',plm.Gname{c},'fwecp',y,c));
            end
        end
    end
    
    % TFCE
    if opts.tfce.do && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        plm.Gtfcemax = cat(3,plm.Gtfcemax{:});
        distmax = max(plm.Gtfcemax,[],3);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},distmax(:,y)), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},'fwecp',y,c));
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
                opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                opts.o,plm.Ykindstr{y},plm.Gname{c},'fwecmp',y,c));
        end
    end
    
    % Cluster extent
    if all(plm.checkcle) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        distmax = max(max(plm.Gclemax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gcle{y,c},distmax), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clustere',plm.Gname{c},'fwecmp',y,c));
            end
        end
    end
    
    % Cluster mass
    if all(plm.checkclm) && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        distmax = max(max(plm.Gclmmax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gclm{y,c},distmax), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'clusterm',plm.Gname{c},'fwecmp',y,c));
            end
        end
    end
    
    % TFCE
    if opts.tfce.do && ...
            (all(plm.Yisvol) || all(plm.Yissrf)) && ...
            numel(unique(plm.rC)) == 1,
        distmax = max(max(plm.Gtfcemax,[],3),[],2);
        for c = 1:plm.nC,
            for y = 1:plm.nY,
                palm_quicksave(1-palm_datacdf(plm.Gtfce{y,c},distmax), ...
                    opts,plm,y,sprintf('%s_%s_%s_%s_mod%d_con%d', ...
                    opts.o,'tfce',plm.Gname{c},'fwecmp',y,c));
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
            opts,plm,[],sprintf('%s_%s_%s_con%d', ...
            opts.o,plm.Ykindstr{1},'npc','fwecp',c));
    end
    
    % Cluster extent NPC
    if opts.clustere_npc.do,
        plm.Tclemax = cat(3,plm.Tclemax{:});
        distmax = max(plm.Tclemax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Tcle{c},distmax), ...
                opts,plm,[],sprintf('%s_%s_%s_con%d', ...
                opts.o,'clustere','npc','fwecp',c));
        end
    end
    
    % Cluster mass NPC
    if opts.clusterm_npc.do,
        plm.Tclmmax = cat(3,plm.Tclmmax{:});
        distmax = max(plm.Tclmmax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Tclm{c},distmax), ...
                opts,plm,[],sprintf('%s_%s_%s_con%d', ...
                opts.o,'clusterm','npc','fwecp',c));
        end
    end
    
    % TFCE NPC
    if opts.tfce_npc.do,
        plm.Ttfcemax = cat(3,plm.Ttfcemax{:});
        distmax = max(plm.Ttfcemax,[],3);
        for c = 1:plm.nC,
            palm_quicksave(1-palm_datacdf(plm.Ttfce{c},distmax), ...
                opts,plm,[],sprintf('%s_%s_%s_con%d', ...
                opts.o,'tfce','npc','fwecp',c));
        end
    end
end

% ==============================================================
% Below are the functions for each of the regression methods:
% ==============================================================

function [Mr,Y] = noz(P,Y,plm)
% Note Y remains unchanged
% This is the same as Draper-Stoneman
Mr = P*plm.tmp.X;

function [Mr,Yr] = exact(P,Y,plm)
% The "exact" method, in which the coefficients for
% the nuisance are known.
Yr = Y - plm.tmp.Z*plm.g;
Mr = P*plm.tmp.X;

function [Mr,Y] = draperstoneman(P,Y,plm)
% Draper and Stoneman (1966) method.
% Note Y remains unchanged
Mr = [P*plm.tmp.X plm.tmp.Z];

function [Mr,Yr] = stillwhite(P,Y,plm)
% A method following the same logic as the one
% proposed by Still and White (1981)
Yr = plm.tmp.Rz*Y;
Mr = P*plm.tmp.X;

function [Mr,Yr] = freedmanlane(P,Y,plm)
% The Freedman and Lane (1983) method.
Mr = [plm.tmp.X plm.tmp.Z];
Yr = (P'*plm.tmp.Rz + plm.tmp.Hz)*Y;

function [Mr,Yr] = terbraak(P,Y,plm)
% The ter Braak (1992) method.
Mr = [plm.tmp.X plm.tmp.Z];
Yr = (P'*plm.tmp.Rm + plm.tmp.Hm)*Y; % original
% Yr = P'*plm.tmp.Rm*Y; % alternative (causes unpermuted statistic to be 0)

function [Mr,Yr] = kennedy(P,Y,plm)
% The Kennedy (1996) method. This method should NEVER be used.
Mr = plm.tmp.Rz*plm.tmp.X;
Yr = P'*plm.tmp.Rz*Y;

function [Mr,Yr] = manly(P,Y,plm)
% The Manly (1997) method.
Mr = [plm.tmp.X plm.tmp.Z];
Yr = P'*Y;

function [Mr,Yr] = huhjhun(P,Y,plm)
% The Huh and Jhun (2001) method, that fixes the issues
% with Kennedy's, but doesn't allow block permutation.
Mr = plm.tmp.Q'*plm.tmp.Rz*plm.tmp.X;
Yr = P'*plm.tmp.Q'*plm.tmp.Rz*Y;

function [Mr,Y] = smith(P,Y,plm)
% The Smith method, i.e., orthogonalization.
% Note Y remains unchanged
Mr = [P*plm.tmp.Rz*plm.tmp.X plm.tmp.Z];

% ==============================================================
% Below are the functions to compute the statistics:
% ==============================================================
function [G,df2] = fastt(M,psi,res,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
% 
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm.m' and 'takeargs.m'
% 
% Outputs:
% G   : t statistic.
% df2 : Degrees of freedom. df1 is 1 for the t statistic.

df2 = plm.tmp.N-plm.tmp.rM;
G = plm.tmp.eC'*psi;
den = sqrt(plm.tmp.eC'*inv(M'*M)*plm.tmp.eC*sum(res.^2)./df2); %#ok inv here
G = G./den;

% ==============================================================
function [G,df2] = fastf(M,psi,res,plm,c)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups = 1
% 
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm.m' and 'takeargs.m'
% 
% Outputs:
% G   : F-statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).

df2 = plm.tmp.N-plm.tmp.rM;
cte = plm.tmp.eC*inv(plm.tmp.eC'*inv(M'*M)*plm.tmp.eC)*plm.tmp.eC'; %#ok inv here

tmp = zeros(size(psi));
for j = 1:size(cte,2),
    tmp(j,:) = sum(bsxfun(@times,psi,cte(:,j)))';
end
G = sum(tmp.*psi);
ete = sum(res.^2);
G = G./ete * (df2/plm.rC(c));

% ==============================================================
function [G,df2] = fastv(M,psi,res,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups > 1
% 
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm.m' and 'takeargs.m'
% 
% Outputs:
% G   : Aspin-Welch v statistic.
% df2 : Degrees of freedom 2. df1 is 1.

r = size(M,2);
m = size(res,2);
[~,~,blkidx] = unique(plm.tmp.VG);

W = zeros(plm.tmp.nVG,m);
dRmb = zeros(plm.tmp.nVG,1);
cte = zeros(r^2,m);
for b = 1:plm.tmp.nVG,
    bidx = blkidx == b;
    dRmb(b) = sum(plm.tmp.dRm(bidx));
    W(b,:) = dRmb(b)./sum(res(bidx,:).^2);
    Mb = M(bidx,:)'*M(bidx,:);
    cte = cte + Mb(:)*W(b,:);
    W(b,:) = W(b,:)*sum(bidx);
end

den = zeros(1,m);
for t = 1:m,
    den(t) = plm.tmp.eC'*inv(reshape(cte(:,t),[r r]))*plm.tmp.eC; %#ok inv here
end
G = plm.tmp.eC'*psi./sqrt(den);

bsum = zeros(1,m);
sW1 = sum(W,1);
for b = 1:plm.tmp.nVG,
    bsum = bsum + bsxfun(@rdivide,(1-W(b,:)./sW1).^2,dRmb(b));
end
df2 = 1/3./bsum;

% ==============================================================
function [G,df2] = fastg(M,psi,res,plm,c)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups > 1
% 
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm.m' and 'takeargs.m'
% 
% Outputs:
% G   : Welch v^2 statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).

r = size(M,2);
m = size(res,2);
[~,~,blkidx] = unique(plm.tmp.VG);

W = zeros(plm.tmp.nVG,m);
dRmb = zeros(plm.tmp.nVG,1);
cte = zeros(r^2,m);
for b = 1:plm.tmp.nVG,
    bidx = blkidx == b;
    dRmb(b) = sum(plm.tmp.dRm(bidx));
    W(b,:) = dRmb(b)./sum(res(bidx,:).^2);
    Mb = M(bidx,:)'*M(bidx,:);
    cte = cte + Mb(:)*W(b,:);
    W(b,:) = W(b,:)*sum(bidx);
end

G = zeros(1,m);
for t = 1:m,
    A = psi(:,t)'*plm.tmp.eC;
    G(t) = A*inv(plm.tmp.eC'*inv(reshape(cte(:,t),[r r]))* ...
        plm.tmp.eC)*A'/plm.rC(c); %#ok inv here
end

bsum = zeros(1,m);
sW1 = sum(W,1);
for b = 1:plm.tmp.nVG,
    bsum = bsum + bsxfun(@rdivide,(1-W(b,:)./sW1).^2,dRmb(b));
end
bsum = bsum/plm.rC(c)/(plm.rC(c)+2);
df2 = 1/3./bsum;
G = G./(1 + 2*(plm.rC(c)-1).*bsum);

% ==============================================================
% Below are the functions to combine statistics:
% ==============================================================
% See the original combine.m for commented lines on implementation

function T = tippett(G,df2,plm,c)
T = max(palm_gcdf(G,plm.rC(c),df2),[],1);

function P = tippettp(T,plm)
P = T.^plm.nY;

% ==============================================================
function T = fisher(G,df2,plm,c)
T = -2*sum(log(1-palm_gcdf(G,plm.rC(c),df2)),1);

function P = fisherp(T,plm)
P = chi2cdf(T,2*plm.nY);

% ==============================================================
function T = pearsondavid(G,df2,plm,c)
P = palm_gcdf(G,plm.rC(c),df2);
T = -2*min(sum(log(1-P),1),sum(log(P),1));

function P = pearsondavidp(T,plm)
P = chi2cdf(T,2*plm.nY);

% ==============================================================
function T = stouffer(G,df2,plm,c)
T = sum(norminv(palm_gcdf(G,plm.rC(c),df2)))/sqrt(plm.nY);

function P = stoufferp(T)
P = normcdf(T);

% ==============================================================
function T = wilkinson(G,df2,plm,c)
T = sum((1-palm_gcdf(G,plm.rC(c),df2)) <= plm.npcparm);

function P = wilkinsonp(T,plm)
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log(f) + lfac(f);
end
lalpha  = log(plm.npcparm);
l1alpha = log(1-plm.npcparm);
P = zeros(size(T));
for k = 1:plm.nY,
    lp1 = lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1);
    lp2 = k*lalpha;
    lp3 = (plm.nY-k)*l1alpha;
    P = P + (k>=T).*exp(lp1+lp2+lp3);
end
P = 1-P;

% ==============================================================
function T = winer(G,df2,plm,c)
df2 = bsxfun(@times,ones(size(G)),df2);
cte = sqrt(sum(df2./(df2-2),1));
T   = sum(tinv(palm_gcdf(G,plm.rC(c),df2),df2))./cte;

function P = winerp(T)
P = normcdf(T);

% ==============================================================
function T = edgington(G,df2,plm,c)
T = -log10(sum(1-palm_gcdf(G,plm.rC(c),df2),1));

function P = edgingtonp(T,plm)
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log(f) + lfac(f);
end
T = 10.^(-T);
fT   = floor(T);
mxfT = max(fT(:));
P = zeros(size(T));
for j = 0:mxfT,
    p1  = (-1)^j;
    lp2 = - lfac(j+1) - lfac(plm.nY-j+1);
    lp3 = plm.nY*log(T-j);
    P = P + (j<=fT).*p1.*exp(lp2+lp3);
end
P = 1-P;

% ==============================================================
function T = mudholkargeorge(G,df2,plm,c)
P = palm_gcdf(G,plm.rC(c),df2);
mhcte = sqrt(3*(5*plm.nY+4)/plm.nY/(5*plm.nY+2))/pi;
T = mhcte*sum(log(P./(1-P)),1);

function P = mudholkargeorgep(T,plm)
P = tcdf(T,5*plm.nY+4);

% ==============================================================
function T = fristonnichols(G,df2,plm,c)
T = min(palm_gcdf(G,plm.rC(c),df2),[],1);

function P = fristonp(T,plm)
P = 1-(1-T).^plm.nY;

function T = nicholsp(T)
% T itself is P, so there is nothing to do.

% ==============================================================
function T = darlingtonhayes(G,df2,plm,c)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= plm.npcparm;
G       = reshape(G(idx),  [plm.npcparm size(G,2)]);
df2     = reshape(df2(idx),[plm.npcparm size(df2,2)]);
P       = palm_gcdf(G,plm.rC(c),df2);
Z       = norminv(P);
T       = mean(Z,1);

% ==============================================================
function T = zaykin(G,df2,plm,c)
P = -log10(1-palm_gcdf(G,plm.rC(c),df2));
P(P < -log10(plm.npcparm)) = 0;
T = sum(P,1);

function P = zaykinp(T,plm)
lT = -T;
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log10(f) + lfac(f);
end
lalpha   = log10(plm.npcparm);
l1alpha  = log10(1-plm.npcparm);
P = zeros(size(lT));
for k = 1:plm.nY,
    lp1 = lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1);
    lp2 = (plm.nY-k)*l1alpha;
    Tsmall = lT <= k*lalpha;
    Tlarge = ~ Tsmall;
    p3 = 0;
    lnum = log10(k*lalpha - lT(Tsmall));
    for j = 1:k,
        p3 = p3 + 10.^(lT(Tsmall) + (j-1).*lnum - lfac(j));
    end
    lp3small = log10(p3);
    lp3large = k*lalpha;
    P(Tsmall) = P(Tsmall) + 10.^(lp1 + lp2 + lp3small);
    P(Tlarge) = P(Tlarge) + 10.^(lp1 + lp2 + lp3large);
end
P = 1-P;

% ==============================================================
function T = dudbridgekoeleman(G,df2,plm,c)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= plm.npcparm;
G       = reshape(G(idx),  [plm.npcparm size(G,2)]);
df2     = reshape(df2(idx),[plm.npcparm size(df2,2)]);
P       = -log10(1-palm_gcdf(G,plm.rC(c),df2));
T       = sum(P,1);

function P = dudbridgekoelemanp(T,plm)
lT = -T;
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log10(f) + lfac(f);
end
P = zeros(size(lT));
lp1 = lfac(plm.nY+1) -         ...
    lfac(plm.npcparm+2) -      ...
    lfac(plm.nY-plm.npcparm) + ...
    log10(plm.npcparm+2);
for v = 1:numel(lT);
    P(v) = quad(@(t)dkint(t,lp1,lT(v),plm.nY,...
        plm.npcparm,lfac(1:plm.npcparm)),eps,1);
end
P = 1-P;

function T = dudbridgekoeleman2(G,df2,plm,c)
df2 = bsxfun(@times,ones(size(G)),df2);
P = -log10(1-palm_gcdf(G,plm.rC(c),df2));
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
P(tmp > plm.npcparm) = 0;
P(P < -log10(plm.npcparm2)) = 0;
T = sum(P,1);

function P = dudbridgekoeleman2p(T,plm)
lT = -T;
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log10(f) + lfac(f);
end
P = zeros(1,size(T,2));
for k = 1:plm.npcparm,
    kk = (plm.nY-k)*log(1-plm.npcparm2);
    if isnan(kk), kk = 0; end
    p1 = exp(lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1) + kk);
    p2 = awtk(lT,plm.npcparm2,k,lfac(1:k));
    P = P + p1.*p2;
end
if k < plm.nY,
    lp1 = lfac(plm.nY+1) -         ...
        lfac(plm.npcparm+2) -      ...
        lfac(plm.nY-plm.npcparm) + ...
        log(plm.npcparm+2);
    for v = 1:numel(lT);
        P(v) = P(v) + ...
            quad(@(t)dkint(t,lp1,lT(v),plm.nY,plm.npcparm, ...
            lfac(1:plm.npcparm)),eps,plm.npcparm2);
    end
end
P = 1-P;

function q = dkint(t,lp1,lT,K,r,lfac)
lp2 = (K-r-1).*log(1-t);
ltr = r.*log(t);
L1  = real(lp1 + lp2 + ltr);
s1  = (lT > ltr).*exp(L1);
j   = (1:r)';
lp3 = lT + (j-1)*log(r*log(t)-lT) ...
    - repmat(lfac(j),[1 numel(t)]);
L2  = real(lp1 + repmat(lp2,[r 1]) + lp3);
s2  = (lT <= ltr).*sum(exp(L2));
q   = s1 + s2;

function A = awtk(lw,t,k,lfac)
ltk = k.*log(t);
tk = real(exp(ltk));
s = (1:k)';
L = bsxfun(@plus,lw,...
    bsxfun(@minus,(s-1)*log(k*log(t)-lw),lfac(s)));
S = sum(real(exp(L)),1);
A = (lw <= ltk).*S + (lw > ltk).*tk;

% ==============================================================
function T = taylortibshirani(G,df2,plm,c)
P = 1-palm_gcdf(G,plm.rC(c),df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum(1-P.*(plm.nY+1)./prank)/plm.nY;

function P = taylortibshiranip(T,plm)
P = normcdf(T,0,1/sqrt(plm.nY));

% ==============================================================
function T = jiang(G,df2,plm,c)
P = 1-palm_gcdf(G,plm.rC(c),df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum((P<=plm.npcparm).*(1-P.*(plm.nY+1)./prank))/plm.nY;

% ==============================================================
% Below are other useful functions:
% ==============================================================
function padj = fastfdr(pval)
V = numel(pval);
[pval,oidx] = sort(pval);
[~,oidxR]   = sort(oidx);
padj = zeros(size(pval));
prev = 1;
for i = V:-1:1,
    padj(i) = min(prev,pval(i)*V/i);
    prev = padj(i);
end
padj = padj(oidxR);

% Finished! :-)