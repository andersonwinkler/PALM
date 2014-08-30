function [opts,plm] = palm_takeargs(varargin)

% Load the defaults
opts = palm_defaults;

% As varargin is actually from another function, fix it.
if nargin == 1,
    if exist(varargin{1},'file'),
        vararginx = palm_configrw(varargin{1});
    else
        error('Unknown option or file not found: %s',varargin{1});
    end
else
    vararginx = varargin;
    idxa = find(strcmpi(vararginx,'-o'));
    if isempty(idxa),
        cfgname = horzcat(opts.o,'_palmconfig.txt');
    else
        cfgname = horzcat(vararginx{idxa+1},'_palmconfig.txt');
    end
    [opth,~,~] = fileparts(cfgname);
    if ~isempty(opth) && ~exist(opth,'dir'),
        mkdir(opth);
    end
    palm_configrw(vararginx,cfgname);
end

% Number of input images/masks/surfaces
Ni = sum(strcmp(vararginx,'-i'));  % number of data inputs
Nm = sum(strcmp(vararginx,'-m'));  % number of masks
Ns = sum(strcmp(vararginx,'-s'));  % number of surfaces
Nd = sum(strcmp(vararginx,'-d'));  % number of design files (currently only one can actually be used)
Nt = sum(strcmp(vararginx,'-t'));  % number of t-contrast files
Nf = sum(strcmp(vararginx,'-f'));  % number of F-test files
opts.i   = cell(Ni,1);  % Input files (to constitute Y later)
opts.m   = cell(Nm,1);  % Mask file(s)
opts.s   = cell(Ns,1);  % Surface file(s)
opts.d   = cell(Nd,1);  % Design file(s)
opts.t   = cell(Nt,1);  % t contrast file(s)
opts.f   = cell(Nf,1);  % F contrast file(s)
opts.eb       = [];     % File with definition of exchangeability blocks
opts.vg       = [];     % File with definition of variance groups
opts.EE       = false;  % To be filled below (don't edit this here!)
opts.ISE      = false;  % To be filled below (don't edit this here!)
opts.within   = false;  % To be filled below (don't edit this here!)
opts.whole    = false;  % To be filled below (don't edit this here!)
opts.conlist  = [];     % File with the list of contrasts to be performed
opts.conskipcount = 0;  % When saving the contrasts, skip how many from 1?
opts.singlevg = true;   % Make sure that sigle VG will be used if nothing is supplied (don't edit this here!)

% These are to be incremented below
a = 1; i = 1; m = 1; d = 1;
t = 1; f = 1; s = 1;

% Remove trailing empty arguments. This is useful for some Octave versions.
while numel(vararginx) > 0 && isempty(vararginx{1}),
    vararginx(1) = [];
end
narginx = numel(vararginx);

% Take the input arguments
while a <= narginx,
    switch vararginx{a},
        case {'-help','-?','-basic','-advanced'},
            
            % Do nothing, as these options are parsed separately,
            % and should anyway be given without any other argument.
            a = a + 1;
            
        case '-i',
            
            % Get the filenames for the data.
            opts.i{i} = vararginx{a+1};
            i = i + 1;
            a = a + 2;
            
        case '-m',
            
            % Get the filenames for the masks, if any.
            opts.m{m} = vararginx{a+1};
            m = m + 1;
            a = a + 2;
            
        case '-s',
            
            % Get the filenames for the surfaces, if any.
            opts.s{s} = vararginx{a+1};
            s = s + 1;
            a = a + 2;
            
        case '-d',
            
            % Get the design matrix file.
            opts.d{d} = vararginx{a+1};
            d = d + 1;
            a = a + 2;
            
        case '-t',
            
            % Get the t contrast files.
            opts.t{t} = vararginx{a+1};
            t = t + 1;
            a = a + 2;
            
        case '-f',
            
            % Get the F contrast files.
            opts.f{f} = vararginx{a+1};
            f = f + 1;
            a = a + 2;
            
        case '-conlist',
            
            % List of contrasts that will actually be performed, from a
            % text file
            opts.conlist = varargin{a+1};
            a = a + 2;
            
        case '-conskipcount',
            
            % Numbers to skip when saving the contrasts
            opts.conskipcount = vararginx{a+1};
            if ischar(opts.conskipcount),
                opts.conskipcount = str2double(opts.concountskip);
            end
            a = a + 2;
            
        case '-fonly',
            
            % Run only the F-contrasts
            opts.fonly = true;
            a = a + 1;
            
        case '-eb',
            
            % Get the exchangeability blocks file.
            opts.eb = vararginx{a+1};
            a = a + 2;
            
        case '-vg'
            
            % Get the variance groups file.
            opts.vg = vararginx{a+1};
            opts.singlevg = false;
            if ischar(opts.vg) && ...
                    ~any(strcmpi(opts.vg,{'auto','automatic','default'})),
                opts.vg = 'auto';
            end
            a = a + 2;
            
        case '-o',
            
            % Output prefix for the files to be saved.
            opts.o = vararginx{a+1};
            a = a + 2;
            
        case '-n',
            
            % Number of permutations
            opts.nP0 = vararginx{a+1};
            if ischar(opts.nP0),
                opts.nP0 = str2double(opts.nP0);
            end
            a = a + 2;
            
        case '-c',
            
            % Threshold for cluster extent, univariate, NPC and MV
            opts.clustere_uni.do  = true;
            opts.clustere_uni.thr = vararginx{a+1};
            if ischar(opts.clustere_uni.thr),
                opts.clustere_uni.thr = str2double(opts.clustere_uni.thr);
            end
            opts.clustere_npc.do  = true;
            opts.clustere_npc.thr = vararginx{a+1};
            if ischar(opts.clustere_npc.thr),
                opts.clustere_npc.thr = str2double(opts.clustere_npc.thr);
            end
            opts.clustere_mv.do   = true;
            opts.clustere_mv.thr  = vararginx{a+1};
            if ischar(opts.clustere_mv.thr),
                opts.clustere_mv.thr = str2double(opts.clustere_mv.thr);
            end
            a = a + 2;
            
        case '-C',
            
            % Threshold for cluster mass, univariate, NPC and MV
            opts.clusterm_uni.do  = true;
            opts.clusterm_uni.thr = vararginx{a+1};
            if ischar(opts.clusterm_uni.thr),
                opts.clusterm_uni.thr = str2double(opts.clusterm_uni.thr);
            end
            opts.clusterm_npc.do  = true;
            opts.clusterm_npc.thr = vararginx{a+1};
            if ischar(opts.clusterm_npc.thr),
                opts.clusterm_npc.thr = str2double(opts.clusterm_npc.thr);
            end
            opts.clusterm_mv.do   = true;
            opts.clusterm_mv.thr  = vararginx{a+1};
            if ischar(opts.clusterm_mv.thr),
                opts.clusterm_mv.thr = str2double(opts.clusterm_mv.thr);
            end
            a = a + 2;
            
        case '-cuni',
            
            % Threshold for cluster extent, univariate
            opts.clustere_uni.do  = true;
            opts.clustere_uni.thr = vararginx{a+1};
            if ischar(opts.clustere_uni.thr),
                opts.clustere_uni.thr = str2double(opts.clustere_uni.thr);
            end
            a = a + 2;
            
        case '-Cuni',
            
            % Threshold for cluster mass, univariate
            opts.clusterm_uni.do  = true;
            opts.clusterm_uni.thr = vararginx{a+1};
            if ischar(opts.clusterm_uni.thr),
                opts.clusterm_uni.thr = str2double(opts.clusterm_uni.thr);
            end
            a = a + 2;
            
        case '-cnpc',
            
            % Threshold for cluster extent, NPC
            opts.NPC = true;
            opts.clustere_npc.do  = true;
            opts.clustere_npc.thr = vararginx{a+1};
            if ischar(opts.clustere_npc.thr),
                opts.clustere_npc.thr = str2double(opts.clustere_npc.thr);
            end
            a = a + 2;
            
        case '-Cnpc',
            
            % Threshold for cluster mass, NPC
            opts.NPC = true;
            opts.clusterm_npc.do  = true;
            opts.clusterm_npc.thr = vararginx{a+1};
            if ischar(opts.clusterm_npc.thr),
                opts.clusterm_npc.thr = str2double(opts.clusterm_npc.thr);
            end
            a = a + 2;
            
        case '-cmv',
            
            % Threshold for cluster extent, MV
            opts.MV = true;
            opts.clustere_mv.do  = true;
            opts.clustere_mv.thr = vararginx{a+1};
            if ischar(opts.clustere_mv.thr),
                opts.clustere_mv.thr = str2double(opts.clustere_mv.thr);
            end
            a = a + 2;
            
        case '-Cmv',
            
            % Threshold for cluster mass, MV
            opts.MV = true;
            opts.clusterm_mv.do  = true;
            opts.clusterm_mv.thr = vararginx{a+1};
            if ischar(opts.clusterm_mv.thr),
                opts.clusterm_mv.thr = str2double(opts.clusterm_mv.thr);
            end
            a = a + 2;
            
        case '-T',
            
            % Do TFCE for univariate, NPC and MV?
            opts.tfce_uni.do = true;
            opts.tfce_npc.do = true;
            opts.tfce_mv.do  = true;
            a = a + 1;
            
        case '-Tuni',
            
            % Do TFCE for uni?
            opts.tfce_uni.do = true;
            a = a + 1;
            
        case '-Tnpc',
            
            % Do TFCE for NPC?
            opts.NPC = true;
            opts.tfce_npc.do = true;
            a = a + 1;
            
        case '-Tmv',
            
            % Do TFCE for MV?
            opts.MV = true;
            opts.tfce_mv.do  = true;
            a = a + 1;
            
        case '-tfce2D',
            
            % Shortcut for -tfce_H 2 -tfce_E 1 -tfce_C 26,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 1;
            opts.tfce.conn   = 26;
            a = a + 1;
            
        case '-tfce_H',
            
            % TFCE H parameter
            opts.tfce.H = vararginx{a+1};
            if ischar(opts.tfce.H),
                opts.tfce.H = str2double(opts.tfce.H);
            end
            a = a + 2;
            
        case '-tfce_E',
            
            % TFCE E parameter
            opts.tfce.E = vararginx{a+1};
            if ischar(opts.tfce.E),
                opts.tfce.E = str2double(opts.tfce.E);
            end
            a = a + 2;
            
        case '-tfce_C',
            
            % TFCE connectivity
            opts.tfce.conn = vararginx{a+1};
            if ischar(opts.tfce.conn),
                opts.tfce.conn = str2double(opts.tfce.conn);
            end
            a = a + 2;
            
        case '-within',
            
            % Define whether should permute blocks as a whole or not
            opts.within = true;
            a = a + 1;
            
        case '-whole',
            
            % Define whether should permute blocks as a whole or not
            opts.whole = true;
            a = a + 1;
            
        case '-ee',
            
            % Exchangeable errors (EE)?
            % If yes, this means permutations.
            opts.EE = true;
            a = a + 1;
            
        case '-ise',
            
            % Independent and symmetric errors (ISE)?
            % If yes, this means sign-flippings.
            opts.ISE = true;
            a = a + 1;
            
        case '-cmcp',
            
            % Define whether Conditional Monte Carlo should be used or not,
            % that is, ignoring repeated elements in the permutation set.
            opts.cmcp = true;
            a = a + 1;
            
        case '-cmcx',
            
            % Define whether repeated rows in X should be ignored or not
            % when defining the permutations, which constitutes another
            % form of CMC
            opts.cmcx = true;
            a = a + 1;
            
        case '-twotail',
            
            % Do a two-tailed test for all t-contrasts?
            opts.twotail = true;
            a = a + 1;
            
        case '-corrmod',
            
            % Correct over modalities.
            opts.corrmod = true;
            a = a + 1;
            
        case '-corrcon',
            
            % Correct over contrasts. In this case, the shuffling method may
            % need to change to ter Braak or Manly, depending on the contrasts.
            opts.corrcon = true;
            a = a + 1;
            
        case '-saveparametric',
            
            % If the user wants to have also the parametric p-values
            opts.savepara = true;
            a = a + 1;
            
        case {'-savecdf','-save1-p'},
            
            % Save 1-p values (CDF) instead of the P-values
            opts.savecdf = true;
            a = a + 1;
            
        case '-logp',
            
            % Convert the P-values or (1-P)-values to -log10 before saving
            opts.savelogp = true;
            a = a + 1;
            
        case '-savemask',
            
            % If the user wants to have also the masks used for each
            % modality
            opts.savemask = true;
            a = a + 1;
            
        case '-rmethod',
            
            % Which method to use for the regression/permutation?
            if nargin > a,
                methlist = {           ...
                    'Draper-Stoneman', ...
                    'Still-White',     ...
                    'Freedman-Lane',   ...
                    'terBraak',        ...
                    'Kennedy',         ... % should never be used
                    'Manly',           ...
                    'Huh-Jhun',        ...
                    'Smith'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Regression/Permutation method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.rmethod = methlist{methidx};
            else
                error([...
                    'The option -rmethod requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
            
        case '-npc',
            
            % Do the non-parametric combination?
            opts.NPC = true;
            if nargin == a,
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1}(1),'-'),
                a = a + 1;
                
            elseif nargin > a,
                
                % Which combining function to use for the combination?
                methlist = {               ...
                    'Tippett',             ...
                    'Fisher',              ...
                    'Pearson-David',       ...
                    'Stouffer',            ...
                    'Wilkinson',           ...
                    'Winer',               ...
                    'Edgington',           ...
                    'Mudholkar-George',    ...
                    'Friston',             ...
                    'Darlington-Hayes',    ...
                    'Zaykin',              ...
                    'Dudbridge-Koeleman',  ...
                    'Dudbridge-Koeleman2', ...
                    'Nichols',             ...
                    'Taylor-Tibshirani',   ...
                    'Jiang'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx);
                    error('Combining method "%s" unknown.',vararginx{a+1});
                elseif any(strcmpi(vararginx{a+1},{...
                        'Wilkinson',       ...
                        'Zaykin',          ...
                        'Jiang'})),
                    if ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 0.05;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        a = a + 3;
                        plm.npcparm = eval(vararginx{a+2});
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif any(strcmpi(vararginx{a+1},{...
                        'Darlington-Hayes',   ...
                        'Dudbridge-Koeleman', ...
                        'Jiang'})),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Friston'),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        plm.npcparm = vararginx{a+2};
                        a = a + 3;
                    end
                elseif strcmpi(vararginx{a+1},'Dudbridge-Koeleman2'),
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        plm.npcparm  = 1;
                        plm.npcparm2 = 0.05;
                        a = a + 2;
                    else
                        if ischar(vararginx{a+2}),
                            plm.npcparm = eval(vararginx{a+2});
                        else
                            plm.npcparm = vararginx{a+2};
                        end
                        if nargin == a + 2 || ...
                                ischar(vararginx{a+3}) && ...
                                strcmpi(vararginx{a+3}(1),'-'),
                            plm.npcparm2 = 0.05;
                        elseif ischar(vararginx{a+3}),
                            plm.npcparm2 = eval(vararginx{a+3});
                        else
                            plm.npcparm2 = vararginx{a+3};
                        end
                        a = a + 4;
                    end
                else
                    a = a + 2;
                end
                opts.npcmethod = methlist{methidx};
            end
            
        case '-mv',
            
            % Compute classic multivariate statistics
            opts.MV = true;
            if nargin == a,
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1}(1),'-'),
                a = a + 1;
                
            elseif nargin > a,
                
                % Which multivariate statistic to use?
                methlist = {     ...
                    'Wilks',     ...
                    'Hotelling', ...
                    'Pillai',    ...
                    'Roy',       ...
                    'Roy_ii',    ...
                    'Roy_iii'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if any(methidx);
                    opts.mvstat = methlist{methidx};
                    a = a + 2;
                else
                    error('Multivariate statistic "%s" unknown.',vararginx{a+1});
                end
            end
            
        case '-fdr',
            
            % Compute FDR-adjusted p-values
            opts.FDR = true;
            a = a + 1;
            
        case '-draft',
            
            % Do a draft scheme
            opts.draft = vararginx{a+1};
            if ischar(opts.draft),
                opts.draft = str2double(opts.draft);
            end
            a = a + 2;
            
        case '-noniiclass',
            
            % Disable using the NIFTI class
            opts.useniiclass = false;
            a = a + 1;
            
        case '-saveperms',
            
            % Allow no use of mask for 4D NIFTI files
            opts.saveperms = true;
            a = a + 1;
            
        case '-savemetrics',
            
            % Save a file with the number of permutations, average
            % Hamming distance, etc.
            opts.savemetrics = true;
            a = a + 1;
            
        case '-inormal',
            
            % Inverse-normal transformation?
            opts.inormal = true;
            
            % Take the parameters given to -inormal
            parms = {};
            if nargin - a >= 1,
                if ~strcmp(vararginx{a+1}(1),'-'),
                    parms{1} = vararginx{a+1};
                end
            end
            if nargin - a >= 2,
                if ~strcmp(vararginx{a+2}(1),'-'),
                    parms{2} = vararginx{a+2};
                end
            end
            a = a + 1 + numel(parms);
            
            % Then modify the variables accordingly
            methlist = {   ...
                'Blom',    ...
                'Tukey',   ...
                'Bliss',   ...
                'Waerden', ...
                'SOLAR'};
            for p = 1:numel(parms),
                methidx = strcmpi(parms{p},methlist);
                if any(methidx);
                    opts.inormal_meth = parms{p};
                elseif any(parms{p},{'quali','qualitative','discrete'}),
                    opts.inormal_quanti = false;
                elseif any(parms{p},{'quanti','quantitative','continuous'}),
                    opts.inormal_quanti = true;
                else
                    error('Parameter "%s" unknown for the -inormal option.',parms{p});
                end
            end
            
        case '-seed'
            
            % Seed for the random number generator
            opts.seed = vararginx{a+1};
            if ischar(opts.seed) && ...
                    ~any(strcmpi(opts.seed,{'shuffle','twist','reset'})),
                opts.seed = str2double(opts.seed);
            end
            a = a + 2;
            
        case '-demean',
            
            % Demean data and design. Additionally, remove
            % a global intercept, if any, from the design.
            opts.demean = true;
            a = a + 1;
            
        case '-vgdemean',
            
            % Demean data and design within VG. Additionally, remove
            % a global intercept, if any, from the design.
            opts.vgdemean = true;
            a = a + 1;
            
        case '-ev4vg',
            
            % Add to the design matrix one EV for each variance group.
            opts.ev4vg = true;
            a = a + 1;
            
        case '-removeignored',
            
            % Remove from the analysis observations that are have their own
            % exclusive regressor and don't belong to any contrast (i.e,
            % always nuisance.
            opts.removeignored = true;
            a = a + 1;
            
        case '-removevgbysize',
            
            % Remove from the analysis observations that are the only
            % in their variance group.
            opts.removevgbysize = vararginx{a+1};
            if ischar(opts.removevgbysize),
                opts.removevgbysize = str2double(opts.removevgbysize);
            end
            a = a + 2;
            
        case '-zstat',
            
            % Convert the statistic for each test to a z-score
            opts.zstat = true;
            a = a + 1;
            
        case '-pearson',
            
            % Compute the Pearson's correlation coefficient (R^2 if rank(C)>1).
            opts.pearson = true;
            a = a + 1;
            
        case '-noranktest',
            
            % Don't test the rank(Y) before doing MANOVA/MANCOVA.
            opts.noranktest = true;
            a = a + 1;
            
        case '-evperdat',
            
            % Use one (single) EV per datum
            opts.evperdat = true;
            a = a + 1;
            
        case '-transposedata',
            
            % Transpose the data if it's 2D?
            opts.transposedata = true;
            a = a + 1;
            
        case '-pmethod1', % not in the help
            
            % Which method to use to partition the model when defining
            % the permutations?
            if nargin > a,
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Winkler',  ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethod1 = methlist{methidx};
            else
                error([...
                    'The option -pmethod1 requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        case '-pmethod2', % not in the help
            
            % Which method to use to partition the model when defining
            % doing the actual regression?
            if nargin > a,
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Winkler',  ...
                    'Ridgway',  ...
                    'none'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethod2 = methlist{methidx};
            else
                error([...
                    'The option -pmethod2 requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        case '-highestH', % not in the help
            
            % Use only the perms with the highest Hamming distance?
            opts.highestH = vararginx{a+1};
            if ischar(opts.highestH),
                opts.highestH = str2double(opts.highestH);
            end
            a = a + 2;
            
        case '-lowestH', % not in the help
            
            % Use only the perms with the lowest Hamming distance?
            opts.lowestH = vararginx{a+1};
            if ischar(opts.lowestH),
                opts.lowestH = str2double(opts.lowestH);
            end
            a = a + 2;
            
        otherwise
            error('Unknown option: ''%s''',vararginx{a});
    end
end

% Until voxelwise regressors are implemented, this should give an error
if Nd > 1,
    error('Only one design file is currently allowed.\n');
end

% A quick check for the case of 1 EV per column in Y.
if opts.evperdat,
    if any([...
            opts.cmcx
            opts.removeignored
            opts.ev4vg
            opts.pearson]),
        error([...
            'The option ''-evperdat'' is incompatible with the options listed below:\n' ...
            '''-igrepx''\n' ...
            '''-removeignored''\n' ...
            '''-ev4vg''\n' ...
            '''-pearson''\n' ...
            'None of these can be enaabled with ''-evperdat''.%s'],'');
    end
    if strcmpi(opts.rmethod,'terBraak'),
        error('The option ''-evperdat'' cannot be used with the ter Braak method (not implemented)');
    end
end

% Only highest or lowest Hamming allowed, not both
if ~isempty(opts.highestH) && ~isempty(opts.lowestH),
    error('Only one of ''-highestH'' or ''-lowestH'' can be used at any given time.');
end
if ~isempty(opts.highestH) && (opts.highestH <= 0 || opts.highestH > 1),
    error('The value given to ''-highestH'' must be >= 0 or < 1');
end
if ~isempty(opts.lowestH)  && (opts.lowestH  <= 0 || opts.lowestH  > 1),
    error('The value given to ''-lowestH'' must be >= 0 or < 1');
end

% This simplifies some tests later
opts.spatial     = false;
opts.spatial_uni = false;
opts.spatial_npc = false;
opts.spatial_mv  = false;
if any([ ...
        opts.clustere_uni.do  ...
        opts.clusterm_uni.do  ...
        opts.tfce_uni.do      ...
        opts.clustere_npc.do  ...
        opts.clusterm_npc.do  ...
        opts.tfce_npc.do      ...
        opts.clustere_mv.do   ...
        opts.clusterm_mv.do   ...
        opts.tfce_mv.do]'),
    opts.spatial = true;
    if any([ ...
            opts.clustere_uni.do  ...
            opts.clusterm_uni.do  ...
            opts.tfce_uni.do]'),
        opts.spatial_uni = true;
    end
    if any([ ...
            opts.clustere_npc.do  ...
            opts.clusterm_npc.do  ...
            opts.tfce_npc.do]'),
        opts.spatial_npc = true;
    end
    if any([ ...
            opts.clustere_mv.do   ...
            opts.clusterm_mv.do   ...
            opts.tfce_mv.do]'),
        opts.spatial_mv = true;
    end
end

if opts.spatial && palm_isoctave,
    pkg load image
end

% No FWER or NPC if using draft mode
if opts.draft,
    if opts.corrmod || opts.corrcon,
        warning('The draft mode does not allow FWER-correction, only FDR.\n%s',''); %#ok
    end
    if opts.NPC,
        warning('The draft mode does not allow NPC.\n%s',''); %#ok
    end
    if opts.clustere_npc.do || opts.clustere_npc.do || opts.tfce_npc.do,
        warning('The draft mode does not allow spatial statistics (cluster or TFCE).\n%s',''); %#ok
    end
    opts.corrmod         = false;
    opts.corrcon         = false;
    opts.NPC             = false;
    opts.clustere_npc.do = false;
    opts.clusterm_npc.do = false;
    opts.tfce_npc.do     = false;
    opts.clustere_mv.do  = false;
    opts.clusterm_mv.do  = false;
    opts.tfce_mv.do      = false;
end

% Some NPC methods don't have an analytical form for the parametric p-value
if opts.NPC && any(strcmpi(opts.npcmethod,{'Darlington-Hayes','Jiang'})),
    plm.nonpcppara = true;
    if opts.savepara,
        warning([...
            'No parametric combination p-value will be saved for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
    end
    if opts.spatial_npc,
        warning([ ...
            'No NPC cluster-level or TFCE statistic will be produced for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
        opts.clustere_npc.do = false;
        opts.clusterm_npc.do = false;
        opts.tfce_npc.do     = false;
        opts.spatial_npc     = false;
    end
else
    plm.nonpcppara = false;
end

% Likewise, some MV methods don't have an analytical form for the parametric p-value
if opts.MV && strcmpi(opts.mvstat,'Roy_iii'),
    plm.nomvppara = true;
    if opts.savepara,
        warning('No parametric p-value will be saved for the Roy_iii method.%s',''); %#ok
    end
    if opts.spatial_mv,
        warning([ ...
            'No multivariate cluster-level or TFCE statistic will be produced\n', ...
            '         for the Roy_iii statistic%s'],'');
        opts.clustere_mv.do  = false;
        opts.clusterm_mv.do  = false;
        opts.tfce_mv.do      = false;
        opts.spatial_mv      = false;
    end
else
    plm.nomvppara = false;
end

% Some more warnings and sanity checks
if opts.pearson && (opts.NPC || opts.MV),
    warning([ ...
        'It''s not possible to compute the Pearson''s r or R^2 together with NPC or\n', ...
        '         multivariate methods. Disabling the options ''-npc'' and ''-mv''.%s'],'');
    opts.NPC = false;
    opts.MV  = false;
end
if opts.pearson && ~ opts.demean,
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the data and the design\n' ...
        '         must be mean centered. Adding option ''-demean''.%s'],'');
    opts.demean = true;
end
if opts.pearson && ~ any(strcmpi(opts.pmethod2,{'beckmann','ridgway'})),
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option ''-pmethod2 Beckmann''.%s'],'');
    opts.pmethod2 = 'beckmann';
end
if opts.demean && opts.vgdemean && ~ opts.pearson,
    warning([...
        'Cannot use the option ''-demean'' together with ''-vgdemean''\n'...
        '         Ignoring the option ''-vgdemean''.%s'],'');
    opts.vgdemean = false;
end
if opts.ev4vg && opts.vgdemean,
    warning([...
        'Cannot use the option ''-ev4vg'' together with ''-vgdemean''\n'...
        '         Ignoring the option ''-ev4vg''%s.'],'');
    opts.ev4vg = false;
end

% Initialize the random number generator
if palm_isoctave,
    if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
        opts.seed = 'reset';
    end
    rand('state',opts.seed); %#ok
else
    if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
        opts.seed = 'shuffle';
    end
    rng(opts.seed);
end

% Read and organise the surfaces. If no surfaces have been loaded, but the
% user wants cluster extent, cluster mass, or TFCE, an error will be
% printed later down in the code.
if opts.spatial && Ns > 0,
    plm.srf = cell(Ns,1);
    for s = 1:Ns,
        plm.srf{s} = palm_miscread(opts.s{s});
    end
end

% There should be no more masks than modalities, and the number of
% masks needs to be either 1 or the same number of modalities.
if Nm > Ni,
    error([...
        'There are more masks supplied with -m (%d masks) than\n'...
        'modalities supplied with -i (%d modalities)'],Nm,Ni);
elseif Nm > 1 && Nm ~= Ni,
    error([...
        'The number of masks supplied with -m (%d masks) is larger than 1,\n'...
        'but still not the same as the number of modalities supplied with\n'...
        'the option -i (%d modalities).'],Nm,Ni);
end

% Read and organise the masks. If there are no masks specified, one for
% each modality will be created after each modality is loaded.
plm.masks = cell(Ni,1);
for m = 1:Nm,
    plm.masks{m} = palm_miscread(opts.m{m},opts.useniiclass);
    if strcmp(plm.masks{m}.readwith,'nifticlass'),
        plm.masks{m}.data = double(plm.masks{m}.data);
    end
    plm.masks{m}.data(isnan(plm.masks{m}.data)) = 0;
    plm.masks{m}.data(isinf(plm.masks{m}.data)) = 0;
    plm.masks{m}.data = logical(plm.masks{m}.data);
end
if Nm == 1,
    for i = 2:Ni,
        plm.masks{i} = plm.masks{1};
    end
end

% Read and organise the data.
plm.Yset     = cell(Ni,1);  % Regressands (Y)
plm.Yisvol   = false(Ni,1); % Is Y a volume image?
plm.Yissrf   = false(Ni,1); % Is Y a surface-based image (DPX)?
plm.Yisvtx   = false(Ns,1); % Is vertexwise?
plm.Yisfac   = false(Ns,1); % is facewise? (this is currently dichotomous with Yisvtx, but later there may be edges/connectivity too)
plm.Yarea    = cell(Ns,1);  % To store area per face or per vertex (used for cluster-level & TFCE).
plm.Ykindstr = cell(Ni,1);  % string to save the files later
for i = 1:Ni,
    
    % Read an initial version
    fprintf('Reading input %d/%d: %s\n',i,Ni,opts.i{i});
    Ytmp = palm_miscread(opts.i{i},opts.useniiclass);
    
    % If this is 4D read with the NIFTI class, it needs a mask now
    if strcmp(Ytmp.readwith,'nifticlass') && ndims(Ytmp.data) == 4,
        if Nm == 0,
            % If a mask hasn't been supplied, make one
            tmpmsk = false(Ytmp.extra.dat.dim(1:3));
            for a = 1:Ytmp.extra.dat.dim(2), % y coords
                for b = 1:Ytmp.extra.dat.dim(3), % z coords
                    I = squeeze(Ytmp.extra.dat(:,a,b,:));
                    inan = any(isnan(I),2);
                    iinf = any(isinf(I),2);
                    icte = sum(diff(I,1,2).^2,2) == 0;
                    tmpmsk(:,a,b) = ~ (inan | iinf | icte);
                end
            end
            tmpmsk = tmpmsk(:)';
            plm.masks{i} = palm_maskstruct(tmpmsk,Ytmp.readwith,Ytmp.extra);
        else
            % If a mask was supplied, check its size
            if any(Ytmp.extra.dat.dim(1:3) ~= size(plm.masks{i}.data)),
                error([...
                    'The size of the data does not match the size of the mask:\n' ...
                    '- Data file %d (%s)\n' ...
                    '- Mask file %d (%s)'],i,opts.i{i},i,opts.m{i})
            end
        end
    end
    
    % Now deal with the data
    if ndims(Ytmp.data) == 2,
        
        % Transpose if that was chosen.
        if opts.transposedata,
            Ytmp.data = Ytmp.data';
        end
        
        % For the first input data, keep the size to
        % compare with the others, then check the size
        if i == 1,
            plm.N = size(Ytmp.data,1);
        end
        if size(Ytmp.data,1) ~= plm.N,
            error([
                'At least two of the input data files do not have\n' ...
                'compatible sizes:\n' ...
                '- File %d (%s) has %d observations\n'   ...
                '- File %d (%s) has %d observations'], ...
                1,opts.i{1},plm.N, ...
                i,opts.i{i},size(Ytmp.data,1));
        end
        
        % Not all later functions are defined for file_array class,
        % so convert to double
        if strcmp(Ytmp.readwith,'nifticlass'),
            Ytmp.data = double(Ytmp.data);
        end
        
        % This should cover the CSV files and DPX 4D files that
        % were converted to CSV with 'dpx2csv' and then transposed.
        plm.Yset{i} = Ytmp.data;
        
    elseif ndims(Ytmp.data) == 4,
        
        % For the first input data, keep the size to
        % compare with the others, then check the size
        if i == 1,
            plm.N = size(Ytmp.data,4);
        end
        if size(Ytmp.data,4) ~= plm.N,
            error([
                'At least two of the input data files do not have\n' ...
                'compatible sizes:\n' ...
                '- File %d (%s) has %d observations\n'   ...
                '- File %d (%s) has %d observations'], ...
                1,opts.i{1},plm.N, ...
                i,opts.i{i},size(Ytmp.data,4));
        end
        
        % Sort out loading for the NIFTI class
        if strcmp(Ytmp.readwith,'nifticlass'),
            tmpmsk = plm.masks{i}.data(:)';
            
            % Read each volume, reshape and apply the mask
            plm.Yset{i} = zeros(plm.N,sum(tmpmsk));
            for n = 1:plm.N,
                tmp = Ytmp.extra.dat(:,:,:,n);
                tmp = tmp(:)';
                plm.Yset{i}(n,:) = tmp(tmpmsk);
            end
        else
            % If not read with the NIFTI class, read all immediately
            plm.Yset{i} = palm_conv4to2(Ytmp.data);
        end
    end
    
    % Check if the size of data is compatible with size of mask.
    % If read with the NIFTI class, this was already taken care of
    % and can be skipped.
    if ~ strcmp(Ytmp.readwith,'nifticlass'),
        if Nm > 0 && size(plm.Yset{i},2) ~= numel(plm.masks{i}.data),
            error([...
                'The size of the data does not match the size of the mask:\n' ...
                '- Data file %d (%s)\n' ...
                '- Mask file %d (%s)'],i,opts.i{i},i,opts.m{i})
        end
    end
    
    % Make mask that removes constant values, Inf and NaN. This will be
    % merged with the user-supplied mask, if any, or will be the sole mask
    % available to select the datapoints of interest.
    if Nm == 0 && ndims(Ytmp.data) == 4 ...
            && strcmp(Ytmp.readwith,'nifticlass'),
        maskydat = true(1,size(plm.Yset{i},2));
    else
        ynan = any(isnan(plm.Yset{i}),1);
        yinf = any(isinf(plm.Yset{i}),1);
        ycte = sum(diff(plm.Yset{i},1,1).^2) == 0;
        maskydat = ~ (ynan | yinf | ycte);
    end
    
    % Now apply the mask created above and the one supplied by the user
    % for each modality. If no masks were supplied, create them, except
    % for the NIFTI class, which should have been created above
    if strcmp(Ytmp.readwith,'nifticlass'),
        plm.masks{i}.data(plm.masks{i}.data) = maskydat(:);
    else
        if Nm == 0,
            plm.masks{i} = palm_maskstruct(maskydat,Ytmp.readwith,Ytmp.extra);
        else
            maskydat = plm.masks{i}.data(:) & maskydat(:);
            plm.masks{i}.data = reshape(maskydat,size(plm.masks{i}.data));
        end
    end
    plm.Yset{i} = plm.Yset{i}(:,maskydat);
    
    % Prepare a string with a representative name for the kind of data,
    % i.e., voxel for volumetric data,
    switch Ytmp.readwith,
        case {'nifticlass','fs_load_nifti','fsl_read_avw',...
                'spm_spm_vol','nii_load_nii','fs_load_mgh'},
            plm.Yisvol(i)   = true;
            plm.Ykindstr{i} = 'vox';
        case 'fs_read_curv',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = 'dpv';
        case 'dpxread',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = 'dpx'; % this may be overriden below if a surface file is supplied
        otherwise
            plm.Ykindstr{i} = 'dat';
    end
    
    % If this is a DPX/curvature file, and if one of the spatial
    % statistics has been invoked, check if surfaces are available
    % and with compatible size, then compute the area (dpv or dpf)
    if plm.Yissrf(i) && opts.spatial,
        if Ns == 0,
            error([ ...
                'To use cluster extent, cluster mass, or TFCE with vertexwise or facewise data\n'...
                'it is necessary to provide the surface files (with the option -s).%s'],'');
        elseif Ns == 1,
            s = 1;
        else
            s = i;
        end
        if size(plm.srf{s}.data.vtx,1) == size(plm.Yset{i},2);
            plm.Yisvtx(i)   = true;
            plm.Yisfac(i)   = false;
            plm.Ykindstr{i} = 'vtx';
        elseif size(plm.srf{s}.data.fac,1) == size(plm.Yset{i},2);
            plm.Yisvtx(i)   = false;
            plm.Yisfac(i)   = true;
            plm.Ykindstr{i} = 'fac';
        end
        plm.Yarea{i} = palm_calcarea(plm.srf{s},plm.Yisvtx(i));
    end
end
plm.nY = numel(plm.Yset);
plm.nmasks = numel(plm.masks);

% Create an intersection mask if NPC or MV is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.NPC || opts.MV || opts.evperdat,
    if plm.nmasks > 1,
        
        % If there is one mask per modality, make an instersection mask.
        maskinter = true(size(plm.masks{1}.data));
        for m = 1:plm.nmasks,
            maskinter = maskinter & plm.masks{m}.data;
        end
        
        % Note that this line below uses Ytmp, which is from the previous loop.
        % This can be used here because with NPC all data has the same size.
        plm.maskinter = palm_maskstruct(maskinter(:)',Ytmp.readwith,Ytmp.extra);
        
        % Apply it to further subselect data points
        for y = 1:plm.nY,
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{y}.data));
        end
    else
        
        % If only one mask was given.
        plm.maskinter = plm.masks{1};
        for y = 1:plm.nY,
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{1}.data));
        end
    end
end

% Make sure that all data have the same size if NPC or MV are to be done
if opts.NPC || opts.MV || opts.evperdat,
    plm.Ysiz = zeros(plm.nY,1);
    siz1 = size(plm.Yset{1});
    for y = 1:plm.nY,
        plm.Ysiz(y) = size(plm.Yset{y},2);
        sizy = size(plm.Yset{y});
        if any(siz1 ~= sizy),
            error('The sizes of some of the imaging modalities don''t match');
        end
    end
end

% Take this opportunity to save the masks if the user requested.
if opts.savemask,
    for y = 1:plm.nmasks,
        M = plm.masks{y};
        if plm.nY == 1,
            M.filename = sprintf('%s_mask',opts.o);
        elseif plm.nmasks == 1,
            M.filename = sprintf('%s_mask_allmods',opts.o);
        else
            M.filename = sprintf('%s_mask_mod%d',opts.o,y);
        end
        M.data = double(M.data);
        palm_miscwrite(M);
    end
    if plm.nY > 1 && (opts.NPC || opts.MV || opts.evperdat),
        M          = plm.maskinter;
        M.filename = sprintf('%s_intersection_mask',opts.o);
        M.data     = double(M.data);
        palm_miscwrite(M);
    end
end

% If MV was selected, make sure that Y is full rank.
if opts.MV && ~ opts.noranktest,
    fprintf('Testing rank of the data (for MV tests). To skip, use -noranktest.\n')
    Y = cat(3,plm.Yset{:});
    Y = permute(Y,[1 3 2]);
    failed = false(1,size(Y,3));
    for v = 1:size(Y,3),
        if rank(Y(:,:,v)) ~= plm.nY;
            failed(v) = true;
        end
    end
    if any(failed),
        fname = sprintf('%s_%s_mv_illconditioned',...
            opts.o,plm.Ykindstr{1});
        palm_quicksave(double(failed),0,opts,plm,[],[],fname);
        error([
            'One or more datapoints have ill-conditioned data. This makes\n' ...
            'it impossible to run multivariate analyses as MANOVA/MANCOVA.\n' ...
            'Please, see these datapoints marked as 1 in the file:\n' ...
            '%s.*\n'],fname); %#ok
    end
end

% Applies an inverse-normal transformation to all
% modalities if the user requested
if opts.inormal,
    for y = 1:plm.nY,
        plm.Yset{y} = palm_inormal( ...
            plm.Yset{y},            ...
            opts.inormal_meth,      ...
            opts.inormal_quanti);
    end
end

% Make the adjustments for the EE and ISE options.
% - if the user gives nothing, its EE by default.
% - if the user gives ISE only, it's ISE only
% - if the user gives EE only, it's EE only
% - if the user gives both, it's both
if ~opts.EE && ~opts.ISE,
    opts.EE  = true;
end

% Read and assemble the design matrix.
if Nd > 0,
    fprintf('Reading design matrix and contrasts.\n');
    Mtmp = palm_miscread(opts.d{1});
    if opts.evperdat,
        if ndims(Mtmp.data) == 2,
            plm.M = Mtmp.data;
        elseif ndims(Mtmp.data) == 4,
            if strcmp(Mtmp.readwith,'nifticlass'),
                tmpmsk = plm.maskinter.data(:)';
                plm.M = zeros(plm.N,sum(tmpmsk));
                for n = 1:plm.N,
                    tmp = Mtmp.extra.dat(:,:,:,n);
                    tmp = tmp(:)';
                    plm.M(n,:) = tmp(tmpmsk);
                end
            else
                plm.M = palm_conv4to2(Mtmp.data);
                plm.M = plm.M(:,plm.maskinter.data(:));
            end
        end
    else
        plm.M = Mtmp.data;
    end
    if size(plm.M,1) ~= plm.N,
        error([
            'The number of rows in the design matrix does not\n' ...
            'match the number of observations in the data.\n' ...
            '- Rows in the matrix: %d\n' ...
            '- Observations in the data: %d'],size(plm.M,1),plm.N);
    end
else
    % If a design matrix has not been specified, use a single column of ones and
    % make sure that ISE only (not EE) is used.
    plm.M    = ones(plm.N,1);
    opts.EE  = false;
    opts.ISE = true;
end
if any(isnan(plm.M(:))) || any(isinf(plm.M(:))),
    error('The design matrix cannot contain NaN or Inf.');
end

% Read and organise the contrasts.
plm.Cset = cell(0);
if Nt || Nf,
    
    % There can't be more F than t contrast files
    if Nf > Nt,
        warning([...
            'More F-contrast files than t-contrast files were supplied.\n'...
            '         The last %d file(s) supplied with -f will be ignored.\n'...
            Nf-Nt]); %#ok
    end
    
    % Each t contrast is treated separately, even if many are
    % specified in a VEST or CSV file.
    tcon = cell(Nt,1);
    c = 1;
    for t = 1:Nt,
        tmp = palm_miscread(opts.t{t});
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            tcon{t} = tmp.data;
            for j = 1:size(tcon{t},1),
                plm.Cset{c} = tcon{t}(j,:)';
                c = c + 1;
            end
        else
            error('Invalid t contrast file: %s',opts.t{t});
        end
    end
    
    % Each F contrast assembles the t contrast from
    % the corresponding loaded t contrast VEST file.
    for f = 1:Nf,
        tmp = palm_miscread(opts.f{f});
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            for j = 1:size(tmp.data,1),
                plm.Cset{c} = tcon{f}(logical(tmp.data(j,:)),:)';
                c = c + 1;
            end
        else
            error('Invalid F contrast file: %s',opts.t{t});
        end
    end
    nC = numel(plm.Cset);
    
    % If only some contrasts are to be run
    if ~isempty(opts.conlist),
        conlist = palm_miscread(opts.conlist);
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            conlist = unique(conlist.data(:));
            if min(conlist) < 1 || max(conlist) > nC,
                error('Contrasts specified with ''-conlist'' are out of bounds.');
            end
            plm.Cset = plm.Cset(conlist);
        else
            error('Invalid file: %s',opts.conlist);
        end
    end
else
    % If no constrasts were at all specified:
    if size(plm.M,2) == 1 || opts.evperdat,
        
        % If there is only 1 regressor, test its effect both
        % positive and negative.
        % The statistic will be t or v, depending on the number of VGs.
        plm.Cset{1} = 1;
        plm.Cset{2} = -1;
        
    else
        % Otherwise, run an F-test over all regressors in the design matrix.
        % The statistic will be F or G, depending on the number of VGs.
        plm.Cset{1} = eye(size(plm.M,2));
    end
end
if opts.fonly,
    for c = numel(plm.Cset):-1:1,
        if rank(plm.Cset{c}) <= 1,
            plm.Cset(c) = [];
        end
    end
end
plm.nC = numel(plm.Cset);
for c = 1:plm.nC,
    if any(isnan(plm.Cset{c}(:))) || any(isinf(plm.Cset{c}(:))),
        error('The constrasts cannot contain NaN or Inf.');
    end
    if opts.evperdat && size(plm.Cset{c},1) ~= 1,
        error('The contrast file must contain a single element when using the option ''-evperdat''.');
    end
end

% Partition the model according to the contrasts and design matrix.
% The partitioning needs to be done now, because some regression methods
% may not be used if correction over contrasts is needed and the relevant
% regressors aren't compatible with synchronised permutations/sign-flips
if ~ opts.cmcx && ~ opts.evperdat,
    plm.Xset = cell(plm.nC,1);
    seqtmp = zeros(plm.N,plm.nC);
    for c = 1:plm.nC,
        plm.Xset{c} = palm_partition(plm.M,plm.Cset{c},opts.pmethod1);
        [~,~,seqtmp(:,c)] = unique(plm.Xset{c},'rows');
    end
    if opts.corrcon && any(sum(diff(seqtmp,1,2).^2,2) ~= 0) ...
            && ~ any(strcmpi(opts.rmethod,{'terBraak','Manly'})),
        warning([ ...
            'You chose to correct over contrasts, but with the contrasts\n' ...
            '         given, this is not possible using the %s method.\n' ...
            '         Using instead the %s method.\n' ...
            '         If, however, you really want to use %s, use the\n' ...
            '         option -igrepx, which will ignore repeated values in\n' ...
            '         the rows of X when permuting.'], ...
            opts.rmethod,opts.rfallback,opts.rmethod);
        opts.rmethod = opts.rfallback;
    end
elseif opts.evperdat
    % FIXME!!
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same large block. Also treat the legacy format of
% a single column for the EBs.
if isempty(opts.eb),
    plm.EB = [];
    if opts.within || opts.whole,
        warning([ ...
            'Options -within and/or -whole ignored, because no file defining\n' ...
            '         the exchangeability blocks was supplied (option -eb).\n' ...
            '         Performing free shuffling.']); %#ok
    end
else
    plm.EB = palm_miscread(opts.eb);
    plm.EB = plm.EB.data;
    if isvector(plm.EB),
        if opts.within && opts.whole, % within + whole block shuffling
            plm.EB = [+ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        elseif opts.whole             % whole-block shuffling
            plm.EB = [+ones(plm.N,1) -plm.EB(:) (1:plm.N)'];
        else                          % within-block shuffling (this is the default, and not meant to be changed)
            plm.EB = [-ones(plm.N,1) +plm.EB(:) (1:plm.N)'];
        end
    elseif opts.within || opts.whole,
        warning([ ...
            'Options -within and/or -whole ignored, as the file defining\n' ...
            '         the exchangeability blocks (option -eb) already defines\n' ...
            '         how the data should be shuffled.%s'],'');
    end
    plm.EB = palm_reindex(plm.EB,'fixleaves');
end

% Load/define the variance groups.
if opts.singlevg,
    % If single VG, it's all ones
    plm.VG = ones(plm.N,1);
elseif strcmpi(opts.vg,'auto'),
    if isempty(plm.EB),
        % If auto, but there are no exchangeability blocks, it's all ones too
        plm.VG = ones(plm.N,1);
    else
        % Generate an initial dependence tree, to be used to define variance groups.
        % The tree used for the permutations later require the design matrix, and
        % varies for each contrast -- all to be taken care of later.
        Ptree  = palm_tree(plm.EB,(1:plm.N)');
        plm.VG = palm_ptree2vg(Ptree);
    end
else
    % The automatic variance groups can be overriden if the user specified
    % a file with the custom definitions.
    plm.VG = palm_miscread(opts.vg);
    plm.VG = plm.VG.data;
end
[tmp,~,plm.VG] = unique(plm.VG);
plm.nVG = numel(tmp);

% MV can't yet be used if nVG>1, although NPC remains an option
if opts.MV && plm.nVG > 1,
    error('There are more than one variance group. MV cannot be used (but NPC can).');
end

% Remove observations marked to be ignored, i.e, those that
% are have their own colum in the design matrix and which
% aren't part of any contrast
if opts.removeignored,
    
    % Indices of the observations to keep
    idx = true(plm.N,1);
    lM  = logical(M);
    F   = find(sum(lM,1) == 1);
    for f = numel(F):-1:1,
        funused = false(plm.nC,1);
        for c = 1:plm.nC,
            funused(c) = all(plm.Cset{c}(F(f),:) == 0);
        end
        if all(funused),
            idx(lM(:,F(f))) = false;
        else
            F(f) = [];
        end
    end
    
    % Modify all data as needed
    for y = 1:plm.nY,
        plm.Yset{y} = plm.Yset{y}(idx,:);
    end
    if ~ isempty(plm.EB),
        plm.EB = plm.EB(idx,:);
    end
    for c = 1:plm.nC,
        plm.Cset{c}(F,:) = [];
    end
    plm.M      = plm.M(idx,:);
    plm.M(:,F) = [];
    plm.N      = sum(idx);
    plm.VG     = plm.VG(idx);
    [tmp,~,plm.VG] = unique(plm.VG(idx));
    plm.nVG    = numel(tmp);
end

% Remove the variance groups with just 1 observation?
if plm.nVG > 1 && ~ opts.removevgbysize && (opts.vgdemean || opts.ev4vg) && ...
        any(sum(bsxfun(@eq,plm.VG,unique(plm.VG)'),1) == 1),
    warning([...
        'The options ''-vgdemean'' and ''-ev4vg'' require that observations\n' ...
        '         in variance groups of size 1 are removed.\n' ...
        '         Enabling the option ''-removevgsize1''%s.'],'');
    opts.removevgbysize = true;
end
if ~ opts.removevgbysize,
    tmpwarned = false;
    for u = 1:plm.nVG,
        if sum((plm.VG == u),1) == 1,
            if ~ tmpwarned,
                warning([...
                    'There are variance groups with just one observation.\n' ...
                    '         Consider using the option ''-removevgsize1'' to improve the\n' ...
                    '         variance estimates (at the cost of reducing sample size).%s'],'');
                tmpwarned = true;
            end
        end
    end
end
if opts.removevgbysize > 0,
    
    % Indices of the observations to keep
    uVG = unique(plm.VG)';
    idxvg = sum(bsxfun(@eq,plm.VG,uVG),1) <= opts.removevgbysize;
    idx   = any(bsxfun(@eq,plm.VG,uVG(~idxvg)),2);
    
    % Modify all data as needed
    for y = 1:plm.nY,
        plm.Yset{y} = plm.Yset{y}(idx,:);
    end
    if ~ isempty(plm.EB),
        plm.EB = plm.EB(idx,:);
    end
    plm.M = plm.M(idx,:);
    plm.N = sum(idx);
    [tmp,~,plm.VG] = unique(plm.VG(idx));
    plm.nVG = numel(tmp);
end

% Add one regressor for each variance group, if requested
if opts.ev4vg,
    Mvg = zeros(plm.N,plm.nVG);
    V = unique(plm.VG);
    for v = 1:plm.nVG,
        Mvg(plm.VG == V(v),v) = 1;
    end
    rM   = round(sum(diag(plm.M*pinv(plm.M))));
    Mnew = horzcat(plm.M,Mvg);
    if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG),
        plm.M = Mnew;
        nadded = plm.nVG;
    else
        Mnew = Mnew(:,1:end-1);
        if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG - 1),
            plm.M = Mnew;
            nadded = plm.nVG - 1;
        else
            error([ ...
                'It was not possible to add one regressor for each variance group\n' ...
                'perhaps because they already exist in the design. Check your design\n' ...
                'matrix and maybe consider including these regressors manually.%s'],'');
        end
    end
    for c = 1:plm.nC,
        plm.Cset{c} = vertcat(plm.Cset{c},...
            zeros(nadded,size(plm.Cset{c},2)));
    end
end

% Remove intercept from the design for the options -demean and -vgdemean
if opts.demean || opts.vgdemean,
    intercp = all(bsxfun(@eq,plm.M(1,:),plm.M),1);
    if any(intercp),
        for c = 1:plm.nC,
            if any(intercp*plm.Cset{c} ~= 0,2),
                error([ ...
                    'Contrast %d (and perhaps others) tests the intercept. This means\n' ...
                    'that the options ''-demean'' and ''-vgdemean'' cannot be used.\n' ...
                    'If ''-demean'' was added to calculate Pearson''s "r" or the "R^2"\n' ...
                    'note that these statistics cannot be computed for constant variables.%s'],c,'');
            else
                plm.Cset{c}(intercp,:) = [];
            end
        end
        plm.M(:,intercp) = [];
    end
end

% Mean center data and design
if opts.demean,
    
    % Demean design
    plm.M = bsxfun(@minus,plm.M,mean(plm.M,1));
    
    % Demean data
    for y = 1:plm.nY,
        plm.Yset{y} = bsxfun(@minus,...
            plm.Yset{y},mean(plm.Yset{y},1));
    end
end

% Mean center data and design, within VG
if opts.vgdemean,
    
    % For each VG
    V = unique(plm.VG);
    for v = 1:plm.nVG,
        vidx = plm.VG == V(v);
        
        % Demean design within VG
        plm.M(vidx,:) = bsxfun(@minus,plm.M(vidx,:),mean(plm.M(vidx),1));
        
        % Demean data within VG
        for y = 1:plm.nY,
            plm.Yset{y}(vidx,:) = bsxfun(@minus,...
                plm.Yset{y}(vidx,:),mean(plm.Yset{y}(vidx,:),1));
        end
    end
end
