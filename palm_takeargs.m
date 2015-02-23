function [opts,plm] = palm_takeargs(varargin)
% Handle the inputs for PALM.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
% http://brainder.org

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
        otmp = opts.o;
    else
        otmp = vararginx{idxa+1};
    end
    if ~strcmp(otmp(end),'_'),
        otmp = horzcat(otmp,'_');
    end
    cfgname = horzcat(otmp,'palmconfig.txt');
    [opth,~,~] = fileparts(cfgname);
    if ~isempty(opth) && ~exist(opth,'dir'),
        mkdir(opth);
    end
    palm_configrw(vararginx,cfgname);
end

% Number of input images/masks/surfaces
% These are NOT meant to be edited. Don't edit this!
Ni   = sum(strcmp(vararginx,'-i'));         % number of data inputs
Nm   = sum(strcmp(vararginx,'-m'));         % number of masks
Ns   = sum(strcmp(vararginx,'-s'));         % number of surfaces
Nd   = sum(strcmp(vararginx,'-d'));         % number of design files
Nt   = sum(strcmp(vararginx,'-t'));         % number of t-contrast files
Nf   = sum(strcmp(vararginx,'-f'));         % number of F-test files
Ncon = sum(strcmp(vararginx,'-con'));       % number of contrast files (t or F, mset format)
Nevd  = sum(strcmp(vararginx,'-evperdat')); % number of EV per datum inputs
opts.i    = cell(Ni,1);   % Input files (to constitute Y later)
opts.m    = cell(Nm,1);   % Mask file(s)
opts.s    = cell(Ns,1);   % Surface file(s)
opts.d    = cell(Nd,1);   % Design file(s)
opts.t    = cell(Nt,1);   % t contrast file(s)
opts.f    = cell(Nf,1);   % F contrast file(s)
opts.Ccon = cell(Ncon,1); % Contrast file(s) (t or F, mset format)
opts.Dcon = cell(Ncon,1); % Contrast file(s) (multivariate, mset format)
opts.eb       = [];       % File with definition of exchangeability blocks
opts.vg       = [];       % File with definition of variance groups
opts.EE       = false;    % To be filled below (don't edit this!)
opts.ISE      = false;    % To be filled below (don't edit this!)
opts.within   = false;    % To be filled below (don't edit this!)
opts.whole    = false;    % To be filled below (don't edit this!)
opts.conskipcount = 0;    % When saving the contrasts, skip how many from 1?
opts.singlevg = true;     % Make sure that sigle VG will be used if nothing is supplied (this is NOT a "default" setting, and it's not a setting at all, but hard coded. Don't edit it!)
opts.subjidx  = [];       % Filename of the indices of subjects to keep
plm.subjidx   = [];       % Indices of subjects to keep

% These are to be incremented below
a = 1; i = 1; m = 1; d = 1;
t = 1; f = 1; s = 1;
con = 1; ev = 1;

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
            
        case '-evperdat',
            
            % Use one EV per datum?
            opts.evperdat = true;
            opts.evdatfile{ev} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-'),
                opts.evpos{ev}(1,1) = 1; % EV position
                opts.evpos{ev}(1,2) = 1; % Design number
                opts.evpos{ev}(1,3) = ev; % File index (this isn't redundant)
                a = a + 2;
            elseif nargin == a + 2 || ...
                    ischar(vararginx{a+3}) && ...
                    strcmpi(vararginx{a+3}(1),'-'),
                opts.evpos{ev}(1,1) = eval(vararginx{a+2}); % EV position
                opts.evpos{ev}(1,2) = 1; % Design number
                opts.evpos{ev}(1,3) = ev; % File index (this isn't redundant)
                a = a + 3;
            else
                opts.evpos{ev}(1,1) = eval(vararginx{a+2}); % EV position
                opts.evpos{ev}(1,2) = eval(vararginx{a+3}); % Design number
                opts.evpos{ev}(1,3) = ev; % File index (this isn't redundant)
                a = a + 4;
            end
            ev = ev + 1;
            
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
            
        case '-con',
            
            % Get the contrast files from an .mset file or
            % pair of files. If a pair, the 1st is for Cset
            % and the second for Dset.
            opts.Ccon{con} = vararginx{a+1};
            if nargin == a + 1 || ...
                    ischar(vararginx{a+2}) && ...
                    strcmpi(vararginx{a+2}(1),'-'),
                opts.Dcon{con} = [];
                a = a + 2;
            else
                opts.Dcon{con} = vararginx{a+2};
                a = a + 3;
            end
            con = con + 1;

        case '-conskipcount',
            
            % Numbers to skip when saving the contrasts
            opts.conskipcount = vararginx{a+1};
            if ischar(opts.conskipcount),
                opts.conskipcount = str2double(opts.concountskip);
            end
            a = a + 2;
            
        case '-tonly',
            
            % Run only the t-contrasts
            opts.tonly = true;
            a = a + 1;
            
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
            
        case '-tfce1D',
            
            % Shortcut for -tfce_H 2 -tfce_E 2 -tfce_C 6,
            % i.e., parameters for TFCE in 2D mode
            opts.tfce.H      = 2;
            opts.tfce.E      = 2;
            opts.tfce.conn   = 6;
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
            
        case '-concordant',
            
            % For the NPC, favour alternatives with the same sign?
            opts.concordant = true;
            a = a + 1;
            
        case '-corrmod',
            
            % Correct over modalities.
            opts.corrmod = true;
            a = a + 1;
            
        case '-corrcon',
            
            % Correct over contrasts.
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
                    'Kennedy',         ... % Kennedy won't be in the help
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
            if nargin == a || (nargin > a && strcmp(vararginx{a+1}(1),'-')),
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
            
        case '-npcmod',
            
            % Correct over modalities.
            opts.NPC    = true;
            opts.npcmod = true;
            a = a + 1;
            
        case '-npccon',
            
            % Correct over contrasts.
            opts.NPC    = true;
            opts.npccon = true;
            a = a + 1;
   
        case '-mv',
            
            % Compute classic multivariate statistics
            if nargin == a,
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1}(1),'-'),
                opts.MV = true;
                a = a + 1;
                
            elseif nargin > a,
                
                % Which multivariate statistic to use?
                methlist = {            ...
                    'auto',             ...
                    'HotellingTsq',     ...
                    'Wilks',            ...
                    'Lawley',           ...
                    'Lawley-Hotelling', ...
                    'Pillai',           ...
                    'Roy',              ...
                    'Roy-ii',           ...
                    'Roy-iii',          ...
                    'CCA'};
                methidx = strcmpi(vararginx{a+1},methlist);
                
                % Check if method exists, and load extra parameters if needed
                if ~any(methidx);
                    error('Multivariate statistic "%s" unknown.',vararginx{a+1});
                elseif strcmpi(vararginx{a+1},'CCA'),
                    opts.MV  = false;
                    opts.CCA = true;
                    if nargin == a + 1 || ...
                            ischar(vararginx{a+2}) && ...
                            strcmpi(vararginx{a+2}(1),'-'),
                        opts.ccaparm = 1;
                        a = a + 2;
                    elseif ischar(vararginx{a+2}),
                        opts.ccaparm = eval(vararginx{a+2});
                        a = a + 3;
                    else
                        opts.ccaparm = vararginx{a+2};
                        a = a + 3;
                    end
                else
                    opts.MV = true;
                    a = a + 2;
                end
                opts.mvstat = methlist{methidx};
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
                    error('Parameter "%s" unknown for the "-inormal" option.',parms{p});
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
        
        case '-transposedata',
            
            % Transpose the data if it's 2D?
            opts.transposedata = true;
            a = a + 1;
            
        case '-inputmv',
            
            % Treat the (sole) input as multivariate, that is,
            % each column is a variable in a multivariate model,
            % as opposed to independent univariate tests.
            % Useful with non-imaging data.
            opts.inputmv = true;
            a = a + 1;
            
        case '-verbosefilenames',
            
            % Don't simplify filenames when otherwise it would be possible
            opts.verbosefilenames = true;
            a = a + 1;
            
        case '-syncperms',
            
            % Synchronize permutations regardless of other options?
            opts.syncperms = true;
            a = a + 1;
            
        case {'-designperinput','-erie'},
            
            % Use one design matrix for each input dataset? This enables
            % syncP regardless.
            opts.designperinput = true;
            opts.syncperms      = true;
            a = a + 1;
            
        case '-subjidx',
                        
            % Indices of the subjects to keep in the design
            opts.subjidx = vararginx{a+1};
            a = a + 2;
            
        case '-quiet',
            
            % Don't print lines indicating progress
            opts.showprogress = false;
            a = a + 1;
            
        case '-nounivariate',
            
            % Save or not univariate tests? Useful with MV/NPC/CCA
            opts.saveunivariate = false;
            a = a + 1;
            
        case '-savedof',
            
            % Save a file with the degrees of freedom?
            opts.savedof = true;
            a = a + 1;

        case '-pmethodp',
            
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
                opts.pmethodp = methlist{methidx};
            else
                error([...
                    'The option "-pmethodp" requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        case '-pmethodr',
            
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
                opts.pmethodr = methlist{methidx};
            else
                error([...
                    'The option "-pmethodr" requires a method to be specified.\n' ...
                    'Consult the documentation.']);
            end
        otherwise
            error('Unknown option: "%s"',vararginx{a});
    end
end

% Check for the existance of other programs for input/output
palm_checkprogs;

% Make sure the NPC options make sense
% - if -npc only, it's -npcmod
% - if -npcmod only, it's -npcmod
% - if -npccon only, it's -npccon
% - if -npcmod and -npccon, it's both
if opts.NPC && ~ opts.npcmod && ~ opts.npccon,
    opts.npcmod = true;
end

% A quick check for the case of 1 EV per column of Y.
if opts.evperdat,
    if any([...
            opts.cmcx
            opts.ev4vg
            opts.pearson]),
        error([...
            'The option "-evperdat" is incompatible with the options listed below:\n' ...
            '"-cmcx"\n' ...
            '"-ev4vg"\n' ...
            '"-pearson"\n' ...
            'None of these can be enabled with "-evperdat".%s'],'');
    end
    if strcmpi(opts.rmethod,'terBraak'),
        error('The option "-evperdat" cannot be used with the ter Braak method (not implemented)');
    end
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
% if opts.spatial_npc && strcmpi(opts.npcmethod,'nichols'),
%     error('The Nichols combination method doesn''t allow spatial statistics (cluster and/or TFCE).')
% end

% Some extra packages for Octave
if opts.spatial && palm_isoctave,
    pkg load image
end

% No FWER or NPC if using draft mode
if opts.draft,
    if opts.corrmod || opts.corrcon,
        error('The draft mode does not allow FWER-correction, only FDR.');
    end
    if opts.NPC,
        error('The draft mode does not allow NPC.');
    end
    if opts.clustere_npc.do || opts.clustere_npc.do || opts.tfce_npc.do,
        error('The draft mode does not allow spatial statistics (cluster or TFCE).'); 
    end
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
if (Nt || Nf) && Ncon,
    error('Cannot mix options "-t" or "-f" with "-con".');
end
if Nt || Nf,
    if Nt > Nd,
        error('More t-contrast files (%d) than valid design files (%d) were supplied.',Nt,Nd);
    end
    if Nt ~= 1 && Nt ~= Nd,
        error(['The number of supplied t-contrast files (option "-t") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
    if Nf > Nt,
        error('More F-contrast files (%d) than t-contrast files (%d) were supplied.',Nf,Nt);
    end
elseif Ncon,
    if Ncon > Nd,
        error('More contrast files (%d) than design files (%d) were supplied.',Nt,Nd);
    end
    if Ncon ~= 1 && Ncon ~= Nd,
        error(['The number of supplied contrast files (option "-con") must be 1 or\n',...
            'the same as the number of design files (option "-d") (%d).'],Nd);
    end
end
if opts.pearson && (opts.NPC || opts.MV),
    warning([ ...
        'It''s not possible to compute the Pearson''s r or R^2 together with NPC or\n', ...
        '         multivariate methods. Disabling the options "-npc" and "-mv".%s'],'');
    opts.NPC = false;
    opts.MV  = false;
end
if opts.pearson && ~ opts.demean,
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if opts.pearson && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'})),
    warning([ ...
        'To compute Pearson''s "r" or the "R^2", the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if opts.CCA && ~ opts.demean,
    warning([ ...
        'To perform CCA, the data and the design\n' ...
        '         must be mean centered. Adding option "-demean".%s'],'');
    opts.demean = true;
end
if opts.CCA && ~ any(strcmpi(opts.pmethodr,{'beckmann','ridgway'})),
    warning([ ...
        'To perform CCA, the design must be\n' ...
        '         partitioned using the Beckmann or Ridgway schemes.'...
        '         Adding the option "-pmethodr Beckmann".%s'],'');
    opts.pmethodr = 'beckmann';
end
if opts.demean && opts.vgdemean && ~opts.pearson && ~opts.CCA,
    warning([...
        'Cannot use the option "-demean" together with "-vgdemean"\n'...
        '         Ignoring the option "-vgdemean".%s'],'');
    opts.vgdemean = false;
end
if opts.ev4vg && opts.vgdemean,
    error('Cannot use the option "-ev4vg" together with "-vgdemean"');
end
if opts.MV && ~ opts.singlevg,
    error('Currently MV is only allowed with a single VG, that is, assuming homoscedasticity.');
end
if opts.designperinput && opts.MV,
    error('It''s not possible to use the option "-onedesignperinput" together with the option "-mv".');
end
if ~opts.saveunivariate && ~opts.MV && ~opts.NPC && ~opts.CCA,
    error('The option "-nounivariate" can only be used with "-mv", "-npc" or "-cca".');
end
if opts.MV && opts.CCA,
    error('Cannot do classical MANOVA at the same time as CCA.');
end
if Ni > 1 && opts.inputmv,
    error('Option "-inputmv" cannot be used with more than one "-i".')
end
if opts.inputmv,
    opts.saveunivariate = false;
end
if opts.corrcon,
    opts.zstat = true;
end
if opts.concordant && ~ opts.NPC,
    error('The option "-concordant" is for use with NPC only.');
end
if opts.concordant && opts.twotail,
    error(['Cannot use "-concordant" together with "-twotail" (inadmissible).\n'...
        'Use either of these, but not both together.%s'],'');
end
if opts.tonly && opts.fonly,
    error('Cannot use "-tonly" together with "-fonly".');
end

% Initialize the random number generator (if nP = 0, no need for that)
if opts.nP0,
    if palm_isoctave,
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
            opts.seed = 'reset';
        end
        rand('state',opts.seed); %#ok
    else
        if any(strcmpi(opts.seed,{'reset','shuffle','twist'})),
            opts.seed = 'shuffle';
        end
        rng('default');
        rng(opts.seed);
    end
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
        'There are more masks supplied with "-m" (%d masks) than\n'...
        'modalities supplied with "-i" (%d modalities)'],Nm,Ni);
elseif Nm > 1 && Nm ~= Ni,
    error([...
        'The number of masks supplied with "-m" (%d masks) is larger than 1,\n'...
        'but still not the same as the number of modalities supplied with\n'...
        'the option "-i" (%d modalities).'],Nm,Ni);
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

% Indices of the subjects that will be kept
if ~isempty(opts.subjidx),
    plm.subjidx = palm_miscread(opts.subjidx);
    plm.subjidx = round(plm.subjidx.data);
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
    
    % Now deal with the actual data
    if ndims(Ytmp.data) == 2,
        
        % Transpose if that was chosen.
        if opts.transposedata,
            Ytmp.data = Ytmp.data';
        end
        
        % Select subjects
        if ~isempty(plm.subjidx),
            Ytmp.data = Ytmp.data(plm.subjidx,:);
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
        
        % Select subjects
        if ~isempty(plm.subjidx),
            Ytmp.data = Ytmp.data(:,:,:,plm.subjidx);
        end
        
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
            % If not read with the NIFTI class, get all immediately
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
    % for this modality. If no masks were supplied, create them, except
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
                'spm_spm_vol','nii_load_nii'},
            plm.Yisvol(i)   = true;
            plm.Ykindstr{i} = '_vox';
        case 'fs_read_curv',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = '_dpv';
        case 'dpxread',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = '_dpx'; % this may be overriden below if a surface file is supplied
        case 'fs_load_mgh',
            if ndims(Ytmp.data) == 4 && ...
                    size(Ytmp.data,2) == 1 && ...
                    size(Ytmp.data,3) == 1,
                plm.Yissrf(i)   = true;
                plm.Ykindstr{i} = '_dpx'; % this may be overriden below if a surface file is supplied
            else
                plm.Yisvol(i)   = true;
                plm.Ykindstr{i} = '_vox';
            end
        otherwise
            plm.Ykindstr{i} = '_dat';
    end
    
    % If this is a DPX/curvature file, and if one of the spatial
    % statistics has been invoked, check if surfaces are available
    % and with compatible size, then compute the area (dpv or dpf)
    if plm.Yissrf(i) && opts.spatial,
        if Ns == 0,
            error([ ...
                'To use cluster extent, cluster mass, or TFCE with vertexwise or facewise data\n'...
                'it is necessary to provide the surface files (with the option "-s").%s'],'');
        elseif Ns == 1,
            s = 1;
        else
            s = i;
        end
        if size(plm.srf{s}.data.vtx,1) == ...
                max(size(plm.masks{i}.data));
            plm.Yisvtx(i)   = true;
            plm.Yisfac(i)   = false;
            plm.Ykindstr{i} = '_dpx';
        elseif size(plm.srf{s}.data.fac,1) == ...
                max(size(plm.masks{i}.data));
            plm.Yisvtx(i)   = false;
            plm.Yisfac(i)   = true;
            plm.Ykindstr{i} = '_dpf';
        else
            error([...
                    'Surface file does not match the input data:\n' ...
                    '- Surface file %d has %d vertices and %d faces (%s)\n' ...
                    '- Input data file %d has %d points (%s)'],...
                    s,size(plm.srf{s}.data.vtx,1),size(plm.srf{s}.data.fac,1),opts.s{s},...
                    i,max(size(plm.masks{i}.data)),opts.i{i});
        end
        plm.Yarea{i} = palm_calcarea(plm.srf{s}.data,plm.Yisvtx(i));
    end
end
plm.nY = numel(plm.Yset); % this is redefined below if opts.inputmv is set.

% Read and organise the EV per datum.
if opts.evperdat,
    plm.EVset = cell(Nevd,1);
    for ev = 1:Nevd,
        
        % Read an initial version
        fprintf('Reading EV for each datum %d/%d: %s\n',ev,Nevd,opts.i{ev});
        EVtmp = palm_miscread(opts.evdatfile{ev},opts.useniiclass);
        
        % If this is 4D read with the NIFTI class, it needs a mask now
        if strcmp(EVtmp.readwith,'nifticlass') && ndims(EVtmp.data) == 4,
            if Nm == 0,
                % If a mask hasn't been supplied, make one
                tmpmsk = false(EVtmp.extra.dat.dim(1:3));
                for a = 1:EVtmp.extra.dat.dim(2), % y coords
                    for b = 1:EVtmp.extra.dat.dim(3), % z coords
                        I = squeeze(EVtmp.extra.dat(:,a,b,:));
                        inan = any(isnan(I),2);
                        iinf = any(isinf(I),2);
                        icte = sum(diff(I,1,2).^2,2) == 0;
                        tmpmsk(:,a,b) = ~ (inan | iinf | icte);
                    end
                end
                tmpmsk = tmpmsk(:)';
                plm.masks{ev} = palm_maskstruct(tmpmsk,EVtmp.readwith,EVtmp.extra);
            else
                % If a mask was supplied, check its size
                if any(EVtmp.extra.dat.dim(1:3) ~= size(plm.masks{ev}.data)),
                    error([...
                        'The size of the data does not match the size of the mask:\n' ...
                        '- Data file %d (%s)\n' ...
                        '- Mask file %d (%s)'],ev,opts.evdatfile{ev},ev,opts.m{ev})
                end
            end
        end
        
        % Now deal with the actual data
        if ndims(EVtmp.data) == 2,
            
            % Transpose if that was chosen.
            if opts.transposedata,
                EVtmp.data = EVtmp.data';
            end
            
            % Select subjects
            if ~isempty(plm.subjidx),
                EVtmp.data = EVtmp.data(plm.subjidx,:);
            end
            
            % Compare size with the others
            if size(EVtmp.data,1) ~= plm.N,
                error([
                    'At least two of the input data files do not have\n' ...
                    'compatible sizes:\n' ...
                    '- File %d (%s) has %d observations\n'   ...
                    '- File %d (%s) has %d observations'], ...
                    1,opts.evdatfile{1},plm.N, ...
                    ev,opts.evdatfile{ev},size(EVtmp.data,1));
            end
            
            % Not all later functions are defined for file_array class,
            % so convert to double
            if strcmp(EVtmp.readwith,'nifticlass'),
                EVtmp.data = double(EVtmp.data);
            end
            
            % This should cover the CSV files and DPX 4D files that
            % were converted to CSV with 'dpx2csv' and then transposed.
            plm.EVset{ev} = EVtmp.data;
            
        elseif ndims(EVtmp.data) == 4,
            
            % Select subjects
            if ~isempty(plm.subjidx),
                EVtmp.data = EVtmp.data(:,:,:,plm.subjidx);
            end
            
            % Compare size with the others
            if size(EVtmp.data,4) ~= plm.N,
                error([
                    'At least two of the input data files do not have\n' ...
                    'compatible sizes:\n' ...
                    '- File %d (%s) has %d observations\n'   ...
                    '- File %d (%s) has %d observations'], ...
                    1,opts.evdatfile{1},plm.N, ...
                    ev,opts.evdatfile{ev},size(EVtmp.data,4));
            end
            
            % Sort out loading for the NIFTI class
            if strcmp(EVtmp.readwith,'nifticlass'),
                tmpmsk = plm.masks{ev}.data(:)';
                
                % Read each volume, reshape and apply the mask
                plm.EVset{ev} = zeros(plm.N,sum(tmpmsk));
                for n = 1:plm.N,
                    tmp = EVtmp.extra.dat(:,:,:,n);
                    tmp = tmp(:)';
                    plm.EVset{ev}(n,:) = tmp(tmpmsk);
                end
            else
                % If not read with the NIFTI class, get all immediately
                plm.EVset{ev} = palm_conv4to2(EVtmp.data);
            end
        end
        
        % Check if the size of data is compatible with size of mask.
        % If read with the NIFTI class, this was already taken care of
        % and can be skipped.
        if ~ strcmp(EVtmp.readwith,'nifticlass'),
            if Nm > 0 && size(plm.EVset{ev},2) ~= numel(plm.masks{ev}.data),
                error([...
                    'The size of the data does not match the size of the mask:\n' ...
                    '- Data file %d (%s)\n' ...
                    '- Mask file %d (%s)'],ev,opts.evdatfile{ev},ev,opts.m{ev})
            end
        end
        
        % Make mask that removes constant values, Inf and NaN. This will be
        % merged with the user-supplied mask, if any, or will be the sole mask
        % available to select the datapoints of interest.
        if Nm == 0 && ndims(EVtmp.data) == 4 ...
                && strcmp(EVtmp.readwith,'nifticlass'),
            maskydat = true(1,size(plm.EVset{ev},2));
        else
            ynan = any(isnan(plm.EVset{ev}),1);
            yinf = any(isinf(plm.EVset{ev}),1);
            ycte = sum(diff(plm.EVset{ev},1,1).^2) == 0;
            maskydat = ~ (ynan | yinf | ycte);
        end
        
        % Now apply the mask created above and the one supplied by the user
        % for this modality. If no masks were supplied, create them, except
        % for the NIFTI class, which should have been created above
        if strcmp(EVtmp.readwith,'nifticlass'),
            plm.masks{Ni+ev}.data(plm.masks{Ni+ev}.data) = maskydat(:);
        else
            if Nm == 0,
                plm.masks{Ni+ev} = palm_maskstruct(maskydat,EVtmp.readwith,EVtmp.extra);
            else
                maskydat = plm.masks{Ni+ev}.data(:) & maskydat(:);
                plm.masks{Ni+ev}.data = reshape(maskydat,size(plm.masks{Ni+ev}.data));
            end
        end
        plm.EVset{ev} = plm.EVset{ev}(:,maskydat);
    end
    plm.nEVdat = numel(plm.EVset);
    
    % Sizes of EV per datum arrays
    plm.EVsiz = zeros(plm.nEVdat,1);
    for ev = 1:plm.nEVdat,
        plm.EVsiz(ev) = size(plm.EVset{ev},2);
    end
end
plm.nmasks = numel(plm.masks);

% Create an intersection mask if NPC or MV is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.npcmod || opts.MV || opts.evperdat,
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
        if opts.evperdat,
            for ev = 1:plm.nEVdat,
                plm.EVset{ev} = plm.EVset{ev}(:,plm.maskinter.data(plm.masks{plm.nY+ev}.data));
            end
        end
    else
        
        % If only one mask was given.
        plm.maskinter = plm.masks{1};
        for y = 1:plm.nY,
            plm.Yset{y} = plm.Yset{y}(:,plm.maskinter.data(plm.masks{1}.data));
        end
        if opts.evperdat,
            for ev = 1:plm.nEVdat,
                plm.EVset{ev} = plm.EVset{ev}(:,plm.maskinter.data(plm.masks{1}.data));
            end
        end
    end
end
clear Ytmp EVtmp

% If the multiple columns of the (sole) input are to be treated
% in a multivariate fashion
if opts.inputmv,
    tmp1 = plm.Yset{1};
    tmp2 = plm.Ykindstr{1};
    nY  = size(tmp1,2);
    plm.Yset     = cell(nY,1);
    plm.Ykindstr = cell(nY,1);
    for y = 1:nY,
        plm.Yset{y}     = tmp1(:,y);
        plm.Ykindstr{y} = tmp2;
    end
    clear tmp1 tmp2 nY;
    plm.nY = numel(plm.Yset);
end

% A variable with the sizes of all modalities will be handy later
plm.Ysiz = zeros(plm.nY,1);
for y = 1:plm.nY,
    plm.Ysiz(y) = size(plm.Yset{y},2);
end
plm.Ycumsiz = vertcat(0,cumsum(plm.Ysiz));

% Make sure that all data have the same size if NPC or MV are to be done
if opts.npcmod || opts.MV || opts.evperdat,
    siz1 = size(plm.Yset{1});
    for y = 1:plm.nY,
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
        if plm.nY == 1 || plm.nmasks == 1,
            M.filename = sprintf('%s_mask',opts.o);
        else
            M.filename = sprintf('%s_mask_i%d',opts.o,y);
        end
        M.data = double(M.data);
        palm_miscwrite(M);
    end
    if plm.nY > 1 && (opts.npcmod || opts.MV || opts.evperdat),
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
        fname = sprintf('%s_mv_illconditioned',opts.o);
        palm_quicksave(double(failed),0,opts,plm,[],[],fname);
        error([
            'One or more datapoints have ill-conditioned data. It is\n' ...
            'not possible to run multivariate analyses as MANOVA/MANCOVA.\n' ...
            'Please, see these datapoints marked as 1 in the file:\n' ...
            '%s.*\n'],fname); %#ok
    end
end

% Applies an inverse-normal transformation to the modalities if the user requested
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

% Read and assemble the design matrices.
if opts.evperdat,
    evpos = cat(1,opts.evpos{:});
    if size(unique(evpos(:,1:2),'rows'),1) ~= size(evpos,1);
        error(['Some EV per datum have been defined for the same\n'...
            'position in the same design matrices.%s'],'');
    end
    desidx = unique(evpos(:,2));
    Ndev = desidx(end);
else
    Ndev = 0;
end
plm.Mset = cell(max(Nd,Ndev),1);
% Read & assemble each design matrix
if Nd > 0,
    fprintf('Reading design matrix and contrasts.\n');
    for m = 1:Nd,
        Mtmp = palm_miscread(opts.d{m});
        plm.Mset{m} = Mtmp.data;
        if size(plm.Mset{m},1) ~= plm.N,
            error([
                'The number of rows in the design matrix does\n' ...
                'not match the number of observations in the data.\n' ...
                '- Rows in the matrix: %d\n' ...
                '- Observations in the data: %d\n' ...
                'In file %s\n'], ...
                size(plm.Mset{m},1),plm.N,opts.d{m});
        end
        if any(isnan(plm.Mset{m}(:))) || any(isinf(plm.Mset{m}(:))),
            error([
                'The design matrix cannot contain NaN or Inf.\n' ...
                'In file %s\n'],opts.d{m});
        end
    end
end
plm.nM = numel(plm.Mset);
% Include the EV per datum
if Ndev > 0,
    for m = desidx',
        evposm = evpos(evpos(:,2) == m,:);
        if isempty(plm.Mset{m}),
            tmp = zeros(plm.N,size(evposm,1),plm.EVsiz(evposm(1,3)));
            for ev = 1:size(evposm,1),
                tmp(:,ev,:) = permute(plm.EVset{evposm(ev,3)},[1 3 2]);
            end
        else
            tmp = zeros(plm.N,size(evposm,1)+size(plm.Mset{m},2),plm.EVsiz(evposm(1,3)));
            idx = true(1,size(evposm,1)+size(plm.Mset{m},2));
            idx(evposm(:,1)) = false;
            tmp(:,idx,:) = repmat(plm.Mset{m},[1 1 plm.EVsiz(evposm(1,3))]);
            for ev = 1:size(evposm,1),
                tmp(:,find(~idx,ev),:) = permute(plm.EVset{evposm(ev,3)},[1 3 2]);
                plm.EVset{evposm(ev,3)} = [];
                plm.EVsiz(evposm(1,3)) = NaN;
            end
        end
        plm.Mset{m} = tmp;
    end
    plm = rmfield(plm,{'EVset','EVsiz'});
    plm.evperdat = false(plm.nM,1);
    plm.evperdat(desidx) = true;
end
% If no design was specified
if Nd == 0 && Ndev == 0,
    plm.Mset{1} = ones(plm.N,1);
    opts.EE  = false;
    opts.ISE = true;
    plm.nM = numel(plm.Mset);
end
% Some related sanity checks
if opts.designperinput && plm.nY ~= plm.nM,
    error(['To use the option "-onedesignperinput", the number of design files must\n' ...
        'match the number of inputs.\n%s'],'');
end
if opts.evperdat,
    for m = 1:plm.nM,
        if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
        for y = loopY,
            if size(plm.Yset{y},2) ~= size(plm.Mset{m},3),
                error(['The size of the data and the size of the EV per datum\n' ...
                    'don''t match.%s'],'')
            end
        end
    end
end

% Read and organise the contrasts for each design.
plm.Cset = cell(plm.nM,1);
plm.Dset = cell(plm.nM,1);
plm.nC   = zeros(plm.nM,1);
plm.nD   = zeros(plm.nM,1);
if Nt || Nf,
    % Load FSL style t and F contrasts

    % Load t and F contrast files
    tcon = cell(Nt,1);
    for t = 1:Nt,
        tmp = palm_miscread(opts.t{t});
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            tcon{t} = tmp.data;
        else
            error('Invalid t contrast file: %s',opts.t{t});
        end
    end
    fcon = cell(Nf,1);
    for f = 1:Nf,
        tmp = palm_miscread(opts.f{f});
        if any(strcmp(tmp.readwith,{'vestread','csvread','load'})),
            fcon{f} = tmp.data;
        else
            error('Invalid F contrast file: %s',opts.f{f});
        end
    end
    
    % For each valid design, assemble the contrasts.
    for m = 1:plm.nM,
        if Nt == 1;
            t = 1;
        else
            t = m;
        end
        c = 1;
        for j = 1:size(tcon{t},1),
            plm.Cset{m}{c} = tcon{t}(j,:)';
            c = c + 1;
        end
        if numel(fcon),
            for j = 1:size(fcon{t},1),
                plm.Cset{m}{c} = tcon{t}(logical(fcon{t}(j,:)),:)';
                c = c + 1;
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        for c = 1:plm.nC(m),
            plm.Dset{m}{c} = eye(plm.nY);
        end
        plm.nD(m) = numel(plm.Dset{m});
    end
    
elseif Ncon,
    % Load MSET style contrasts
    
    % Load all contrast pairs
    Ccon = cell(Ncon,1);
    Dcon = cell(Ncon,1);
    for con = 1:Ncon,
        tmp = palm_miscread(opts.Ccon{con});
        if strcmpi(tmp.readwith,'mset'),
            Ccon{con} = tmp.data;
        else
            error(['Files given to the option "-con" must be in .mset format.' ...
                'For .csv or .con files, use "-t"; for .fts files, use "-f".%s'],'');
        end
        if isempty(opts.Dcon{con}),
            for c = 1:numel(Ccon{con}),
                Dcon{con}{c} = eye(plm.nY);
            end
        else
            tmp = palm_miscread(opts.Dcon{con});
            if strcmpi(tmp.readwith,'mset'),
                Dcon{con} = tmp.data;
            else
                error(['Files given to the option "-con" must be in .mset format.' ...
                    'For .csv or .con files, use "-t"; for .fts files, use "-f".%s'],'');
            end
        end
    end
    
    % Assign the contrast sets to the design matrices
    for m = 1:plm.nM,
        if Ncon == 1;
            con = 1;
        else
            con = m;
        end
        plm.Cset{m} = Ccon{con};
        plm.Dset{m} = Dcon{con};
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
        for c = 1:plm.nC(m),
            plm.Cset{m}{c} = plm.Cset{m}{c}';
        end
        for d = 1:plm.nD(m),
            plm.Dset{m}{d} = plm.Dset{m}{d}';
        end
    end
else
    % If no constrasts were at all specified:
    for m = 1:plm.nM,
        if size(plm.Mset{m},2) == 1,
            
            % If there is only 1 regressor, test its effect both
            % positive and negative.
            % The statistic will be t or v, depending on the number of VGs.
            plm.Cset{m}{1} = 1;
            plm.Cset{m}{2} = -1;
            plm.Dset{m}{1} = eye(plm.nY);
            plm.Dset{m}{2} = eye(plm.nY);
        else
            % Otherwise, run an F-test over all regressors in the design matrix.
            % The statistic will be F or G, depending on the number of VGs.
            plm.Cset{m}{1} = eye(size(plm.Mset{m},2));
            plm.Dset{m}{1} = eye(plm.nY);
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% If only the t or F tests are to be performed
if opts.tonly,
    for m = 1:plm.nM,
        for c = plm.nC(m):-1:1,
            if rank(plm.Cset{m}{c}) > 1,
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
elseif opts.fonly,
    for m = 1:plm.nM,
        for c = plm.nC(m):-1:1,
            if rank(plm.Cset{m}{c}) <= 1,
                plm.Cset{m}(c) = [];
                plm.Dset{m}(c) = [];
            end
        end
        plm.nC(m) = numel(plm.Cset{m});
        plm.nD(m) = numel(plm.Dset{m});
    end
end

% Some more sanity checks
for m = 1:plm.nM,
    for c = 1:plm.nC(m),
        if any(isnan(plm.Cset{m}{c}(:))) || any(isinf(plm.Cset{m}{c}(:))),
            error('The constrasts cannot contain NaN or Inf.');
        end
        if size(plm.Cset{m}{c},1) ~= size(plm.Mset{m},2),
            error('The size of one or more contrasts don''t match the size of the respective design matrix.')
        end
    end
    for c = 1:plm.nD(m),
        if any(isnan(plm.Dset{m}{c}(:))) || any(isinf(plm.Dset{m}{c}(:))),
            error('The constrasts cannot contain NaN or Inf.');
        end
    end
end
if opts.MV
    if any(plm.nC ~= plm.nD),
        error('The number of contrasts on the L and R sides of the GLM must be the same');
    end
    for m = 1:plm.nM,
        for d = 1:plm.nD(m),
            if size(plm.Dset{m}{d},1) ~= plm.nY,
                error('The size of one or more MV contrasts don''t match the number of modalities.');
            end
        end
    end
end
for m = 1:plm.nM,
    for c = 1:plm.nC(m),
        if opts.concordant && rank(plm.Cset{m}{c}) > 1,
            error(['Cannot use the "-concordant" option with F-tests (inadmissible).\n'...
                'Use "-tonly" to run just the t-tests.%s'],'');
        end
    end
end


% Partition the model according to the contrasts and design matrix.
% The partitioning needs to be done now, because of the need for
% synchronised permutations/sign-flips
if opts.cmcx,
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM,
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m),
            plm.seq{m}{c} = (1:plm.N)';
        end
    end
else
    seqtmp = zeros(plm.N,sum(plm.nC));
    j = 1;
    plm.seq = cell(plm.nM,1);
    for m = 1:plm.nM,
        plm.seq{m} = cell(plm.nC(m),1);
        for c = 1:plm.nC(m),
            Xtmp = palm_partition(plm.Mset{m}(:,:,1),plm.Cset{m}{c},opts.pmethodp);
            [~,~,plm.seq{m}{c}] = unique(Xtmp,'rows');
            seqtmp(:,j) = plm.seq{m}{c};
            j = j + 1;
        end
    end
    tmp = sum(diff(seqtmp,1,2).^2,2);
    if (opts.corrcon || opts.npccon) && any(tmp(:) ~= 0),
        warning([ ...
            'You chose to correct over contrasts, or run NPC between contrasts, but with\n' ...
            '         the design(s) contrasts given, it is not possible to run synchronized permutations\n' ...
            '         without ignoring repeated elements in the design matrix (or matrices).\n' ...
            '         To solve this, adding the option -cmcx automatically.%s\n'],'');
        opts.cmcx = true;
    end
    if opts.corrcon || opts.npccon,
        opts.syncperms = true;
    end
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same large block. Also treat the legacy format of
% a single column for the EBs.
if isempty(opts.eb),
    plm.EB = [];
    if opts.within || opts.whole,
        error([ ...
            'Options -within and/or -whole require a file defining\n' ...
            '         the exchangeability blocks (option -eb).\n%s'],'');
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

% Remove the variance groups with tiny sample sizes?
if plm.nVG > 1 && ~ opts.removevgbysize && (opts.vgdemean || opts.ev4vg) && ...
        any(sum(bsxfun(@eq,plm.VG,unique(plm.VG)'),1) == 1),
    warning([...
        'The options "-vgdemean" and "-ev4vg" require that observations\n' ...
        '         in variance groups of size 1 are removed.\n' ...
        '         Enabling the option "-removevgbysize 1"%s.'],'');
    opts.removevgbysize = 1;
end
if ~ opts.removevgbysize,
    tmpwarned = false;
    for u = 1:plm.nVG,
        if sum((plm.VG == u),1) == 1,
            if ~ tmpwarned,
                warning([...
                    'There are variance groups with just one observation.\n' ...
                    '         Consider using the option "-removevgbysize 1" to improve the\n' ...
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
    for m = 1:plm.nM,
        plm.Mset{m} = plm.Mset{m}(idx,:);
    end
    plm.N = sum(idx);
    [tmp,~,plm.VG] = unique(plm.VG(idx));
    plm.nVG = numel(tmp);
end

% Add one regressor for each variance group if requested
if opts.ev4vg,
    for m = 1:plm.nM,
        Mvg = zeros(plm.N,plm.nVG);
        V = unique(plm.VG);
        for v = 1:plm.nVG,
            Mvg(plm.VG == V(v),v) = 1;
        end
        rM   = round(sum(diag(plm.Mset{m}*pinv(plm.Mset{m}))));
        Mnew = horzcat(plm.Mset{m},Mvg);
        if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG),
            plm.Mset{m} = Mnew;
            nadded      = plm.nVG;
        else
            Mnew = Mnew(:,1:end-1);
            if round(sum(diag(Mnew*pinv(Mnew)))) == (rM + plm.nVG - 1),
                plm.Mset{m} = Mnew;
                nadded      = plm.nVG - 1;
            else
                error([ ...
                    'It was not possible to add one regressor for each variance group\n' ...
                    'perhaps because they already exist in the design. Check your design\n' ...
                    'matrix and maybe consider including these regressors manually.%s'],'');
            end
        end
        for c = 1:plm.nC(m),
            plm.Cset{m}{c} = vertcat(plm.Cset{m}{c},...
                zeros(nadded,size(plm.Cset{m}{c},2)));
        end
    end
end

% Remove intercept from the design for the options -demean and -vgdemean
if opts.demean || opts.vgdemean,
    for m = 1:plm.nM,
        intercp = all(bsxfun(@eq,plm.Mset{m}(1,:),plm.Mset{m}),1);
        if any(intercp),
            for c = 1:plm.nC(m),
                if any(intercp*plm.Cset{m}{c}~=0,2),
                    error([ ...
                        'Contrast %d (and perhaps others) tests the intercept. This means\n' ...
                        'that the options "-demean" and "-vgdemean" cannot be used.\n' ...
                        'If "-demean" was added to calculate Pearson''s "r" or the "R^2"\n' ...
                        'note that these statistics cannot be computed for constant variables.%s'],c,'');
                else
                    plm.Cset{m}{c}(intercp,:) = [];
                end
            end
            plm.Mset{m}(:,intercp) = [];
        end
    end
end

% Mean center data and design (-demean)
if opts.demean,
    for m = 1:plm.nM,
        plm.Mset{m} = bsxfun(@minus,plm.Mset{m},mean(plm.Mset{m},1));
    end
    for y = 1:plm.nY,
        plm.Yset{y} = bsxfun(@minus,plm.Yset{y},mean(plm.Yset{y},1));
    end
end

% Mean center data and design, within VG
if opts.vgdemean,
    
    % For each VG
    V = unique(plm.VG);
    for v = 1:plm.nVG,
        vidx = plm.VG == V(v);
        
        % Demean design within VG
        for m = 1:plm.nM,
            plm.Mset{m}(vidx,:) = bsxfun(@minus,...
                plm.Mset{m}(vidx,:),mean(plm.Mset{m}(vidx,:),1));
        end
        
        % Demean data within VG
        for y = 1:plm.nY,
            plm.Yset{y}(vidx,:) = bsxfun(@minus,...
                plm.Yset{y}(vidx,:),mean(plm.Yset{y}(vidx,:),1));
        end
    end
end
