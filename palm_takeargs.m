function [opts,plm] = palm_takeargs(varargin)

% As varargin is actually from another function, fix it.
varargin = varargin{1};
nargin = numel(varargin);

% Check if the number of input images/lists
% match the number of masks.
Ni = sum(strcmp(varargin,'-i'));
Nm = sum(strcmp(varargin,'-m'));
Ns = sum(strcmp(varargin,'-s'));
Nt = sum(strcmp(varargin,'-t'));
Nf = sum(strcmp(varargin,'-f'));

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

% Define some defaults and organise all as a struct
opts.i               = cell(Ni,1);      % Input files (to constitute Y later)
opts.m               = cell(Nm,1);      % Mask file(s)
opts.s               = cell(Nm,1);      % Surface file(s)
opts.t               = cell(Nt,1);      % t contrast file(s)
opts.f               = cell(Nf,1);      % F contrast file(s)
opts.nP0             = 10000;           % Number of permutations
opts.eb              = [];              % File with definition of exchangeability blocks
opts.vg              = [];              % File with definition of variance groups
opts.PB              = false;           % Whole block permutation?
opts.lx              = true;            % Lexicographic permutations?
opts.EE              = true;            % Exchangeable errors?
opts.ISE             = false;           % Independent and symmetric errors?
opts.pmethod         = 'Beckmann';      % Method to partition the model.
opts.rmethod         = 'Freedman-Lane'; % Regression/permutation method.
opts.rfallback       = 'terBraak';      % Regression/permutation method if correcting over contrasts
opts.NPC             = false;           % Do non-parametric combination?
opts.cmethod         = 'Tippett';       % Combination method.
opts.cfallback       = 'Fisher';        % ...
opts.savepara        = false;           % Save parametric p-values too?
opts.corrmod         = false;           % FWER correction over modalities?
opts.corrcon         = false;           % FWER correction over contrasts?
opts.savemask        = false;           % Save the masks?
opts.FDR             = false;           % FDR adjustment?
opts.draft           = false;           % Run a draft scheme
opts.clustere_t.do   = false;           % Do cluster extent for the t-stat?
opts.clusterm_t.do   = false;           % Do cluster mass for the t-stat?
opts.clustere_F.do   = false;           % Do cluster extent for the F-stat?
opts.clusterm_F.do   = false;           % Do cluster mass for the F-stat?
opts.tfce.do         = false;           % Do TFCE?
opts.tfce.H          = 2;               % TFCE H parameter
opts.tfce.E          = 0.5;             % TFCE E parameter
opts.tfce.conn       = 6;               % TFCE connectivity neighbourhood
opts.clustere_npc.do = false;           % Do cluster extent for the NPC z-stat?
opts.clusterm_npc.do = false;           % Do cluster mass for the NPC z-stat?
opts.tfce_npc.do     = false;           % Do TFCE for NPC?
opts.tfce_npc.H      = 2;               % TFCE H parameter for NPC
opts.tfce_npc.E      = 0.5;             % TFCE E parameter for NPC
opts.tfce_npc.conn   = 6;               % TFCE connectivity neighbourhood for NPC

% These are to be incremented below
a = 1; i = 1; m = 1;
t = 1; f = 1; s = 1; 

% Take the input arguments
while a <= nargin,
    switch varargin{a},
        case '-i',
            
            % Get the filenames for the data.
            opts.i{i} = varargin{a+1};
            i = i + 1;
            a = a + 2;
            
        case '-m',
            
            % Get the filenames for the masks, if any.
            opts.m{m} = varargin{a+1};
            m = m + 1;
            a = a + 2;
            
        case '-s',
            
            % Get the filenames for the surfaces, if any.
            opts.s{s} = varargin{a+1};
            s = s + 1;
            a = a + 2;
            
        case '-d',
            
            % Get the design matrix file.
            opts.d = varargin{a+1};
            a = a + 2;
                        
        case '-t',
            
            % Get the t contrast files.
            opts.t{t} = varargin{a+1};
            t = t + 1;
            a = a + 2;
            
        case '-f',
            
            % Get the F contrast files.
            opts.f{f} = varargin{a+1};
            f = f + 1;
            a = a + 2;
            
        case '-eb',
            
            % Get the exchangeability blocks file.
            opts.eb = varargin{a+1};
            a = a + 2;
            
        case '-vg',
            
            % Get the variance groups file.
            opts.vg = varargin{a+1};
            a = a + 2;
            
        case '-o',
            
            % Output prefix for the files to be saved.
            opts.o = varargin{a+1};
            a = a + 2;
            
        case '-n',
            
            % Number of permutations
            opts.nP0 = varargin{a+1};
            if ischar(opts.nP0),
                opts.nP0 = str2double(opts.nP0);
            end
            a = a + 2;
            
        case '-c',
            
            % Threshold for cluster extent, t-stat
            opts.clustere_t.do  = true;
            opts.clustere_t.thr = varargin{a+1};
            if ischar(opts.clustere_t.thr),
                opts.clustere_t.thr = str2double(opts.clustere_t.thr);
            end
            a = a + 2;
            
        case '-C',
            
            % Threshold for cluster mass, t-stat
            opts.clusterm_t.do  = true;
            opts.clusterm_t.thr = varargin{a+1};
            if ischar(opts.clusterm_t.thr),
                opts.clusterm_t.thr = str2double(opts.clusterm_t.thr);
            end
            a = a + 2;
            
        case '-F',
            
            % Threshold for cluster extent, F-stat
            opts.clustere_F.do  = true;
            opts.clustere_F.thr = varargin{a+1};
            if ischar(opts.clustere_F.thr),
                opts.clustere_F.thr = str2double(opts.clustere_F.thr);
            end
            a = a + 2;
            
        case '-S',
            
            % Threshold for cluster mass, F-stat
            opts.clusterm_F.do  = true;
            opts.clusterm_F.thr = varargin{a+1};
            if ischar(opts.clusterm_F.thr),
                opts.clusterm_F.thr = str2double(opts.clusterm_F.thr);
            end
            a = a + 2;
            
        case '-T',
            
            % Do TFCE?
            opts.tfce.do   = true;
            opts.tfce.H    = 2;
            opts.tfce.E    = 0.5;
            opts.tfce.conn = 6;
            a = a + 1;
            
        case '-T2',
            
            % Do TFCE in 2D mode?
            opts.tfce.do   = true;
            opts.tfce.H    = 2;
            opts.tfce.E    = 1;
            opts.tfce.conn = 26;
            a = a + 1;
            
        case '-tfce_H',
            
            % TFCE H parameter
            opts.tfce.H = varargin{a+1};
            if ischar(opts.tfce.H),
                opts.tfce.H = str2double(opts.tfce.H);
            end
            a = a + 2;
            
        case '-tfce_E',
            
            % TFCE E parameter
            opts.tfce.E = varargin{a+1};
            if ischar(opts.tfce.E),
                opts.tfce.E = str2double(opts.tfce.E);
            end
            a = a + 2;
            
        case '-tfce_C',
            
            % TFCE connectivity
            opts.tfce.conn = varargin{a+1};
            if ischar(opts.tfce.conn),
                opts.tfce.conn = str2double(opts.tfce.conn);
            end
            a = a + 2;
            
        case '-cnpc',
            
            % Threshold for cluster extent, NPC, z-stat
            opts.clustere_npc.do  = true;
            opts.clustere_npc.thr = varargin{a+1};
            if ischar(opts.clustere_npc.thr),
                opts.clustere_npc.thr = str2double(opts.clustere_npc.thr);
            end
            a = a + 2;
            
        case '-Cnpc',
            
            % Threshold for cluster mass, NPC, z-stat
            opts.clusterm_npc.do  = true;
            opts.clusterm_npc.thr = varargin{a+1};
            if ischar(opts.clusterm_npc.thr),
                opts.clusterm_npc.thr = str2double(opts.clusterm_npc.thr);
            end
            a = a + 2;
             
        case '-Tnpc',
            
            % Do TFCE?
            opts.tfce_npc.do   = true;
            opts.tfce_npc.H    = 2;
            opts.tfce_npc.E    = 0.5;
            opts.tfce_npc.conn = 6;
            a = a + 1;
            
        case '-T2npc',
            
            % Do TFCE in 2D mode?
            opts.tfce_npc.do   = true;
            opts.tfce_npc.H    = 2;
            opts.tfce_npc.E    = 1;
            opts.tfce_npc.conn = 26;
            a = a + 1;
            
        case '-pb',
            
            % Define whether should permute blocks as a whole (-pb) or not
            opts.PB = true;
            a = a + 1;
            
        case '-ise',
            
            % Independent and symmetric errors (ISE)?
            % If yes, this means sign-flippings.
            % Note that EE is the default. If ISE is specified alone, then EE
            % is disabled, unless it is also specified with -ee, in which case
            % both permutations and sign-flippings are performed.
            if strncmpi(varargin{a+1},'y',1),
                opts.ISE = true;
            elseif strncmpi(varargin{a+1},'n',1),
                opts.ISE = false;
            else
                error('Unrecognised option ''%s'' for %s.',varargin{a+1},varargin{a});
            end
            a = a + 2;
            
        case '-ee',
            
            % Exchangeable errors (EE)?
            % If yes, this means permutations.
            if strncmpi(varargin{a+1},'y',1),
                opts.EE = true;
            elseif strncmpi(varargin{a+1},'n',1),
                opts.EE = false;
            else
                error('Unrecognised option ''%s'' for %s.',varargin{a+1},varargin{a});
            end
            a = a + 2;
            
        case '-pmethod',
            
            % Which method to use for to partition the model?
            methlist = {    ...
                'Guttman',  ...
                'Beckmann', ...
                'Ridgway'};
            methidx = strcmpi(varargin{a+1},methlist);
            if ~any(methidx);
                error('Partition method "%s" unknown.',varargin{a+1});
            else
                a = a + 2;
            end
            opts.pmethod = methlist{methidx};
            
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
            
        case '-savemask',
            
            % If the user wants to have also the masks used for each
            % modality
            opts.savemask = true;
            a = a + 1;
            
        case '-rmethod',
            
            % Which method to use for the regression/permutation?
            methlist = {           ...
                'Draper-Stoneman', ...
                'Still-White',     ...
                'Freedman-Lane',   ...
                'terBraak',        ...
                'Kennedy',         ...
                'Manly',           ...
                'Huh-Jhun',        ...
                'Smith'};
            methidx = strcmpi(varargin{a+1},methlist);
            if ~any(methidx);
                error('Regression/Permutation method "%s" unknown.',varargin{a+1});
            else
                a = a + 2;
            end
            opts.rmethod = methlist{methidx};
            
        case '-npc',
            
            % Do the non-parametric combination?
            opts.NPC = true;
            a = a + 1;
            
        case '-cmethod',
            
            % Which combining function to use for the combination?
            opts.NPC = true;
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
            methidx = strcmpi(varargin{a+1},methlist);
            
            % Check if method exists, and load extra parameters if needed
            if ~any(methidx);
                error('Combining method "%s" unknown.',varargin{a+1});
            elseif any(strcmpi(varargin{a+1},{...
                    'wilkinson',          ...
                    'darlington-hayes',   ...
                    'zaykin',             ...
                    'dudbridge-koeleman', ...
                    'jiang'})),
                if ischar(varargin{a+2}),
                    plm.npcparm = eval(varargin{a+2});
                else
                    plm.npcparm = varargin{a+2};
                end
                a = a + 3;
            elseif any(strcmpi(varargin{a+1},'dudbridge-koeleman2')),
                if ischar(varargin{a+2}),
                    plm.npcparm = eval(varargin{a+2});
                else
                    plm.npcparm = varargin{a+2};
                end
                if ischar(varargin{a+3}),
                    plm.npcparm2 = eval(varargin{a+3});
                else
                    plm.npcparm2 = varargin{a+3};
                end
                a = a + 4;
            else
                a = a + 2;
            end
            opts.cmethod = methlist{methidx};
            
        case '-draft'
            
            % Do a draft scheme
            opts.draft = true;
            a = a + 1;
            
        case '-fdr'
            
            % Compute FDR-adjusted p-values
            opts.FDR = true;
            a = a + 1;
            
        otherwise
            error('Unknown option %s',varargin{a});
    end
end

% Some obvious sanity check.
% The fallback is to EE and not ISE because ISE can always be present
% if no design matrix is given (in which case it's all ones).
if ~opts.ISE && ~opts.EE
    warning([...
        'You chose not to use either EE (exchangeable errors) or ISE (independent and\n'...
        'symmetric errors), and this is not possible. Changing assumptions to EE.%s'],'');
end

% Make sure opts.NPC is marked as true if any NPC spatial statistic was
% selected.
if any([...
        opts.clustere_npc.do ...
        opts.clusterm_npc.do ...
        opts.tfce_npc.do]),
    opts.NPC = true;
end

% No FWER or NPC if using draft mode
if opts.draft,
    if opts.corrmod || opts.corrcon,
        warning('The draft mode does not allow FWER-correction.\n%s',''); %#ok
    end
    if opts.NPC,
        warning('The draft mode does not allow NPC.\n%s',''); %#ok
    end
    opts.corrmod         = false;
    opts.corrcon         = false;
    opts.NPC             = false;
    opts.clustere_npc.do = false;
    opts.clusterm_npc.do = false;
    opts.tfce_npc.do     = false;
end

% Some NPC methods don't have an analytical form for the parametric p-value
if any(strcmpi(opts.cmethod,{'darlington-hayes','jiang'})),
    plm.nonpcppara = true;
    if opts.savepara,
        warning([...
            'No parametric combination p-value will be saved for the\n', ...
            '         methods Darlington-Hayes or Jiang%s'],'');
    end
    if any([ ...
            opts.clustere_npc.do   ...
            opts.clusterm_npc.do   ...
            opts.tfce_npc.do]'),
        warning([ ...
            'No NPC cluster-level or TFCE statistic will be produced for the\n', ...
            '         methods Darlington-Hayes or Jiang%s'],'');
        opts.clustere_npc.do = false;
        opts.clusterm_npc.do = false;
        opts.tfce_npc.do     = false;
    end
else
    plm.nonpcppara = false;
end
    
% ==============================
% FIXME --- add a check for the minimally sufficient number of
% parameters, and for options that together make no sense. 
% ==============================

% Strings for the filenames
if opts.NPC,
    plm.margstr = '_marg';
    plm.npcstr  = '_npc';
else
    plm.margstr = '';
end

% Read and organise the surfaces. If no surfaces have been loaded, but the
% user wants cluster extent, cluster mass, or TFCE, an error will be
% printed later down in the code.
if any([ ...
        opts.clustere_t.do     ...
        opts.clustere_F.do     ...
        opts.clusterm_t.do     ...
        opts.clusterm_F.do     ...
        opts.clustere_npc.do   ...
        opts.clusterm_npc.do   ...
        opts.tfce.do           ...
        opts.tfce_npc.do]') && ...
        Ns > 0,
    plm.srf = cell(Ns,1);
    for s = 1:Ns,
        plm.srf{s} = palm_miscread(opts.s{s});
    end
end

% Read and organise the masks. If there are no masks specified, one for
% each modality will be created after each modality is loaded.
if Nm == 0,
    plm.masks = cell(Ni,1);
else
    plm.masks = cell(Nm,1);
end
for m = 1:Nm,
    plm.masks{m} = palm_miscread(opts.m{m});
    plm.masks{m}.data(isnan(plm.masks{m}.data)) = 0;
    plm.masks{m}.data(isinf(plm.masks{m}.data)) = 0;
    plm.masks{m}.data = logical(plm.masks{m}.data);
end

% Read and organise the data.
plm.Yset     = cell(Ni,1);  % Regressands (Y)
plm.Yisvol   = zeros(Ni,1); % Is Y a volume image?
plm.Yissrf   = zeros(Ni,1); % Is Y a surface-based image (DPX)?
plm.Yisvtx   = false(Ns,1); % Is vertexwise?
plm.Yisfac   = false(Ns,1); % is facewise? (this is currently dichotomous with Yisvtx, but later there may be edges/connectivity too)
pln.Yarea    = cell(Ns,1);  % To store area per face or per vertex (used for cluster-level & TFCE inferences).
plm.Ykindstr = cell(Ni,1);  % string to save the files later
for i = 1:Ni,
    Ytmp = palm_miscread(opts.i{i});
    
    if ndims(Ytmp.data) == 2,
        
        % For the first input data, keep the size to
        % compare with the others.
        if i == 1,
            N = size(Ytmp.data,1);
        end
        
        % This should cover the CSV files and DPX 4D files
        % that were converted to CSV with 'dpx2csv'.
        if size(Ytmp.data,1) == N,
            plm.Yset{i} = Ytmp.data;
        else
            error([
                'At least two of the input data files do not have\n' ...
                'compatible sizes:\n' ...
                '- File %d (%s) has %d observations\n'   ...
                '- File %d (%s) has %d observations'], ...
                1,opts.i{1},N, ...
                i,opts.i{i},size(Ytmp.data,1));
        end
        
    elseif ndims(Ytmp.data) == 4,
        
        % For the first input data, keep the size to
        % compare with the others.
        if i == 1,
            N = size(Ytmp.data,4);
        end
        
        % This should cover the other cases,
        % i.e., 4D NIFTI and 4D MGH/MGZ files.
        if size(Ytmp.data,4) == N,
            plm.Yset{i} = palm_conv4to2(Ytmp.data);
        else
            error([
                'At least two of the input data files do not have\n' ...
                'compatible sizes:\n' ...
                '- File %d (%s) has %d observations\n'   ...
                '- File %d (%s) has %d observations'], ...
                1,opts.i{1},N, ...
                i,opts.i{i},size(Ytmp.data,1));
        end
    end
    
    % Prepare a string with a representative name for the kind of data,
    % i.e., voxel for volumetric data, 
    switch Ytmp.readwith,
        case {'fs_load_nifti','fsl_read_avw',...
                'spm_spm_vol','nii_load_nii','fs_load_mgh'},
            plm.Yisvol(i)   = true;
            plm.Ykindstr{i} = 'vox';
        case 'fs_read_curv',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = 'dpv';
        case 'dpxread',
            plm.Yissrf(i)   = true;
            plm.Ykindstr{i} = 'dpx';
        otherwise
            plm.Ykindstr{i} = 'dat';
    end
    
    % If this is a DPX/curvature file, and if one of the spatial
    % statistics has been invoked, check if surfaces are available
    % and with compatible size, then compute the area (dpv or dpf)
    if plm.Yissrf(i) && ...
            any([ ...
            opts.clustere_t.do     ...
            opts.clustere_F.do     ...
            opts.clusterm_t.do     ...
            opts.clusterm_F.do     ...
            opts.clustere_npc.do   ...
            opts.clusterm_npc.do   ...
            opts.tfce.do           ...
            opts.tfce_npc.do]'),
        if Ns == 0,
            error([ ...
                'To use cluster extent, cluster mass, or TFCE with vertexwise or facewise data\n'...
                'it is necessary to provide the surface files (with the surface geometry).%s'],'');
        elseif Ns == 1,
            s = 1;
        else
            s = i;
        end
        if size(plm.srf{s}.data.vtx,1) == size(plm.Yset{i},2);
            plm.Yisvtx(i) = true;
            plm.Yisfac(i) = false;
        elseif size(plm.srf{s}.data.fac,1) == size(plm.Yset{i},2);
            plm.Yisvtx(i) = false;
            plm.Yisfac(i) = true;
        end
        pln.Yarea{i} = palm_calcarea(plm.srf{s},plm.Yisvtx(i));
    end
    
    % Check if size of data is compatible with size of mask
    if Nm == 1, m = 1; elseif Nm > 1, m = i; end
    if Nm ~= 0 && size(plm.Yset{i},2) ~= numel(plm.masks{m}.data),
        error([...
            'The size of the data does not match the size of the mask:\n' ...
            '- Data file %d (%s)\n' ...
            '- Mask file %d (%s)'],i,opts.i{i},m,opts.m{m})
    end
    
    % Make an initial mask removing constant values, Inf and NaN
    ynan = any(isnan(plm.Yset{i}),1);
    yinf = any(isinf(plm.Yset{i}),1);
    ycte = sum(diff(plm.Yset{i},1,1).^2) == 0;
    maskydat = ~ (ynan | yinf | ycte);
    
    % Now apply the mask created above and the one supplied by the user
    % for each modality. If no masks were supplied, create them.
    if Nm == 1,
        maskydat          = plm.masks{1}.data(:) & maskydat(:);
        plm.masks{1}.data = reshape(maskydat,size(plm.masks{1}.data));
    elseif Nm > 1,
        maskydat          = plm.masks{i}.data(:) & maskydat(:);
        plm.masks{i}.data = reshape(maskydat,size(plm.masks{i}.data));        
    elseif Nm == 0,
        plm.masks{i}      = palm_maskstruct(maskydat,Ytmp.readwith,Ytmp.extra);
    end
    plm.Yset{i} = plm.Yset{i}(:,maskydat);
end
plm.nY = numel(plm.Yset);
plm.nmasks = numel(plm.masks);

% Create an intersection mask if NPC is to be done, and further apply
% to the data that was previously masked above, as needed.
if opts.NPC && plm.nmasks > 1,
    
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
    
elseif opts.NPC,
    
    % If only one mask was given.
    plm.maskinter = plm.masks{1};
end

% Take this opportunity to save the masks if the user requested.
if opts.savemask,
    for y = 1:plm.nmasks,
        if plm.nY == 1,
            ystr = '';
        elseif plm.nmasks == 1,
            ystr = '_allmods';
        else
            ystr = '_mod%d';
        end
        M          = plm.masks{y};
        M.filename = sprintf(horzcat('%s_mask',ystr),opts.o,y);
        M.data     = double(M.data);
        palm_miscwrite(M);
    end
    if opts.NPC,
        M          = plm.maskinter;
        M.filename = sprintf('%s_npc_mask',opts.o);
        M.data     = double(M.data);
        palm_miscwrite(M);
    end
end

% Read the design matrix.
if isfield(opts,'d'),
    plm.M = palm_miscread(opts.d);
    plm.M = plm.M.data;
    if size(plm.M,1) ~= N,
        error([
            'The number of rows in the design matrix does not\n' ...
            'match the number of observations in the data.\n' ...
            '- Rows in the matrix: %d\n' ...
            '- Observations in the data: %d'],size(plm.M,1),N);
    end
else
    % If a design matrix is not specified, use a single column of ones and
    % make sure that ISE only (not EE) is used.
    plm.M    = ones(N,1);
    opts.EE  = false;
    opts.ISE = true;
end

% Read and organise the contrasts.
plm.Cset = cell(0);
if Nt || Nf,
    c = 1;
    
    % Each t contrast is treated separately, even if many are
    % specified in a VEST file. For CSV, it's one contrast per file.
    for t = 1:numel(opts.t),
        tmp = palm_miscread(opts.t{t});
        if strcmp(tmp.readwith,'csvread'),
            plm.Cset{c} = tmp.data;
            c = c + 1;
        elseif strcmp(tmp.readwith,'vestread'),
            for j = 1:size(tmp.data,1),
                plm.Cset{c} = tmp.data(j,:)';
                c = c + 1;
            end
            tcon = tmp.data;
        end
    end
    
    % The F contrasts, if from CSV, are also one per file. If read
    % from a VEST file, each F contrast assembles the t contrast from
    % the last loaded t contrast VEST file.
    for f = 1:numel(opts.f),
        tmp = palm_miscread(opts.f{f});
        if strcmp(tmp.readwith,'csvread'),
            plm.Cset{c} = tmp.data;
            c = c + 1;
        elseif strcmp(tmp.readwith,'vestread'),
            for j = 1:size(tmp.data,1),
                plm.Cset{c} = tcon(logical(tmp.data(j,:)),:)';
                c = c + 1;
            end
        end
    end
else
    % If no constrasts were at all specified, run an F-test over all
    % regressors in the design matrix. If there is only 1 regressor, this
    % will be a t-test. If there are exchangeability blocks, the statistic
    % then will be either a v^2 or v.
    plm.Cset{1} = eye(size(plm.M,2));
end
plm.nC = numel(plm.Cset);

% Give a warning and change to ter Braak or Manly if correcting over
% contrasts except in the particular case with nC == 2 and one contrast
% is the mirror of the other.
if opts.corrcon && ~any(strcmpi(opts.rmethod,{'terbraak','manly'})),
    if plm.nC ~= 2 || any(plm.Cset{1} ~= -plm.Cset{2}),
        warning([ ...
            'You chose to correct over contrasts, but with the contrasts\n' ...
            '         given, this is not possible using the %s method.\n' ...
            '         Using instead the %s method.'],opts.rmethod,opts.rfallback);
        opts.rmethod = opts.rfallback;
    end
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same big block.
% Also, take the opportunity to define the variance groups.
if isempty(opts.eb),
    plm.EB  = ones(N,1);
    plm.VG  = plm.EB;
else
    plm.EB = palm_miscread(opts.eb);
    plm.EB = plm.EB.data;
    if isvector(plm.EB),
        plm.EB = plm.EB(:);
    else
        error('The exchangeability blocks file must contain a vector%s','');
    end
    [plm.VG,plm.EB] = palm_eb2vg(plm.EB,opts.PB);
end

% Make sure that outputs have been set before moving on
if ~ isfield(opts,'o') || isempty(opts.o),
    error('You must supply an output string with the option -o.');
end
