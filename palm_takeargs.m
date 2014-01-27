function [opts,plm] = palm_takeargs(varargin)

% Load the defaults
opts = palm_defaults;

% As vararginx is actually from another function, fix it.
if nargin == 1,
    vararginx = palm_configrw(varargin{1});
else
    vararginx = varargin;
    idxa = find(strcmpi(vararginx,'-o'));
    if isempty(idxa),
        cfgname = horzcat(opts.o,'_palmconfig.txt');
    else
        cfgname = horzcat(vararginx{idxa+1},'_palmconfig.txt');
    end
    palm_configrw(vararginx,cfgname);
end

% Check if the number of input images/lists
% match the number of masks.
Ni = sum(strcmp(vararginx,'-i'));  % number of data inputs
Nm = sum(strcmp(vararginx,'-m'));  % number of masks
Ns = sum(strcmp(vararginx,'-s'));  % number of surfaces
Nt = sum(strcmp(vararginx,'-t'));  % number of t-contrast files
Nf = sum(strcmp(vararginx,'-f'));  % number of F-test files

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

opts.i  = cell(Ni,1);  % Input files (to constitute Y later)
opts.m  = cell(Nm,1);  % Mask file(s)
opts.s  = cell(Nm,1);  % Surface file(s)
opts.t  = cell(Nt,1);  % t contrast file(s)
opts.f  = cell(Nf,1);  % F contrast file(s)
opts.eb = [];          % File with definition of exchangeability blocks
opts.vg = [];          % File with definition of variance groups

% These are to be incremented below
a = 1; i = 1; m = 1;
t = 1; f = 1; s = 1; 

% Take the input arguments
while a <= nargin,
    switch vararginx{a},
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
            opts.d = vararginx{a+1};
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
            
        case '-eb',
            
            % Get the exchangeability blocks file.
            opts.eb = vararginx{a+1};
            a = a + 2;
            
        case '-vg',
            
            % Get the variance groups file.
            opts.vg = vararginx{a+1};
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
            
            % Threshold for cluster extent, t-stat
            opts.clustere_t.do  = true;
            opts.clustere_t.thr = vararginx{a+1};
            if ischar(opts.clustere_t.thr),
                opts.clustere_t.thr = str2double(opts.clustere_t.thr);
            end
            a = a + 2;
            
        case '-C',
            
            % Threshold for cluster mass, t-stat
            opts.clusterm_t.do  = true;
            opts.clusterm_t.thr = vararginx{a+1};
            if ischar(opts.clusterm_t.thr),
                opts.clusterm_t.thr = str2double(opts.clusterm_t.thr);
            end
            a = a + 2;
            
        case '-F',
            
            % Threshold for cluster extent, F-stat
            opts.clustere_F.do  = true;
            opts.clustere_F.thr = vararginx{a+1};
            if ischar(opts.clustere_F.thr),
                opts.clustere_F.thr = str2double(opts.clustere_F.thr);
            end
            a = a + 2;
            
        case '-S',
            
            % Threshold for cluster mass, F-stat
            opts.clusterm_F.do  = true;
            opts.clusterm_F.thr = vararginx{a+1};
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
            
        case '-cnpc',
            
            % Threshold for cluster extent, NPC, z-stat
            opts.clustere_npc.do  = true;
            opts.clustere_npc.thr = vararginx{a+1};
            if ischar(opts.clustere_npc.thr),
                opts.clustere_npc.thr = str2double(opts.clustere_npc.thr);
            end
            a = a + 2;
            
        case '-Cnpc',
            
            % Threshold for cluster mass, NPC, z-stat
            opts.clusterm_npc.do  = true;
            opts.clusterm_npc.thr = vararginx{a+1};
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
            
        case '-sb',
            
            % Define whether should permute blocks as a whole (-pb) or not
            opts.SB = true;
            a = a + 1;
            
        case '-ee',
            
            % Exchangeable errors (EE)?
            % If yes, this means permutations.
            opts.EE = true;
            a = a + 1;
            
        case '-ise',
            
            % Independent and symmetric errors (ISE)?
            % If yes, this means sign-flippings.            
            opts.ISE= true;
            a = a + 1;
            
        case '-cmc',
            
            % Define whether Conditional Monte Carlo should be used or not
            opts.CMC = true;
            a = a + 1;
            
        case '-igrepx',
            
            % Define whether repeated rows in X should be ignored or not
            % when defining the permutations
            opts.igrepx = true;
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
            
        case '-save1-p',
            
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
                error([....
                    'The option -rmethod requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end

            
        case '-npc',
            
            % Do the non-parametric combination?
            opts.NPC = true;
            if nargin == a,
                a = a + 1;
                
            elseif nargin > a && strcmp(vararginx{a+1},'-'),
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
                        'Wilkinson',          ...
                        'Darlington-Hayes',   ...
                        'Zaykin',             ...
                        'Dudbridge-Koeleman', ...
                        'Jiang'})),
                    if ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                    else
                        plm.npcparm = vararginx{a+2};
                    end
                    a = a + 3;
                elseif strcmpi(vararginx{a+1},'Friston'),
                    if ischar(vararginx{a+2}) && ...
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
                    if ischar(vararginx{a+2}),
                        plm.npcparm = eval(vararginx{a+2});
                    else
                        plm.npcparm = vararginx{a+2};
                    end
                    if ischar(vararginx{a+3}),
                        plm.npcparm2 = eval(vararginx{a+3});
                    else
                        plm.npcparm2 = vararginx{a+3};
                    end
                    a = a + 4;
                else
                    a = a + 2;
                end
                opts.cmethod = methlist{methidx};
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
            
        case '-inormal',
            
            % Inverse-normal transformation?
            opts.inormal = true;
            a = a + 1;
            
        case '-seed'
            
            % Seed for the random number generator
            opts.seed = vararginx{a+1};
            if ischar(opts.seed) && ~ strcmpi(opts.seed,'shuffle'),
                opts.seed = str2double(opts.seed);
            end
            a = a + 2;
            
        case '-shrink',
            
            % Remove from the analysis observations that are 0 in X and
            % belong to VGs not being considered for this contrast
            opts.shrink = true;
            a = a + 1;
            
        case '-zstat',
            
            % Convert the statistic for each test (not NPC) to a z-score
            opts.zstat = true;
            a = a + 1;
            
        case '-pmethod', % removed from the help
            
            % Which method to use for to partition the model?
            if nargin > a,
                methlist = {    ...
                    'Guttman',  ...
                    'Beckmann', ...
                    'Ridgway'};
                methidx = strcmpi(vararginx{a+1},methlist);
                if ~any(methidx);
                    error('Partition method "%s" unknown.',vararginx{a+1});
                else
                    a = a + 2;
                end
                opts.pmethod = methlist{methidx};
            else
                error([....
                    'The option -pmethod requires a method to be specified.\n'...
                    'Consult the documentation.']);
            end
            
        otherwise
            error('Unknown option: ''%s''',vararginx{a});
    end
end

% Some obvious sanity check.
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
end

% Some NPC methods don't have an analytical form for the parametric p-value
if any(strcmpi(opts.cmethod,{'Darlington-Hayes','Jiang'})),
    plm.nonpcppara = true;
    if opts.savepara,
        warning([...
            'No parametric combination p-value will be saved for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
    end
    if any([ ...
            opts.clustere_npc.do   ...
            opts.clusterm_npc.do   ...
            opts.tfce_npc.do]'),
        warning([ ...
            'No NPC cluster-level or TFCE statistic will be produced for the\n', ...
            '         Darlington-Hayes or Jiang methods%s'],'');
        opts.clustere_npc.do = false;
        opts.clusterm_npc.do = false;
        opts.tfce_npc.do     = false;
    end
else
    plm.nonpcppara = false;
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
    plm.masks{m} = palm_miscread(opts.m{m},opts.useniiclass);
    if strcmp(plm.masks{m}.readwith,'nifticlass'),
        plm.masks{m}.data = double(plm.masks{m}.data);
    end
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
plm.Yarea    = cell(Ns,1);  % To store area per face or per vertex (used for cluster-level & TFCE inferences).
plm.Ykindstr = cell(Ni,1);  % string to save the files later

for i = 1:Ni,
    
    % Read a temporary version
    fprintf('Reading input %d/%d: %s\n',i,Ni,opts.i{i});
    Ytmp = palm_miscread(opts.i{i},opts.useniiclass);
    
    % If this is 4D, read with NIFTI, needs a mask now
    if strcmp(Ytmp.readwith,'nifticlass') && ndims(Ytmp.data) == 4,
        if Nm == 0 
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
            % Otherwise, check the size of the one supplied
            if Nm == 1, m = 1; elseif Nm > 1, m = i; end
            if any(Ytmp.extra.dat.dim(1:3) ~= size(plm.masks{m}.data)),
                error([...
                    'The size of the data does not match the size of the mask:\n' ...
                    '- Data file %d (%s)\n' ...
                    '- Mask file %d (%s)'],i,opts.i{i},m,opts.m{m})
            end
        end
    end
    
    % Now deal with the data
    if ndims(Ytmp.data) == 2,
        
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
                i,opts.i{i},size(Ytmp.data,1));
        end
        
        % Sort out loading for the NIFTI class
        if strcmp(Ytmp.readwith,'nifticlass'),
            if Nm == 1, m = 1; elseif Nm ~= 1, m = i; end
            tmpmsk = plm.masks{m}.data(:)';
            
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
    % If read with the NIFTI class, this was already taken care of.
    if ~ strcmp(Ytmp.readwith,'nifticlass'),
        if Nm == 1, m = 1; elseif Nm > 1, m = i; end
        if Nm > 0 && size(plm.Yset{i},2) ~= numel(plm.masks{m}.data),
            error([...
                'The size of the data does not match the size of the mask:\n' ...
                '- Data file %d (%s)\n' ...
                '- Mask file %d (%s)'],i,opts.i{i},m,opts.m{m})
        end
    end
    
    % Make a mask removing constant values, Inf and NaN. This will be
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
        if Nm == 1, m = 1; else m = i; end
        plm.masks{m}.data(plm.masks{m}.data) = maskydat(:);
    else
        if Nm == 0,
            plm.masks{i} = palm_maskstruct(maskydat,Ytmp.readwith,Ytmp.extra);
        else
            if Nm == 1, m = 1; else m = i; end
            maskydat = plm.masks{m}.data(:) & maskydat(:);
            plm.masks{m}.data = reshape(maskydat,size(plm.masks{m}.data));
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
    if plm.Yissrf(i) && ...
            any([ ...
            opts.clustere_t.do   ...
            opts.clustere_F.do   ...
            opts.clusterm_t.do   ...
            opts.clusterm_F.do   ...
            opts.clustere_npc.do ...
            opts.clusterm_npc.do ...
            opts.tfce.do         ...
            opts.tfce_npc.do]'),
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
            plm.Yisvtx(i) = true;
            plm.Yisfac(i) = false;
            plm.Ykindstr{i} = 'vtx';
        elseif size(plm.srf{s}.data.fac,1) == size(plm.Yset{i},2);
            plm.Yisvtx(i) = false;
            plm.Yisfac(i) = true;
            plm.Ykindstr{i} = 'fac';
        end
        plm.Yarea{i} = palm_calcarea(plm.srf{s},plm.Yisvtx(i));
    end
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
    if opts.NPC,
        M          = plm.maskinter;
        M.filename = sprintf('%s_npc_mask',opts.o);
        M.data     = double(M.data);
        palm_miscwrite(M);
    end
end

% Applies an inverse-normal transformation to all
% modalities if the user requested
if opts.inormal,
    for y = 1:plm.nY,
        plm.Yset{y} = palm_inormal(plm.Yset{y},'Waerden');
    end
end

% Make the adjustments for the EE and ISE options.
% - if the user gives nothing, its EE by default.
% - if the user gives ISE only, it's ISE only
% - if the user gives EE only, it's EE only
% - if the user gives both, it's both
opts.EE  = [];
opts.ISE = [];
if isempty(opts.EE) && isempty(opts.ISE),
    opts.EE  = true;
    opts.ISE = false;
elseif opts.ISE && isempty(opts.EE),
    opts.EE  = false;
elseif opts.EE  && isempty(opts.ISE),
    opts.ISE = false;
end

% Read the design matrix.
if isfield(opts,'d'),
    plm.M = palm_miscread(opts.d);
    plm.M = plm.M.data;
    if size(plm.M,1) ~= plm.N,
        error([
            'The number of rows in the design matrix does not\n' ...
            'match the number of observations in the data.\n' ...
            '- Rows in the matrix: %d\n' ...
            '- Observations in the data: %d'],size(plm.M,1),plm.N);
    end
else
    % If a design matrix is not specified, use a single column of ones and
    % make sure that ISE only (not EE) is used.
    plm.M    = ones(plm.N,1);
    opts.EE  = false;
    opts.ISE = true;
end

% Read and organise the contrasts.
plm.Cset = cell(0);
if Nt || Nf,

    % There can't be more F than t contrast files
    if Nf > Nt,
        warning([...
            'You specified more F-contrast files than t-contrast files\n'...
            'The last %d file(s) supplied with -f will be ignored.\n'...
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

% Partition the model according to the contrasts and design matrix.
% The partitioning needs to be done now, because some regression methods
% may not be used if correction over contrasts is needed and the relevant
% regressors aren't compatible with synchronised permutations/sign-flips
if ~ opts.igrepx,
    plm.Xset = cell(plm.nC,1);
    seqtmp = zeros(plm.N,plm.nC);
    for c = 1:plm.nC,
        plm.Xset{c} = palm_partition(plm.M,plm.Cset{c},'Guttman');
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
end
if opts.corrcon && opts.shrink,
    warning([ ...
        'You chose to correct over contrasts and also to remove from\n' ...
        '         the design the observations that are not considered\n' ...
        '         for inference on each of them. This is impossible.\n' ...
        '         Switching off the option ''-shrink''']); %#ok<WNTAG>
    opts.shrink = false;
end

% Read the exchangeability blocks. If none is specified, all observations
% are assumed to be in the same big block. Also treat the legacy format of
% a single column for the EBs.
if isempty(opts.eb),
    %plm.EB = [ones(plm.N,1) (1:plm.N)'];
    plm.EB = [];
else
    plm.EB = palm_miscread(opts.eb);
    plm.EB = plm.EB.data;
    if isvector(plm.EB),
        if opts.SB,  % whole-block shuffling
            plm.EB = [+ones(plm.N,1) plm.EB(:)];
        else         % within-block shuffling
            plm.EB = [-ones(plm.N,1) plm.EB(:) (1:plm.N)'];
        end
    end
    plm.EB = palm_reindex(plm.EB,'fixleaves'); 
end

% Load/define the variance groups. They depend only on permutations, not
% on sign flips
if isempty(opts.vg),
    % Generate an initial dependence tree, to be used to define variance groups.
    % The tree used for the permutations later require the design matrix, and
    % varies for each contrast -- all to be taken care of later.
    if isempty(plm.EB),
        plm.VG = ones(plm.N,1);
    else
        Ptree = palm_tree(plm.EB,(1:plm.N)');
        plm.VG = palm_ptree2vg(Ptree);
    end
else
    % The automatic variance groups can be overriden if the user specified
    % a file with the custom definitions.
    plm.VG = palm_miscread(opts.vg);
end
plm.nVG = numel(unique(plm.VG));