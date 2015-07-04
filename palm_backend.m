function palm_backend(varargin)
% This is the core PALM function.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Take the arguments. Save a small log if needed.
%clear global plm opts % comment for debugging
% global plm opts; % uncomment for debugging
ticI = tic;
[opts,plm] = palm_takeargs(varargin{:});

% Variables to store stuff for later.
plm.X        = cell(plm.nM,1); % effective regressors
plm.Z        = cell(plm.nM,1); % nuisance regressors
plm.eCm      = cell(plm.nM,1); % effective contrast (for Mp)
plm.eCx      = cell(plm.nM,1); % effective contrast (for the effective regressors only)
plm.eC       = cell(plm.nM,1); % final effective contrast (depends on the method)
plm.Mp       = cell(plm.nM,1); % partitioned model, joined
plm.rmethod  = cell(plm.nM,1); % regression method (to allow the case of 'noz')
plm.nEV      = cell(plm.nM,1); % number of regressors
plm.Hm       = cell(plm.nM,1); % hat (projection) matrix
plm.Rm       = cell(plm.nM,1); % residual forming matrix
plm.dRm      = cell(plm.nM,1); % diagonal elements of the residual forming matrix
plm.Gname    = cell(plm.nM,1); % name of the statistic for each contrast
plm.rM       = cell(plm.nM,1); % rank of the design matrix
prepglm      = cell(plm.nM,1); % to store the function that prepares data for the regression
fastpiv      = cell(plm.nM,1); % to store the function that calculates the pivotal statistic
G            = cell(plm.nY,1); % to store G at each permutation
df2          = cell(plm.nY,1); % to store df2 at each permutation
plm.Gpperm   = cell(plm.nY,1); % counter, for the permutation p-value
plm.G        = cell(plm.nY,1); % for the unpermuted G (and to be saved)
plm.df2      = cell(plm.nY,1); % for the unpermuted df2 (and to be saved)
plm.Gmax     = cell(plm.nY,1); % to store the max statistic (Y collapses)
plm.nP       = cell(plm.nM,1); % number of permutations for each contrast

% Complete the lower levels (for the inner loops) for the variables above
for m = 1:plm.nM,
    plm.X                    {m}    = cell (plm.nC(m),1);
    plm.Z                    {m}    = cell (plm.nC(m),1);
    plm.eCm                  {m}    = cell (plm.nC(m),1);
    plm.eCx                  {m}    = cell (plm.nC(m),1);
    plm.eC                   {m}    = cell (plm.nC(m),1);
    plm.Mp                   {m}    = cell (plm.nC(m),1);
    plm.rmethod              {m}    = cell (plm.nC(m),1);
    plm.nEV                  {m}    = zeros(plm.nC(m),1);
    plm.Hm                   {m}    = cell (plm.nC(m),1);
    plm.Rm                   {m}    = cell (plm.nC(m),1);
    plm.dRm                  {m}    = cell (plm.nC(m),1);
    plm.Gname                {m}    = cell (plm.nC(m),1);
    plm.rM                   {m}    = zeros(plm.nC(m),1);
    prepglm                  {m}    = cell (plm.nC(m),1);
    fastpiv                  {m}    = cell (plm.nC(m),1);
    plm.nP                   {m}    = zeros(plm.nC(m),1);
end
for y = 1:plm.nY,
    G                     {y}       = cell(plm.nM,1);
    df2                   {y}       = cell(plm.nM,1);
    plm.Gpperm            {y}       = cell(plm.nM,1);
    plm.G                 {y}       = cell(plm.nM,1);
    plm.df2               {y}       = cell(plm.nM,1);
    plm.Gmax              {y}       = cell(plm.nM,1);
    if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
    for m = loopM,
        G                 {y}{m}    = cell(plm.nC(m),1);
        df2               {y}{m}    = cell(plm.nC(m),1);
        plm.Gpperm        {y}{m}    = cell(plm.nC(m),1);
        plm.G             {y}{m}    = cell(plm.nC(m),1);
        plm.df2           {y}{m}    = cell(plm.nC(m),1);
        plm.Gmax          {y}{m}    = cell(plm.nC(m),1);
        for c = 1:plm.nC(m),
            plm.G         {y}{m}{c} = zeros(1,plm.Ysiz(y));
            plm.Gpperm    {y}{m}{c} = zeros(1,plm.Ysiz(y));
        end
    end
end
if opts.approx.negbin,
    plm.Gppermp                     = plm.Gpperm; % number of perms done, for the negative binomial mode
end
if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
    plm.Gperms                      = plm.G;      % The 1xYsiz(1) will be replaced by nPxYsiz(1) later below
end
if opts.approx.lowrank,
    Gperms                          = G;          % To store the initial full permutations
    plm.effrank                     = plm.rC;     % To store the effective rank of the set of permutations
    plm.moments                     = G;
    plm.U                           = G;
end
if opts.savemetrics,
    plm.metr                        = plm.Gname;  % to store permutation metrics
end

% Spatial stats, univariate
if opts.cluster.uni.do,
    plm.Gclu                        = cell(plm.nY,1); % to store cluster statistic
    plm.Gclumax                     = cell(plm.nY,1); % for the max cluster statistic
    for y = 1:plm.nY,
        plm.Gclu          {y}       = cell(plm.nM,1);
        plm.Gclumax       {y}       = cell(plm.nM,1);
        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
        for m = loopM,
            plm.Gclu      {y}{m}    = cell(plm.nC(m),1);
            plm.Gclumax   {y}{m}    = cell(plm.nC(m),1);
        end
    end
end
if opts.tfce.uni.do,
    Gtfce                           = cell(plm.nY,1); % to store TFCE at each permutation
    plm.Gtfcepperm                  = cell(plm.nY,1); % counter, for the TFCE p-value
    plm.Gtfce                       = cell(plm.nY,1); % to store TFCE statistic
    plm.Gtfcemax                    = cell(plm.nY,1); % for the max TFCE statistic
    for y = 1:plm.nY,
        Gtfce             {y}       = cell(plm.nM,1);
        plm.Gtfcepperm    {y}       = cell(plm.nM,1);
        plm.Gtfce         {y}       = cell(plm.nM,1);
        plm.Gtfcemax      {y}       = cell(plm.nM,1);
        if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
        for m = loopM,
            Gtfce         {y}{m}    = cell(plm.nC(m),1);
            plm.Gtfcepperm{y}{m}    = cell(plm.nC(m),1);
            plm.Gtfce     {y}{m}    = cell(plm.nC(m),1);
            plm.Gtfcemax  {y}{m}    = cell(plm.nC(m),1);
        end
    end
    if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
        plm.Gtfceperms              = plm.Gtfce; % The 1xYsiz(1) will be replaced by nPxYsiz(1) later below
    end
end

% Variables for NPC
if opts.NPC,
    plm.npcstr = '_npc_';                      % default string for the filenames.
    if       opts.npcmod && ~ opts.npccon,
        Gnpc                        = cell(1);              % to store the G-stats ready for NPC
        df2npc                      = cell(1);              % to store the df2 ready for NPC
        Gnpc  {1}                   = zeros(plm.nY,plm.Ysiz(1));
        df2npc{1}                   = zeros(plm.nY,plm.Ysiz(1));
        T                           = cell(plm.nM,1);  % to store T at each permutation
        plm.Tpperm                  = cell(plm.nM,1);  % counter, for the combined p-value
        Tppara                      = cell(plm.nM,1);  % for the combined parametric p-value
        plm.Tmax                    = cell(plm.nM,1);  % to store the max combined statistic
    elseif  ~ opts.npcmod && opts.npccon,
        Gnpc                        = cell(plm.nY,1);
        df2npc                      = cell(plm.nY,1);
        if opts.designperinput,
            for y = 1:plm.nY,
                Gnpc  {y}           = zeros(plm.nC(y),plm.Ysiz(y));
                df2npc{y}           = zeros(plm.nC(y),plm.Ysiz(y));
            end
        else
            for y = 1:plm.nY,
                Gnpc  {y}           = zeros(sum(plm.nC),plm.Ysiz(y));
                df2npc{y}           = zeros(sum(plm.nC),plm.Ysiz(y));
            end
        end
        T                           = cell(plm.nY,1);
        plm.Tpperm                  = cell(plm.nY,1);
        Tppara                      = cell(plm.nY,1);
        plm.Tmax                    = cell(plm.nY,1);
    elseif   opts.npcmod && opts.npccon,
        Gnpc                        = cell(1);
        df2npc                      = cell(1);
        if opts.designperinput,
            Gnpc  {1}               = zeros(plm.nY*plm.nC(1),plm.Ysiz(1));
            df2npc{1}               = zeros(plm.nY*plm.nC(1),plm.Ysiz(1));
        else
            Gnpc  {1}               = zeros(plm.nY*sum(plm.nC),plm.Ysiz(1));
            df2npc{1}               = zeros(plm.nY*sum(plm.nC),plm.Ysiz(1));
        end
        T                           = cell(1);
        plm.Tpperm                  = cell(1);
        Tppara                      = cell(1);
        plm.Tmax                    = cell(1);
    end
    
    % Spatial stats, NPC
    if opts.cluster.npc.do,
        plm.Tclu                    = plm.Tmax;   % to store cluster NPC statistic
        plm.Tclumax                 = plm.Tmax;   % for the max cluster NPC
    end
    if opts.tfce.npc.do,
        Ttfce                       = plm.Tmax;   % to store TFCE at each permutation
        plm.Ttfcepperm              = plm.Tmax;   % counter, for the TFCE p-value
        plm.Ttfce                   = plm.Tmax;   % for the unpermuted TFCE
        plm.Ttfcemax                = plm.Tmax;   % to store the max TFCE statistic
        if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
            plm.Ttfceperms          = plm.Ttfce;
        end
    end
end

% Variables for MV
if opts.MV,
    plm.mvstr                       = '_mv';           % default string for the filenames.
    Q                               = cell(plm.nM,1);  % to store MV G at each permutation
    plm.Qname                       = cell(plm.nM,1);
    Qdf2                            = cell(plm.nM,1);  % to store MV df2 at each permutation
    plm.Qpperm                      = cell(plm.nM,1);  % counter, for the MV permutation p-value
    Qppara                          = cell(plm.nM,1);  % for the MV parametric p-value
    fastmv                          = cell(plm.nM,1);
    pparamv                         = cell(plm.nM,1);
    plm.Qmax                        = cell(plm.nM,1);  % to store the max multivariate statistic
    plm.mvrev                       = cell(plm.nM,1);  % reverse the direction in which a stat is significant?
    
    % Spatial stats, multivariate
    if opts.cluster.mv.do,
        plm.Qclu                    = cell(plm.nM,1);  % to store cluster MV statistic
        plm.Qclumax                 = cell(plm.nM,1);  % for the max cluster MV
    end
    if opts.tfce.mv.do,
        Qtfce                       = cell(plm.nM,1);  % to store TFCE at each permutation
        plm.Qtfcepperm              = cell(plm.nM,1);  % counter, for the TFCE p-value
        plm.Qtfce                   = cell(plm.nM,1);  % for the unpermuted TFCE
        plm.Qtfcemax                = cell(plm.nM,1);  % to store the max TFCE statistic
    end
    
    % Lower levels of these variables
    for m = 1:plm.nM,
        Q{m}                        = cell(plm.nC(m),1);
        plm.Qname{m}                = cell(plm.nC(m),1);
        Qdf2{m}                     = cell(plm.nC(m),1);
        plm.Qpperm{m}               = cell(plm.nC(m),1);
        Qppara{m}                   = cell(plm.nC(m),1);
        fastmv{m}                   = cell(plm.nC(m),1);
        pparamv{m}                  = cell(plm.nC(m),1);
        plm.Qmax{m}                 = cell(plm.nC(m),1);
        plm.mvrev{m}                = cell(plm.nC(m),1);
        if opts.cluster.mv.do,
            plm.Qclu{m}             = cell(plm.nC(m),1);
            plm.Qclumax{m}          = cell(plm.nC(m),1);
        end
        if opts.tfce.mv.do,
            Qtfce{m}                = cell(plm.nC(m),1);
            plm.Qtfcepperm{m}       = cell(plm.nC(m),1);
            plm.Qtfce{m}            = cell(plm.nC(m),1);
            plm.Qtfcemax{m}         = cell(plm.nC(m),1);
        end
        for c = 1:plm.nC(m),
            plm.Q{m}{c}             = zeros(1,plm.Ysiz(1));
            plm.Qpperm{m}{c}        = zeros(1,plm.Ysiz(1));
        end
    end
    if opts.approx.negbin,
        plm.Qppermp                 = plm.Qpperm; % number of perms done, for the negative binomial mode
    end
    if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
        plm.Qperms                  = plm.Q; % The 1xYsiz(1) will be replaced by nPxYsiz(1) later below
        if opts.tfce.mv.do,
            plm.Qtfceperms          = plm.Qtfce;
        end
        if opts.approx.lowrank,
            Qperms                  = Q;          % To store the initial full permutations
            plm.momentsq            = Q;
            plm.Uq                  = Q;
        end
    end
end

% Variables for CCA
if opts.CCA,
    plm.mvstr = ''; % default string for the filenames.
end

% Functions and strings for spatial statistics
switch opts.cluster.stat,
    case 'extent',
        clusterfunc = @palm_clustere;
        opts.cluster.str = '_clustere';
    case 'mass',
        clusterfunc = @palm_clusterm;
        opts.cluster.str = '_clusterm';
    case 'density',
        clusterfunc = @palm_clusterd;
        opts.cluster.str = '_clusterd';
    case 'tippett',
        clusterfunc = @palm_clustert;
        opts.cluster.str = '_clustert';
end
switch opts.tfce.stat,
    case 'tfce',
        tfcefunc = @palm_tfce;
        opts.tfce.str = '_tfce';
    case 'density',
        tfcefunc = @palm_tfde;
        opts.tfce.str = '_tfde';
end

% Inital strings to save the file names later.
plm.ystr = cell(plm.nY,1);      for y = 1:plm.nY,      plm.ystr{y} = ''; end
plm.mstr = cell(plm.nM,1);      for m = 1:plm.nM,      plm.mstr{m} = ''; end
plm.cstr = cell(max(plm.nC),1); for c = 1:max(plm.nC), plm.cstr{c} = ''; end

% Create the function handles for the NPC.
if opts.NPC,
    
    plm.Tname = lower(opts.npcmethod);
    switch plm.Tname,
        
        case 'tippett',
            if opts.concordant,
                fastnpc = @(G,df1,df2)min( ...
                    tippett( G,df1,df2),   ...
                    tippett(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)tippett(G,df1,df2);
            end
            pparanpc    = @(T,nG)tippettp(T,nG);
            plm.npcrev  = true;
            
        case 'fisher',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max( ...
                    fisher( G,df1,df2),    ...
                    fisher(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)fisher(G,df1,df2);
            end
            pparanpc    = @(T,nG)fisherp(T,nG);
            plm.npcrev  = false;
            
        case 'stouffer',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max( ...
                    stouffer( G,df1,df2),  ...
                    stouffer(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)stouffer(G,df1,df2);
            end
            pparanpc    = @(T,nG)stoufferp(T,nG);
            plm.npcrev  = false;
            
        case 'wilkinson',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(             ...
                    wilkinson( G,df1,df2,plm.npcparm), ...
                    wilkinson(-G,df1,df2,plm.npcparm));
            else
                fastnpc = @(G,df1,df2)wilkinson(G,df1,df2,plm.npcparm);
            end
            pparanpc    = @(T,nG)wilkinsonp(T,nG,plm.npcparm);
            plm.npcrev  = false;
            
        case 'winer',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max( ...
                    winer( G,df1,df2),     ...
                    winer(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)winer(G,df1,df2);
            end
            pparanpc    = @(T,nG)winerp(T,nG);
            plm.npcrev  = false;
            
        case 'edgington',
            if opts.concordant,
                fastnpc = @(G,df1,df2)min( ...
                    edgington( G,df1,df2), ...
                    edgington(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)edgington(G,df1,df2);
            end
            pparanpc    = @(T,nG)edgingtonp(T,nG);
            plm.npcrev  = true;
            
        case 'mudholkar-george',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(       ...
                    mudholkargeorge( G,df1,df2), ...
                    mudholkargeorge(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)mudholkargeorge(G,df1,df2);
            end
            pparanpc    = @(T,nG)mudholkargeorgep(T,nG);
            plm.npcrev  = false;
            
        case 'friston',
            if opts.concordant,
                fastnpc = @(G,df1,df2)min( ...
                    friston( G,df1,df2),   ...
                    friston(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)friston(G,df1,df2);
            end
            pparanpc    = @(T,nG)fristonp(T,nG,plm.npcparm);
            plm.npcrev  = true;
            
        case 'darlington-hayes',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(                   ...
                    darlingtonhayes( G,df1,df2,plm.npcparm), ...
                    darlingtonhayes(-G,df1,df2,plm.npcparm));
            else
                fastnpc = @(G,df1,df2)darlingtonhayes(G,df1,df2,plm.npcparm);
            end
            plm.npcrev  = false;
            
        case 'zaykin',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(          ...
                    zaykin( G,df1,df2,plm.npcparm), ...
                    zaykin(-G,df1,df2,plm.npcparm));
            else
                fastnpc = @(G,df1,df2)zaykin(G,df1,df2,plm.npcparm);
            end
            pparanpc    = @(T,nG)zaykinp(T,nG,plm.npcparm);
            plm.npcrev  = false;
            
        case 'dudbridge-koeleman',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(                     ...
                    dudbridgekoeleman( G,df1,df2,plm.npcparm), ...
                    dudbridgekoeleman(-G,df1,df2,plm.npcparm));
            else
                fastnpc = @(G,df1,df2)dudbridgekoeleman(G,df1,df2,plm.npcparm);
            end
            pparanpc    = @(T,nG)dudbridgekoelemanp(T,nG,plm.npcparm);
            plm.npcrev  = false;
            
        case 'dudbridge-koeleman2',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(                                   ...
                    dudbridgekoeleman2( G,df1,df2,plm.npcparm,plm.npcparm2), ...
                    dudbridgekoeleman2(-G,df1,df2,plm.npcparm,plm.npcparm2));
            else
                fastnpc = @(G,df1,df2)dudbridgekoeleman2(G,df1,df2,plm.npcparm,plm.npcparm2);
            end
            pparanpc    = @(T,nG)dudbridgekoeleman2p(T,nG,plm.npcparm,plm.npcparm2);
            plm.npcrev  = false;
            
        case 'taylor-tibshirani',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(        ...
                    taylortibshirani( G,df1,df2), ...
                    taylortibshirani(-G,df1,df2));
            else
                fastnpc = @(G,df1,df2)taylortibshirani(G,df1,df2);
            end
            pparanpc    = @(T,nG)taylortibshiranip(T,nG);
            plm.npcrev  = false;
            
        case 'jiang',
            if opts.concordant,
                fastnpc = @(G,df1,df2)max(         ...
                    jiang( G,df1,df2,plm.npcparm), ...
                    jiang(-G,df1,df2,plm.npcparm));
            else
                fastnpc = @(G,df1,df2)jiang(G,df1,df2,plm.npcparm);
            end
            plm.npcrev  = false;
    end
    
    % For the NPC methods in which the most significant stats are the
    % smallest, rather than the largest, use reverse comparisons.
    if plm.npcrev,
        npcrel  = @le;
        npcextr = @min;
    else
        npcrel  = @ge;
        npcextr = @max;
    end
end

% Create the permutation set, while taking care of the synchronized
% permutations (see the inner loop below)
if opts.syncperms,
    if ~ opts.approx.noperm,
        if isempty(plm.EB),
            if opts.savemetrics,
                [plm.Pset,plm.nP{1}(1),plm.metr{1}{1}] = ...
                    palm_shuffree(plm.seq{1}{1},opts.nP0, ...
                    opts.cmcp,opts.EE,opts.ISE,false);
            else
                [plm.Pset,plm.nP{1}(1)] = ...
                    palm_shuffree(plm.seq{1}{1},opts.nP0, ...
                    opts.cmcp,opts.EE,opts.ISE,false);
            end
        else
            if opts.savemetrics,
                [plm.Pset,plm.nP{1}(1),plm.metr{1}{1}] = ...
                    palm_shuftree(opts,plm,1,1);
            else
                [plm.Pset,plm.nP{1}(1)] = ...
                    palm_shuftree(opts,plm,1,1);
            end
        end
        fprintf('Building null distribution.\n');
    else
        fprintf('Doing the approximation without permutations.\n');
    end

    % This is for the negative binomial mode
    if opts.approx.negbin && ~ opts.saveunivariate,
        dothisY = false(plm.nY,1);
    else
        dothisY = true(plm.nY,1);
    end
    dotheMVorCCA = true;
    P_outer = 1:plm.nP{1}(1);
else
    P_outer = 1;
end

% To calculate progress
if opts.syncperms,
    ProgressNum = 0;
    if opts.designperinput,
        ProgressDen = sum(plm.nC) * plm.nP{1}(1);
    else
        ProgressDen = sum(plm.nC) * plm.nP{1}(1) * plm.nY;
    end
else
    ProgressCon = 0;
end
tocI = toc(ticI);
fprintf('Elapsed time parsing inputs: ~ %g s\n',tocI);

% For each permutation (outer loop)
ticP = tic;
for po = P_outer,
    
    % For each design matrix
    for m = 1:plm.nM,
        
        % String with the counter
        if po == 1 && (plm.nM > 1 || opts.verbosefilenames),
            plm.mstr{m} = sprintf('_d%d',m);
        end
        
        % For each contrast
        for c = 1:plm.nC(m),
            
            % This is for the negative binomial mode
            if ~ opts.syncperms,
                if opts.approx.negbin && ~ opts.saveunivariate,
                    dothisY = false(plm.nY,1);
                else
                    dothisY = true(plm.nY,1);
                end
                dotheMVorCCA = true;
            end
            
            if po == 1,
                
                % String with the counter
                if (plm.nC(m) > 1 || opts.verbosefilenames),
                    ctmp = c + opts.conskipcount;
                    plm.cstr{c} = sprintf('_c%d',ctmp);
                end
                
                if opts.evperdat,
                    
                    % Partition the model, now using the method chosen by the user
                    nT = size(plm.Mset{m},3);
                    plm.X{m}{c} = zeros(plm.N,rank(plm.Cset{m}{c}),nT);
                    plm.Z{m}{c} = zeros(plm.N,size(plm.Mset{m},2)-rank(plm.Cset{m}{c}),nT);
                    for t = 1:nT,
                        [plm.X{m}{c}(:,:,t),plm.Z{m}{c}(:,:,t),plm.eCm{m}{c},plm.eCx{m}{c}] = ...
                            palm_partition(plm.Mset{m}(:,:,t),plm.Cset{m}{c},opts.pmethodr);
                    end
                    plm.Mp{m}{c} = cat(2,plm.X{m}{c},plm.Z{m}{c}); % partitioned design matrix, joined
                    
                    % Some methods don't work well if Z is empty, and there is no point in
                    % using any of them all anyway.
                    if all(any(logical(plm.Cset{m}{c}),2)),
                        plm.rmethod{m}{c} = 'noz';
                    else
                        plm.rmethod{m}{c} = opts.rmethod;
                    end
                    
                    % As above, but for evperdat
                    plm.Hm{m}{c} = zeros(plm.N,plm.N,nT);
                    plm.Rm{m}{c} = zeros(plm.N,plm.N,nT);
                    for t = 1:nT,
                        plm.Hm {m}{c}(:,:,t) = plm.Mp{m}{c}(:,:,t)*pinv(plm.Mp{m}{c}(:,:,t));
                        plm.Rm {m}{c}(:,:,t) = eye(plm.N) - plm.Hm{m}{c}(:,:,t);
                        plm.dRm{m}{c}(:,t)   = diag(plm.Rm{m}{c}(:,:,t)); % this is used for the pivotal statistic
                    end
                    
                else
                    
                    % Partition the model, now using the method chosen by the user
                    [plm.X{m}{c},plm.Z{m}{c},plm.eCm{m}{c},plm.eCx{m}{c}] = ...
                        palm_partition(plm.Mset{m},plm.Cset{m}{c},opts.pmethodr);
                    plm.Mp{m}{c} = horzcat(plm.X{m}{c},plm.Z{m}{c}); % partitioned design matrix, joined
                    
                    % To avoid rank deficiency issues after partitioning, remove
                    % columns that are all equal to zero. This won't be
                    % done for the evperdat case because it's too slow and
                    % can make the designs too different if EVs are dropped
                    % from just some of the tests.
                    idx = all(plm.X{m}{c} == 0,1);
                    plm.X  {m}{c}(:,idx) = [];
                    plm.eCx{m}{c}(idx,:) = [];
                    idx = all(plm.Mp{m}{c} == 0,1);
                    plm.Mp {m}{c}(:,idx) = [];
                    plm.eCm{m}{c}(idx,:) = [];
                    
                    % Some methods don't work well if Z is empty, and there is no point in
                    % using any of them all anyway.
                    if isempty(plm.Z{m}{c}),
                        plm.rmethod{m}{c} = 'noz';
                    else
                        plm.rmethod{m}{c} = opts.rmethod;
                    end
                    
                    % Residual-forming matrix. This is used by the ter Braak method and
                    % also to compute some of the stats later. Note that, even though the
                    % residual-forming matrix changes at every permutation, the trace
                    % for each VG remains unchanged, hence it's not necessary to recompute
                    % it for every permutation.
                    plm.Hm {m}{c} = plm.Mp{m}{c}*pinv(plm.Mp{m}{c});
                    plm.Rm {m}{c} = eye(plm.N) - plm.Hm{m}{c};
                    plm.dRm{m}{c} = diag(plm.Rm{m}{c}); % this is used for the pivotal statistic
                    plm.rM {m}(c) = plm.N - round(sum(plm.dRm{m}{c})); % this is faster than rank(M)
                    
                end
                
                % Some other variables to be used in the function handles below.
                plm.nEV{m}(c) = size(plm.Mp{m}{c},2);    % number of regressors in the design
                plm.rM {m}(c) = plm.N - round(sum(plm.dRm{m}{c}(:,1,1))); % this is faster than rank(M)
                
                % MV/CCA
                %%% DOUBLE-CHECK THE DEGREES-OF-FREEDOM!!
                if opts.MV,
                    
                    % Make the 3D dataset
                    if opts.approx.negbin,
                        plm.Yq{m}{c} = cat(3,plm.Yset{:});
                    end
                    
                    % Define the functions for the stats. Note that none is
                    % available if nVG > 1, and this should have been
                    % checked when taking the arguments.
                    if plm.rC{m}(c) == 1 && any(strcmpi(opts.mvstat,{'auto','hotellingtsq'})),
                        plm.Qname{m}{c} = '_hotellingtsq';
                        pparamv  {m}{c} = @(Q)fasttsqp(Q,plm.N-plm.rM{m}(c),plm.nY);
                        plm.mvrev{m}{c} = false;
                    else
                        switch lower(opts.mvstat),
                            
                            case {'wilks','auto'},
                                plm.Qname{m}{c} = '_wilks';
                                plm.qfun        = @(H,E)wilks(H,E);
                                pparamv{m}{c}   = @(Q)wilksp(Q, ...
                                    plm.rC{m}(c),plm.N-plm.rM{m}(c),plm.nY);
                                plm.mvrev{m}{c} = true;
                                
                            case {'lawley','lawley-hotelling'},
                                plm.Qname{m}{c} = '_lawley-hotelling';
                                plm.qfun        = @(H,E)lawley(H,E);
                                pparamv{m}{c}   = @(Q)lawleyp(Q, ...
                                    plm.rC{m}(c),plm.N-plm.rM{m}(c),plm.nY);
                                plm.mvrev{m}{c} = false;
                                
                            case 'pillai',
                                plm.Qname{m}{c} = '_pillai';
                                plm.qfun        = @(H,E)pillai(H,E);
                                pparamv{m}{c}   = @(Q)pillaip(Q, ...
                                    plm.rC{m}(c),plm.N-plm.rM{m}(c),plm.nY);
                                plm.mvrev{m}{c} = false;
                                
                            case {'roy-ii','roy'},
                                plm.Qname{m}{c} = '_roy-ii';
                                plm.qfun        = @(H,E)roy_ii(H,E);
                                pparamv{m}{c}   = @(Q)roy_iip(Q, ...
                                    plm.rC{m}(c),plm.N-plm.rM{m}(c),plm.nY);
                                plm.mvrev{m}{c} = false;
                                
                            case 'roy-iii',
                                plm.Qname{m}{c} = '_roy-iii';
                                plm.qfun        = @(H,E)roy_iii(H,E);
                                plm.mvrev{m}{c} = false;
                        end
                    end
                    
                    % For the MV methods in which the most significant stats are the
                    % smallest, rather than the largest, use reverse comparisons.
                    if plm.mvrev{m}{c},
                        mvrel  = @le;
                        mvextr = @min;
                    else
                        mvrel  = @ge;
                        mvextr = @max;
                    end
                    
                elseif opts.CCA,
                    
                    % Output string
                    plm.Qname{m}{c} = sprintf('_cca%d',opts.ccaparm);
                    
                    % For CCA, this obviously false, but need to have when
                    % saving, do make the distributions correct
                    plm.mvrev{m}{c} = false;
                end
                
                % Effective rank of the matrix nP by Ysiz(y) used for the
                % low rank approximation. This number varies according to
                % the permutation and regression strategies, but it's
                % roughly as below:
                if opts.approx.lowrank,
                    plm.effrank{m}(c) = plm.rC{m}(c)*plm.N;
                end
                
                % Decide which method is going to be used for the regression and
                % permutations, compute some useful matrices for later and create
                % the appropriate function handle to prepare for the model fit.
                % Each of these small functions is a replacement for the generic
                % prototype function 'permglm.m', which is far slower.
                % Note that this swich needs to remain inside the for-loops over
                % designs and contrasts, because they vary. Nonetheless,
                % this all runs just for the 1st permutation.
                isterbraak = false;
                if opts.evperdat,
                    
                    % Pick the regression/permutation method
                    switch lower(plm.rmethod{m}{c}),
                        
                        case 'noz',
                            prepglm{m}{c}   = @(P,Y)noz3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCx{m}{c};
                            
                        case 'exact',
                            prepglm{m}{c}   = @(P,Y)exact3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCx{m}{c};
                            
                        case 'draper-stoneman',
                            prepglm{m}{c}   = @(P,Y)draperstoneman3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCm{m}{c};
                            
                        case 'still-white',
                            plm.Rz {m}{c}   = zeros(plm.N,plm.N,size(plm.Mset{m},3));
                            for t = 1:size(plm.Mset{m},3),
                                plm.Rz{m}{c}(:,:,t) = eye(plm.N) - plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                            end
                            prepglm{m}{c}   = @(P,Y)stillwhite3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCx{m}{c};
                            
                        case 'freedman-lane',
                            plm.Hz {m}{c}   = zeros(plm.N,plm.N,size(plm.Mset{m},3));
                            plm.Rz {m}{c}   = plm.Hz{m}{c};
                            for t = 1:size(plm.Mset{m},3),
                                plm.Hz{m}{c}(:,:,t) = plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                                plm.Rz{m}{c}(:,:,t) = eye(plm.N) - plm.Hz{m}{c}(:,:,t);
                            end
                            prepglm{m}{c}   = @(P,Y)freedmanlane3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCm{m}{c};
                            
                        case 'terbraak',
                            isterbraak      = true;
                            prepglm{m}{c}   = @(P,Y)terbraak3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCm{m}{c};
                            
                        case 'kennedy',
                            plm.Rz {m}{c}   = zeros(plm.N,plm.N,size(plm.Mset{m},3));
                            for t = 1:size(plm.Mset{m},3),
                                plm.Rz{m}{c}(:,:,t) = eye(plm.N) - plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                            end
                            prepglm{m}{c}   = @(P,Y)kennedy3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCx{m}{c};
                            
                        case 'manly',
                            prepglm{m}{c}   = @(P,Y)manly(P,Y,m,c,plm); % same as the usual Manly
                            plm.eC {m}{c}   = plm.eCm{m}{c};
                            
                        case 'huh-jhun',
                            plm.Rz {m}{c}   = zeros(plm.N,plm.N,size(plm.Mset{m},3));
                            for t = 1:size(plm.Mset{m},3),
                                plm.Rz{m}{c}(:,:,t) = eye(plm.N) - plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                                [Q,D]       = schur(plm.Rz{m}{c}(:,:,t));
                                D           = abs(diag(D)) < 10*eps;
                                Q(:,D)      = [];
                                if t == 1,
                                    plm.hj{m}{c} = zeros([size(Q) size(plm.Mset{m},3)]);
                                end
                                plm.hj{m}{c}(:,:,t) = Q;
                            end
                            prepglm{m}{c}   = @(P,Y)huhjhun3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCx{m}{c};
                            
                        case 'smith',
                            plm.Rz {m}{c}   = zeros(plm.N,plm.N,size(plm.Mset{m},3));
                            for t = 1:size(plm.Mset{m},3),
                                plm.Rz{m}{c}(:,:,t) = eye(plm.N) - plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                            end
                            prepglm{m}{c}   = @(P,Y)smith3d(P,Y,m,c,plm);
                            plm.eC {m}{c}   = plm.eCm{m}{c};
                    end
                    
                    % Pick a name for the function that will compute the statistic
                    % and the name to save the files later.
                    if opts.pearson,
                        if     plm.rC{m}(c) == 1,
                            fastpiv  {m}{c} = @(M,psi,Y)fastr3d(M,psi,Y,m,c,plm);
                        elseif plm.rC{m}(c) >  1,
                            fastpiv  {m}{c} = @(M,psi,Y)fastrsq3d(M,psi,Y,m,c,plm);
                        end
                    else
                        if     plm.rC{m}(c) == 1 && plm.nVG == 1,
                            fastpiv  {m}{c} = @(M,psi,res)fastt3d(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) >  1 && plm.nVG == 1,
                            fastpiv  {m}{c} = @(M,psi,res)fastf3d(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) == 1 && plm.nVG >  1,
                            fastpiv  {m}{c} = @(M,psi,res)fastv3d(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) >  1 && plm.nVG >  1,
                            fastpiv  {m}{c} = @(M,psi,res)fastg3d(M,psi,res,m,c,plm);
                        end
                    end
                    
                    % MV/CCA/Noperm
                    if opts.MV && ~ opts.approx.noperm,
                        if plm.rC{m}(c) == 1 && any(strcmpi(opts.mvstat,{'auto','hotellingtsq'})),
                            fastmv{m}{c} = @(M,psi,res)fasttsq3d(M,psi,res,m,c,plm);
                        else
                            fastmv{m}{c} = @(M,psi,res)fastq3d(M,psi,res,m,c,plm);
                        end
                    end
                    if opts.CCA || opts.approx.noperm,
                        
                        % Residual forming matrix (Z only)
                        plm.Rz{m}{c} = zeros(plm.N,plm.N,plm.Ysiz{1});
                        if strcmpi(plm.rmethod{m}{c},'noz'),
                            plm.Rz{m}{c} = bsxfun(@plus,eye(plm.N),plm.Rz{m}{c});
                        elseif ~any(strcmpi(plm.rmethod{m}{c},{ ...
                                'still-white','freedman-lane',  ...
                                'kennedy','huh-jhun','smith'})),
                            I = eye(N);
                            for t = 1:plm.Ysiz{1},
                                plm.Rz{m}{c}(:,:,t) = I - plm.Z{m}{c}(:,:,t)*pinv(plm.Z{m}{c}(:,:,t));
                            end
                            clear I;
                        end
                        
                        % Make the 3D dataset & residualise wrt Z
                        plm.Yq{m}{c} = cat(3,plm.Yset{:});
                        plm.Yq{m}{c} = permute(plm.Yq{m}{c},[1 3 2]);
                        for t = 1:plm.Ysiz(1),
                            plm.Yq{m}{c}(:,:,t) = plm.Rz{m}{c}(:,:,t)*plm.Yq{m}{c}(:,:,t);
                        end
                    end
                    
                else % i.e., if not evperdat
                    
                    % Pick the regression/permutation method
                    switch lower(plm.rmethod{m}{c}),
                        
                        case 'noz',
                            prepglm{m}{c}      = @(P,Y)noz(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCx{m}{c};
                            
                        case 'exact',
                            prepglm{m}{c}      = @(P,Y)exact(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCx{m}{c};
                            
                        case 'draper-stoneman',
                            prepglm{m}{c}      = @(P,Y)draperstoneman(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCm{m}{c};
                            
                        case 'still-white',
                            plm.Rz{m}{c}       = eye(plm.N) - plm.Z{m}{c}*pinv(plm.Z{m}{c});
                            prepglm{m}{c}      = @(P,Y)stillwhite(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCx{m}{c};
                            
                        case 'freedman-lane',
                            plm.Hz{m}{c}       = plm.Z{m}{c}*pinv(plm.Z{m}{c});
                            plm.Rz{m}{c}       = eye(plm.N) - plm.Hz{m}{c};
                            prepglm{m}{c}      = @(P,Y)freedmanlane(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCm{m}{c};
                            plm.effrank{m}(c)  = plm.rC{m}(c)*(plm.N-sum(diag(plm.Hz{m}{c})));
                            
                        case 'terbraak',
                            isterbraak         = true;
                            prepglm{m}{c}      = @(P,Y)terbraak(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCm{m}{c};
                            
                        case 'kennedy',
                            plm.Rz {m}{c}      = eye(plm.N) - plm.Z{m}{c}*pinv(plm.Z{m}{c});
                            prepglm{m}{c}      = @(P,Y)kennedy(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCx{m}{c};
                            
                        case 'manly',
                            prepglm{m}{c}      = @(P,Y)manly(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCm{m}{c};
                            
                        case 'huh-jhun',
                            plm.Rz {m}{c}      = eye(plm.N) - plm.Z{m}{c}*pinv(plm.Z{m}{c});
                            [plm.hj{m}{c},D]   = schur(plm.Rz{m}{c});
                            D                  = abs(diag(D)) < 10*eps;
                            plm.hj {m}{c}(:,D) = [];
                            prepglm{m}{c}      = @(P,Y)huhjhun(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCx{m}{c};
                            
                        case 'smith',
                            plm.Rz {m}{c}      = eye(plm.N) - plm.Z{m}{c}*pinv(plm.Z{m}{c});
                            prepglm{m}{c}      = @(P,Y)smith(P,Y,m,c,plm);
                            plm.eC {m}{c}      = plm.eCm{m}{c};
                    end
                    
                    % Pick a name for the function that will compute the statistic
                    % and the name to save the files later.
                    if opts.pearson,
                        if     plm.rC{m}(c) == 1,
                            fastpiv  {m}{c} = @(M,psi,Y)fastr(M,psi,Y,m,c,plm);
                        elseif plm.rC{m}(c) >  1,
                            fastpiv  {m}{c} = @(M,psi,Y)fastrsq(M,psi,Y,m,c,plm);
                        end
                    else
                        if     plm.rC{m}(c) == 1 && plm.nVG == 1,
                            fastpiv  {m}{c} = @(M,psi,res)fastt(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) >  1 && plm.nVG == 1,
                            fastpiv  {m}{c} = @(M,psi,res)fastf(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) == 1 && plm.nVG >  1,
                            fastpiv  {m}{c} = @(M,psi,res)fastv(M,psi,res,m,c,plm);
                        elseif plm.rC{m}(c) >  1 && plm.nVG >  1,
                            fastpiv  {m}{c} = @(M,psi,res)fastg(M,psi,res,m,c,plm);
                        end
                    end
                    
                    % MV/CCA/Noperm
                    if opts.MV  && ~ opts.approx.noperm,
                        if plm.rC{m}(c) == 1 && any(strcmpi(opts.mvstat,{'auto','hotellingtsq'})),
                            fastmv{m}{c} = @(M,psi,res)fasttsq(M,psi,res,m,c,plm);
                        else
                            fastmv{m}{c} = @(M,psi,res)fastq(M,psi,res,m,c,plm);
                        end
                    end
                    if opts.CCA || opts.approx.noperm,
                        
                        % Residual forming matrix (Z only)
                        if strcmpi(plm.rmethod{m}{c},'noz'),
                            plm.Rz{m}{c} = eye(plm.N);
                        elseif ~any(strcmpi(plm.rmethod{m}{c},{ ...
                                'still-white','freedman-lane',  ...
                                'kennedy','huh-jhun','smith'})),
                            plm.Rz{m}{c} = eye(plm.N) - plm.Z{m}{c}*pinv(plm.Z{m}{c});
                        end
                        
                        % Make the 3D dataset & residualise wrt Z
                        plm.Yq{m}{c} = cat(3,plm.Yset{:});
                        for y = 1:plm.nY,
                            plm.Yq{m}{c}(:,:,y) = plm.Rz{m}{c}*plm.Yq{m}{c}(:,:,y);
                        end
                        plm.Yq{m}{c} = permute(plm.Yq{m}{c},[1 3 2]);
                    end
                end
                
                % Pick a name to save the files later.
                if opts.pearson || opts.approx.noperm,
                    if     plm.rC{m}(c) == 1,
                        plm.Gname{m}{c} = '_rstat';
                    elseif plm.rC{m}(c) >  1,
                        plm.Gname{m}{c} = '_rsqstat';
                    end
                else
                    if     plm.rC{m}(c) == 1 && plm.nVG == 1,
                        plm.Gname{m}{c} = '_tstat';
                    elseif plm.rC{m}(c) >  1 && plm.nVG == 1,
                        plm.Gname{m}{c} = '_fstat';
                    elseif plm.rC{m}(c) == 1 && plm.nVG > 1,
                        plm.Gname{m}{c} = '_vstat';
                    elseif plm.rC{m}(c) >  1 && plm.nVG > 1,
                        plm.Gname{m}{c} = '_gstat';
                    end
                end
            end
            
            % Create the permutation set, while taking care of the synchronized
            % permutations (see the outer loop above)
            if opts.syncperms,
                P_inner = po;
                plm.nP{m}(c) = plm.nP{1}(1);
            else
                if ~ opts.approx.noperm,
                    if isempty(plm.EB),
                        if opts.savemetrics,
                            [plm.Pset,plm.nP{m}(c),plm.metr{m}{c}] = ...
                                palm_shuffree(plm.seq{m}{c},opts.nP0, ...
                                opts.cmcp,opts.EE,opts.ISE,false);
                        else
                            [plm.Pset,plm.nP{m}(c)] = ...
                                palm_shuffree(plm.seq{m}{c},opts.nP0, ...
                                opts.cmcp,opts.EE,opts.ISE,false);
                        end
                    else
                        if opts.savemetrics,
                            [plm.Pset,plm.nP{m}(c),plm.metr{m}{c}] = ...
                                palm_shuftree(opts,plm,m,c);
                        else
                            [plm.Pset,plm.nP{m}(c)] = ...
                                palm_shuftree(opts,plm,m,c);
                        end
                    end
                    fprintf('Building null distribution.\n');
                else
                    fprintf('Doing the approximation without permutations.\n');
                end
                P_inner = 1:plm.nP{m}(c);
            end
            
            if po == 1 && ~ opts.approx.noperm,
                % If the user wants to save the permutations, save the vectors now.
                % This has 3 benefits: (1) the if-test below will run just once, rather
                % than many times inside the loop, (2) if the user only wants the
                % vectors, not the images, he/she can cancel immediately after the
                % text file has been created and (3) having all just as a single big
                % file is more convenient than hundreds of small ones.
                if opts.saveperms,
                    % It's faster to write directly as below than using dlmwrite and
                    % palm_swapfmt.m
                    fid = fopen(sprintf('%s%s%s_permidx.csv',opts.o,plm.mstr{m},plm.cstr{c}),'w');
                    for p = 1:plm.nP{m}(c),
                        fprintf(fid,'%d,',palm_perm2idx(plm.Pset{p})');
                        fseek(fid,-1,'eof');
                        fprintf(fid,'\n');
                    end
                    fclose(fid);
                end
                
                % If the user requests, save the permutation metrics
                if opts.savemetrics,
                    fid = fopen(sprintf('%s%s%s_metrics.csv',opts.o,plm.mstr{m},plm.cstr{c}),'w');
                    fprintf(fid,[ ...
                        'Log of max number of permutations given the tree (W),%f\n' ...
                        'Log of max number of permutations if unrestricted (W0),%f\n' ...
                        'Huberman & Hogg complexity (tree only),%d\n' ...
                        'Huberman & Hogg complexity (tree & design),%d\n' ...
                        'Average Hamming distance (tree only),%f\n' ...
                        'Average Hamming distance (tree & design),%f\n' ...
                        'Average Euclidean distance (tree only),%f\n' ...
                        'Average Euclidean distance (tree & design),%f\n' ...
                        'Average Spearman correlation,%f\n'], plm.metr{m}{c});
                    fclose(fid);
                end
                
                % Some vars for later
                if isterbraak, psi0 = cell(plm.nY,1); end
                if opts.approx.negbin, ysel = cell(plm.nY,1); end
                if opts.npcmod && ~ opts.npccon,
                    plm.Tmax{m}{c} = zeros(plm.nP{m}(c),1);
                end
                if opts.MV && ~ opts.approx.noperm,
                    plm.Qmax{m}{c} = zeros(plm.nP{m}(c),1);
                    if opts.approx.lowrank,
                        psiq = zeros(plm.nEV{m}(c),plm.nsel(1),plm.nY);
                        resq = zeros(plm.N,plm.nsel(1),plm.nY);
                    elseif ~ opts.approx.negbin,
                        psiq = zeros(plm.nEV{m}(c),plm.Ysiz(1),plm.nY);
                        resq = zeros(plm.N,plm.Ysiz(1),plm.nY);
                    end
                end
            end
            
            if ~ opts.syncperms,
                ProgressNum = 0;
            end
            
            % Whimsical permutation p-vals using no permutations at all
            if opts.approx.noperm,
                
                % Prepare W for the univariate (no for-loop here):
                if opts.saveunivariate,
                    normY = sum(plm.Yq{m}{c}.^2,1).^.5;
                    u     = bsxfun(@rdivide,plm.Yq{m}{c},normY);
                    
                    % Also, take care of the file names
                    if (plm.nY > 1 || opts.verbosefilenames),
                        if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
                        for y = loopY,
                            plm.ystr{y} = sprintf('_m%d',y);
                        end
                    end
                end
                
                % For four cases below, compute the statistic trace(AW),
                % the moments of the permutation distribution, and fit a gamma:
                % Because this all runs very quickly, the unpermuted
                % statistics can be saved later.
                if opts.evperdat,
                    if opts.saveunivariate,
                        if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
                        for t = 1:plm.Ysiz(y),
                            RzX = plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t);
                            A   = RzX*pinv(RzX);
                            for y = loopY,
                                W   = u(:,y,t)*u(:,y,t)';
                                plm.G{y}{m}{c}(1,t) = trace(A*W);
                                [m1,m2,m3] = palm_moments(A,W,plm.N);
                                plm.Gpperm{y}{m}{c}(1,t) = ...
                                    palm_gamma(plm.G{y}{m}{c}(1,t),m1,m2,m3,false);
                                
                                % If rank(C)=1, the test can be two-tailed,
                                % under the assumption that it's symmetric
                                if opts.ISE && plm.rC{m}(c) == 1,
                                    [M,Y] = prepglm{m}{c}(eye(plm.N),plm.Yset{y});
                                    psi = zeros(size(M,2),plm.Ysiz(y));
                                    for tt = 1:plm.Ysiz(y),
                                        psi(:,tt) = M(:,:,tt)\Y(:,tt);
                                    end
                                    sgn = sign(plm.eC{m}{c}'*psi);
                                    plm.G{y}{m}{c} = plm.G{y}{m}{c}.^.5.*sgn;
                                    plm.Gpperm{y}{m}{c} = plm.Gpperm{y}{m}{c}./2;
                                    plm.Gpperm{y}{m}{c}(sgn<0) = 1 - plm.Gpperm{y}{m}{c}(sgn<0);
                                end
                            end
                        end
                    end
                    if opts.MV,
                        for t = 1:plm.Ysiz(1),
                            RzX     = plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t);
                            A       = RzX*pinv(RzX);
                            [u,~,~] = svd(plm.Yq{m}{c}(:,:,t),'econ');
                            W       = u*u';
                            plm.Q{m}{c}(1,t) = trace(A*W);
                            [m1,m2,m3] = palm_moments(A,W,plm.N);
                            plm.Qpperm{m}{c}(1,t) = ...
                                palm_gamma(plm.Q{m}{c}(1,t),m1,m2,m3,plm.mvrev{m}{c});
                        end
                    end
                else
                    RzX   = plm.Rz{m}{c}*plm.X{m}{c};
                    A     = RzX*pinv(RzX);
                    if opts.saveunivariate,
                        if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
                        for y = loopY,
                            for t = 1:plm.Ysiz(y),
                                
                                W = u(:,y,t)*u(:,y,t)';
                                plm.G{y}{m}{c}(1,t) = trace(A*W);
                                [m1,m2,m3] = palm_moments(A,W,plm.N);
                                plm.Gpperm{y}{m}{c}(1,t) = ...
                                    palm_gamma(plm.G{y}{m}{c}(1,t),m1,m2,m3,false);
                            end
                            
                            % If rank(C)=1, the test can be two-tailed,
                            % under the assumption that it's symmetric
                            if opts.ISE && plm.rC{m}(c) == 1,
                                [M,Y] = prepglm{m}{c}(eye(plm.N),plm.Yset{y});
                                sgn = sign(plm.eC{m}{c}'*(M\Y));
                                plm.G{y}{m}{c} = plm.G{y}{m}{c}.^.5.*sgn;
                                plm.Gpperm{y}{m}{c} = plm.Gpperm{y}{m}{c}./2;
                                plm.Gpperm{y}{m}{c}(sgn<0) = 1 - plm.Gpperm{y}{m}{c}(sgn<0);
                            end
                        end
                    end
                    if opts.MV,
                        for t = 1:plm.Ysiz(1),
                            [u,~,~] = svd(plm.Yq{m}{c}(:,:,t),'econ');
                            W       = u*u';
                            plm.Q{m}{c}(1,t) = trace(A*W);
                            [m1,m2,m3] = palm_moments(A,W,plm.N);
                            plm.Qpperm{m}{c}(1,t) = ...
                                palm_gamma(plm.Q{m}{c}(1,t),m1,m2,m3,plm.mvrev{m}{c});
                        end
                    end
                end
                
                % Ensure that not even the univariate will be done
                P_inner = [];
            end
            
            % For each permutation (inner loop)
            for p = P_inner,
                
                % For each input dataset
                if opts.designperinput, loopY = m; else loopY = 1:plm.nY; end
                for y = loopY,
                    
                    % This if is for the negative binomial mode (if one
                    % modality finished, others my keep running)
                    if dothisY(y),
                        
                        % Some feedback
                        ProgressNum = ProgressNum + 1;
                        if opts.showprogress,
                            if opts.syncperms,
                                fprintf('%0.3g%%\t [Shuffling %d/%d, Design %d/%d, Contrast %d/%d, Modality %d/%d]\n', ...
                                    100*ProgressNum/ProgressDen,...
                                    p,plm.nP{m}(c),m,plm.nM,c,plm.nC(m),y,plm.nY);
                            else
                                fprintf('%0.3g%%\t [Design %d/%d, Contrast %d/%d, Shuffling %d/%d, Modality %d/%d]\n', ...
                                    100*(ProgressNum/plm.nP{m}(c)/plm.nY + ProgressCon)/sum(plm.nC),...
                                    m,plm.nM,c,plm.nC(m),p,plm.nP{m}(c),y,plm.nY);
                            end
                        end

                        % String for the modality index & var for the max:
                        if p == 1,
                            plm.Gmax{y}{m}{c} = zeros(plm.nP{m}(c),1);
                            if (plm.nY > 1 || opts.verbosefilenames),
                                plm.ystr{y} = sprintf('_m%d',y);
                            end
                        end
                        
                        % Shuffle the data and/or design.
                        if opts.approx.negbin,
                            if p == 1,
                                ysel{y} = true(1,plm.Ysiz(y));
                            end
                            [M,Y] = prepglm{m}{c}(plm.Pset{p},plm.Yset{y}(:,ysel{y}));
                        elseif opts.approx.lowrank,
                            if p <= 2*plm.effrank{m}(c),
                                ysel{y} = true(1,plm.Ysiz(y));
                            else
                                if opts.MV,
                                    if y == 1,
                                        ysel{y} = randperm(plm.Ysiz(y));
                                        ysel{y} = ysel{y}(1:plm.nsel(y));
                                    else
                                        ysel{y} = ysel{1};
                                    end
                                else
                                    ysel{y} = randperm(plm.Ysiz(y));
                                    ysel{y} = ysel{y}(1:plm.nsel(y));
                                end
                            end
                            [M,Y] = prepglm{m}{c}(plm.Pset{p},plm.Yset{y}(:,ysel{y}));
                        else
                            [M,Y] = prepglm{m}{c}(plm.Pset{p},plm.Yset{y});
                        end
                        
                        % Do the GLM fit.
                        if opts.evperdat,
                            if opts.approx.lowrank,
                                psi = zeros(size(M,2),plm.nsel(y));
                                res = zeros(size(Y));
                                for t = 1:plm.nsel(y),
                                    psi(:,ysel{y}(t)) = M(:,:,ysel{y}(t))\Y(:,ysel{y}(t));
                                    res(:,ysel{y}(t)) = Y(:,ysel{y}(t)) - M(:,:,ysel{y}(t))*psi(:,ysel{y}(t));
                                end
                            else
                                psi = zeros(size(M,2),plm.Ysiz(y));
                                res = zeros(size(Y));
                                for t = 1:plm.Ysiz(y),
                                    psi(:,t) = M(:,:,t)\Y(:,t);
                                    res(:,t) = Y(:,t) - M(:,:,t)*psi(:,t);
                                end
                            end
                        else
                            psi = M\Y;
                            res = Y - M*psi;
                        end
                        
                        % Unless this is negbin mode, there is no need to fit
                        % again for the MV later
                        if opts.MV && ~ opts.approx.negbin,
                            if opts.approx.lowrank && y == 1,
                                if p == 1,
                                    psiq = zeros(size(psi,1),size(psi,2),plm.nY);
                                    resq = zeros(plm.N,size(psi,2),plm.nY);
                                elseif p == 2*plm.effrank{m}(c) + 1,
                                    psiq = zeros(size(psi,1),plm.nsel(1),plm.nY);
                                    resq = zeros(plm.N,plm.nsel(1),plm.nY);
                                end
                            end
                            psiq(:,:,y) = psi;
                            resq(:,:,y) = res;
                        end
                        
                        % ter Braak permutes under alternative.
                        if isterbraak,
                            if p == 1,
                                psi0{y} = psi;
                            else
                                psi = psi - psi0{y};
                            end
                        end
                        
                        % Compute the pivotal statistic.
                        if opts.pearson,
                            G  {y}{m}{c} = fastpiv{m}{c}(M,psi,Y);
                            df2{y}{m}{c} = NaN;
                        else
                            [G{y}{m}{c},df2{y}{m}{c}] = fastpiv{m}{c}(M,psi,res);
                        end
                        
                        % Save the unpermuted statistic if not z-score
                        if ~ opts.zstat
                            if p == 1,
                                if opts.saveunivariate,
                                    palm_quicksave(G{y}{m}{c},0,opts,plm,y,m,c, ...
                                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                                    
                                    % Save also the degrees of freedom for the unpermuted
                                    if opts.savedof,
                                        if numel(df2{y}{m}{c}) == 1,
                                            savedof(plm.rC{m}(c),df2{y}{m}{c}, ...
                                                horzcat(sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}),'_dof.txt'));
                                        else
                                            savedof(plm.rC{m}(c),mean(df2{y}{m}{c}), ...
                                                horzcat(sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}),'_meandof.txt'));
                                            palm_quicksave(df2{y}{m}{c},0,opts,plm,y,m,c, ...
                                                horzcat(sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}),'_dof'));
                                        end
                                    end
                                end
                            end
                            
                            % Save the stats for each permutation if that was asked
                            if opts.saveperms && ~ opts.approx.negbin && ~ opts.approx.lowrank,
                                palm_quicksave(G{y}{m}{c},0,opts,plm,y,m,c, ...
                                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}),sprintf('_perm%06d',p)));
                            end
                        end
                        
                        % Convert to Z, but make sure that the rank is changed
                        % just once, regardless
                        if p == 1 && (opts.designperinput || y == 1),
                            plm.rC0{m}(c) = plm.rC{m}(c);
                            plm.rC {m}(c) = 0;
                        end
                        G{y}{m}{c} = palm_gtoz(G{y}{m}{c},plm.rC0{m}(c),df2{y}{m}{c});
                        
                        % Save the unpermuted statistic if z-score
                        if opts.zstat
                            if p == 1 && y == 1,
                                plm.Gname{m}{c} = sprintf('_z%s',plm.Gname{m}{c}(2:end));
                                if opts.saveunivariate,
                                    palm_quicksave(G{y}{m}{c},0,opts,plm,y,m,c, ...
                                        sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}));
                                end
                            end
                            
                            % Save the stats for each permutation if that was asked
                            if opts.saveperms && ~ opts.approx.negbin,
                                palm_quicksave(G{y}{m}{c},0,opts,plm,y,m,c, ...
                                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{y},plm.Gname{m}{c},plm.ystr{y},plm.mstr{m},plm.cstr{c}),sprintf('_perm%06d',p)));
                            end
                        end
                        
                        % Draft mode
                        if opts.approx.negbin,
                            
                            % Remove the sign if this is a two-tailed test. This
                            % makes no difference if rank(C) > 1
                            % This needs to be here, inside the if-condition) because of the
                            % lowrank approximation stuff
                            if opts.twotail, G{y}{m}{c} = abs(G{y}{m}{c}); end
                            
                            % In p = 1, there is no counter being incremented (stays at 0) and the number
                            % of permutations performed stays also at 0. In other words, in the negbin mode,
                            % the first permutation is entirely ignored, so that the Haldane equation can 
                            % be applied directly.
                            % If the number of desired exceedances isn't found, this means the extra first
                            % permutation needs to be counted later, and the p-val is then computed as usual.
                            if p == 1,
                                % In the first permutation, keep G and df2,
                                % and start the counter.
                                % In the "negbin" mode, the plm.Gpperm variable isn't really a
                                % counter, but the number of permutations performed until
                                % a certain number of exceedances were found.
                                plm.G      {y}{m}{c} = G  {y}{m}{c};
                                plm.df2    {y}{m}{c} = df2{y}{m}{c};
                                plm.Gppermp{y}{m}{c} = zeros(size(G{y}{m}{c}));
                            else
                                % Otherwise, store the permutation in which a larger
                                % statistic happened, and remove this voxel/vertex/face
                                % from further runs.
                                % These increments only happen for p >= 2, such that the first permutation
                                % isn't counted (needs to be taken care of later).
                                plm.Gpperm{y}{m}{c}(ysel{y}) = plm.Gpperm{y}{m}{c}(ysel{y}) + ...
                                    (G{y}{m}{c} >= plm.G{y}{m}{c}(ysel{y}));
                                plm.Gppermp{y}{m}{c}(ysel{y}) = p - 1; % the -1 is to ignore the 1st perm (unpermuted)
                                ysel{y} = plm.Gpperm{y}{m}{c} < opts.approx.negbin;
                            end
                        else
                            
                            % Low rank approximation
                            if opts.approx.lowrank,
                                if p == 1,
                                    Gperms{y}{m}{c} = zeros(2*plm.effrank{m}(c),plm.Ysiz(y));
                                end
                                if p < 2*plm.effrank{m}(c),
                                    
                                    % Initial permutations (done fully)
                                    Gperms{y}{m}{c}(p,:) = G{y}{m}{c};
                                    
                                elseif p == 2*plm.effrank{m}(c),
                                    
                                    % Generate a new basis
                                    plm.U{y}{m}{c} = palm_lowrank(Gperms{y}{m}{c});
                                    
                                    % Reconstruct past permutations in this
                                    % new basis, and save the moments
                                    [Gperms{y}{m}{c},plm.moments{y}{m}{c}] = ...
                                        palm_lowrank(Gperms{y}{m}{c},plm.U{y}{m}{c},plm.nsel(y));
                                    if opts.twotail,
                                        Gperms{y}{m}{c} = abs(Gperms{y}{m}{c});
                                    end
                        
                                    % Redefine the previous G and counter, now
                                    % using the reconstructed stats
                                    plm.G{y}{m}{c} = Gperms{y}{m}{c}(1,:);
                                    plm.Gpperm{y}{m}{c} = sum(bsxfun(@ge,Gperms{y}{m}{c},plm.G{y}{m}{c}),1);
                                    plm.Gmax  {y}{m}{c}(1:p) = max(Gperms{y}{m}{c},[],2);
                                    Gperms{y}{m}{c} = []; % free up a bit of memory
                                    
                                else
                                    
                                    % Once a basis is known, use it.
                                    G{y}{m}{c} = palm_lowrank(G{y}{m}{c},plm.U{y}{m}{c},...
                                        ysel{y},plm.moments{y}{m}{c});
                                end
                            end
                            
                            % This needs to be here, inside the if-condition) because of the
                            % lowrank approximation stuff
                            if opts.twotail, G{y}{m}{c} = abs(G{y}{m}{c}); end
                            
                            % In the first permutation, keep G and df2,
                            % and start the counter.
                            if p == 1,
                                plm.G  {y}{m}{c} = G  {y}{m}{c};
                                plm.df2{y}{m}{c} = df2{y}{m}{c};
                            end
                            if ~ opts.approx.lowrank || p > 2*plm.effrank{m}(c),
                                plm.Gpperm{y}{m}{c}    = plm.Gpperm{y}{m}{c} + (G{y}{m}{c} >= plm.G{y}{m}{c});
                                plm.Gmax  {y}{m}{c}(p) = max(G{y}{m}{c},[],2);
                            end
                            
                            % Tail and gamma approximations
                            if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                if p == 1,
                                    plm.Gperms{y}{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(y));
                                end
                                plm.Gperms{y}{m}{c}(p,:) = G{y}{m}{c};
                            end
 
                            % Cluster statistic is here
                            if opts.cluster.uni.do,
                                if p == 1,
                                    plm.Gclumax{y}{m}{c} = zeros(plm.nP{m}(c),1);
                                    [plm.Gclumax{y}{m}{c}(p),plm.Gclu{y}{m}{c}] = clusterfunc( ...
                                        G{y}{m}{c},y,opts.cluster.uni.thr,opts,plm);
                                else
                                    plm.Gclumax{y}{m}{c}(p) = clusterfunc( ...
                                        G{y}{m}{c},y,opts.cluster.uni.thr,opts,plm);
                                end
                            end
                            
                            % TFCE is here
                            if opts.tfce.uni.do,
                                Gtfce{y}{m}{c} = tfcefunc(G{y}{m}{c},y,opts,plm);
                                if p == 1,
                                    plm.Gtfcemax  {y}{m}{c} = zeros(plm.nP{m}(c),1);
                                    plm.Gtfce     {y}{m}{c} = Gtfce{y}{m}{c};
                                    plm.Gtfcepperm{y}{m}{c} = zeros(size(G{y}{m}{c}));
                                end
                                plm.Gtfcepperm{y}{m}{c} = plm.Gtfcepperm{y}{m}{c} + ...
                                    (Gtfce{y}{m}{c} >= plm.Gtfce{y}{m}{c});
                                plm.Gtfcemax{y}{m}{c}(p) = max(Gtfce{y}{m}{c},[],2);
                                
                                % Tail and gamma approximations
                                if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                    if p == 1,
                                        plm.Gtfceperms{y}{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(y));
                                    end
                                    plm.Gtfceperms{y}{m}{c}(p,:) = Gtfce{y}{m}{c};
                                end
                            end
                        end
                        if opts.approx.negbin && ~ any(ysel{y}),
                            dothisY(y) = false;
                        end
                    end
                end
                
                % NPC for Y only is here
                if opts.npcmod && ~ opts.npccon && ~ opts.designperinput,
                    % Note that one design per input can't be run here,
                    % because this is the inner loop, and can't combine the
                    % designs here; needs to be outside.
                    if opts.showprogress,
                        fprintf('\t [Combining modalities, Design %d/%d, Contrast %d/%d]\n',m,plm.nM,c,plm.nC(m));
                    end
                    
                    % Just a feedback message for some situations.
                    if opts.showprogress && ...
                            ~ opts.zstat && ...
                            p == 1 && ...
                            opts.savepara && ...
                            ~ plm.nonpcppara && ...
                            ~ opts.spatial.npc && ...
                            any(strcmpi(opts.npcmethod,{ ...
                            'dudbridge-koeleman', ...
                            'dudbridge-koeleman2'})),
                        fprintf('(1st perm is slower) ');
                    end
                    
                    % Compute the combined statistic
                    for y = 1:plm.nY,
                        Gnpc  {1}(y,:) = G  {y}{m}{c};
                        df2npc{1}(y,:) = df2{y}{m}{c};
                    end
                    T{m}{c} = fastnpc(Gnpc{1},0,df2npc{1});
                    
                    % Since computing the parametric p-value for some methods
                    % can be quite slow, it's faster to run all these checks
                    % to ensure that 'pparanpc' runs just once.
                    if opts.zstat || ...
                            opts.spatial.npc || ...
                            (p == 1 && opts.savepara && ~ plm.nonpcppara),
                        Tppara{m}{c} = pparanpc(T{m}{c},plm.nY);
                        
                        % Adjust the concordant signs for the concordant
                        % test. This is the same as 1-(1-P)^2, but
                        % preserves the precision.
                        if opts.concordant,
                            Tppara{m}{c} = 2*Tppara{m}{c}-Tppara{m}{c}.^2;
                        end
                        
                        % Reserve the p-parametric to save later.
                        if p == 1,
                            plm.Tppara{m}{c} = Tppara{m}{c};
                        end
                    end
                    
                    % Convert T to zstat if that was asked (note that at this point,
                    % G was already converted to z before making T).
                    if opts.zstat,
                        T{m}{c} = -erfinv(2*Tppara{m}{c}-1)*sqrt(2);
                        if p == 1 && m == 1 && c == 1,
                            plm.npcstr = horzcat('_z',plm.npcstr(2:end));
                        end
                    end
                    
                    % Save the NPC Statistic (this is inside the loop because
                    % of the two-tailed option)
                    if p == 1,
                        palm_quicksave(T{m}{c},0,opts,plm,[],m,c, ...
                            sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,plm.mstr{m},plm.cstr{c}));
                    end
                    
                    % If the user wants to save the NPC statistic for each
                    % permutation, save it now.
                    if opts.saveperms,
                        palm_quicksave(T{m}{c},0,opts,plm,[],m,c, ...
                            horzcat(sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,plm.mstr{m},plm.cstr{c}),sprintf('_perm%06d',p)));
                    end
                    
                    % Increment counters
                    if p == 1,
                        plm.T{m}{c} = T{m}{c};
                        plm.Tpperm{m}{c} = zeros(size(T{m}{c}));
                    end
                    plm.Tpperm{m}{c} = plm.Tpperm{m}{c} + ...
                        bsxfun(npcrel,T{m}{c},plm.T{m}{c});
                    plm.Tmax{m}{c}(p) = npcextr(T{m}{c},[],2);
                    
                    % Tail and gamma approximations
                    if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                        if p == 1,
                            plm.Tperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                        end
                        plm.Tperms{m}{c}(p,:) = T{m}{c};
                    end
                    
                    % Be sure to use z-scores for the spatial statistics, converting
                    % it if not already.
                    if opts.spatial.npc,
                        if ~ opts.zstat,
                            
                            % Avoid infinities
                            Tppara{m}{c} = Tppara{m}{c}*0.999999999999999 + 1e-15;
                            T{m}{c}      = -erfinv(2*Tppara{m}{c}-1)*sqrt(2);
                        end
                    end
                    
                    % Cluster statistic NPC is here
                    if opts.cluster.npc.do,
                        if p == 1,
                            [plm.Tclumax{m}{c}(p),plm.Tclu{m}{c}] = ...
                                clusterfunc(T{m}{c},1,opts.cluster.npc.thr,opts,plm);
                        else
                            plm.Tclumax{m}{c}(p) = ...
                                clusterfunc(T{m}{c},1,opts.cluster.npc.thr,opts,plm);
                        end
                    end
                    
                    % TFCE NPC is here
                    if opts.tfce.npc.do,
                        Ttfce{m}{c} = tfcefunc(T{m}{c},1,opts,plm);
                        if p == 1,
                            plm.Ttfcemax  {m}{c} = zeros(plm.nP{m}(c),1);
                            plm.Ttfce     {m}{c} = Ttfce{m}{c};
                            plm.Ttfcepperm{m}{c} = zeros(size(T{m}{c}));
                        end
                        plm.Ttfcepperm{m}{c} = plm.Ttfcepperm{m}{c} + ...
                            (Ttfce{m}{c} >= plm.Ttfce{m}{c});
                        plm.Ttfcemax{m}{c}(p) = max(Ttfce{m}{c},[],2);
                        
                        % Tail and gamma approximations
                        if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                            if p == 1,
                                plm.Ttfceperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                            end
                            plm.Ttfceperms{m}{c}(p,:) = Ttfce{m}{c};
                        end
                    end
                end
                
                % MANOVA/MANCOVA is here. CCA is in the elseif below
                if opts.MV && ~ opts.approx.noperm,
                    
                    % This if is for the negative binomial mode.
                    if dotheMVorCCA,
                        if opts.showprogress,
                            if opts.syncperms,
                                fprintf('\t [Shuffling %d/%d, Design %d/%d, Contrast %d/%d, Multivariate]\n', ...
                                    p,plm.nP{m}(c),m,plm.nM,c,plm.nC(m));
                            else
                                fprintf('\t [Design %d/%d, Contrast %d/%d, Shuffling %d/%d, Multivariate]\n', ...
                                    m,plm.nM,c,plm.nC(m),p,plm.nP{m}(c));
                            end
                        end
                        
                        % Shuffle the data and/or design.
                        if opts.approx.negbin,
                            if p == 1,
                                yselq = true(1,size(plm.Yq{m}{c},2),1);
                            end
                            if any(yselq),
                                psiq = zeros(plm.nEV{m}(c),sum(yselq),plm.nY);
                                resq = zeros(plm.N,size(psiq,2),plm.nY);
                                for y = 1:plm.nY,
                                    [M,Y] = prepglm{m}{c}(plm.Pset{p},plm.Yq{m}{c}(:,yselq,y));
                                    psiq(:,:,y) = M\Y;
                                    resq(:,:,y) = Y - M*psiq(:,:,y);
                                end
                            end
                        end
                        
                        % ter Braak permutes under alternative.
                        if isterbraak,
                            if p == 1,
                                psiq0 = psiq;
                            else
                                psiq  = psiq - psiq0;
                            end
                        end
                        
                        % Compute the pivotal multivariate statistic.
                        Q{m}{c} = fastmv{m}{c}(M,psiq,resq);
                        
                        % Since computing the parametric p-value for some methods
                        % can be quite slow, it's faster to run all these checks
                        % to ensure that 'pparamv' runs just once.
                        if opts.zstat            || ...
                                opts.spatial.mv  || ...
                                (p == 1          && ...
                                opts.savepara    && ...
                                ~ plm.nomvppara),
                            Qppara{m}{c} = pparamv{m}{c}(Q{m}{c});
                            
                            % Reserve the p-parametric to save later.
                            if p == 1,
                                plm.Qppara{m}{c} = Qppara{m}{c};
                            end
                        end
                        
                        % Convert to zstat if that was asked
                        if opts.zstat,
                            Q{m}{c} = -erfinv(2*Qppara{m}{c}-1)*sqrt(2);
                            if p == 1 && m == 1 && c == 1,
                                plm.mvstr = horzcat('_z',plm.mvstr(2:end));
                            end
                        end
                        
                        % Save the MV statistic
                        if p == 1,
                            palm_quicksave(Q{m}{c},0,opts,plm,[],m,c, ...
                                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}));
                        end
                        
                        % Draft mode
                        if opts.approx.negbin,
                            % In p = 1, there is no counter being incremented (stays at 0) and the number
                            % of permutations performed stays also at 0. In other words, in the negbin mode,
                            % the first permutation is entirely ignored, so that the Haldane equation can 
                            % be applied directly.
                            % If the number of desired exceedances isn't found, this means the extra first
                            % permutation needs to be counted later, and the p-val is then computed as usual.
                            if p == 1,
                                % In the first permutation, keep Q and Qdf2,
                                % and start the counter.
                                % In the "negbin" mode, the plm.Qpperm variable isn't a counter,
                                % but the number of permutations until a statistic larger than
                                % the unpermuted was found.
                                plm.Q      {m}{c} = Q{m}{c};
                                plm.Qdf2   {m}{c} = Qdf2{m}{c};
                                plm.Qpperm {m}{c} = zeros(size(Q{m}{c}));
                                plm.Qppermp{m}{c} = zeros(size(Q{m}{c}));
                                
                            else
                                % Otherwise, store the permutation in which a larger
                                % statistic happened, and remove this voxel/vertex/face
                                % from further runs.
                                % These increments only happen for p >= 2, such that the first permutation
                                % isn't counted (needs to be taken care of later).
                                plm.Qpperm{m}{c}(yselq) = plm.Qpperm{m}{c}(yselq) + ...
                                    mvrel(Q{m}{c},plm.Q{m}{c}(yselq));
                                plm.Qppermp{m}{c}(yselq) = p - 1; % the -1 is to ignore the 1st perm (unpermuted)
                                yselq = plm.Qpperm{m}{c} < opts.approx.negbin;
                            end
                        else
                            
                            % Low rank approximation
                            if opts.approx.lowrank,
                                if p == 1,
                                    Qperms{m}{c} = zeros(2*plm.effrank{m}(c),plm.Ysiz(1));
                                end
                                if p < 2*plm.effrank{m}(c),
                                    
                                    % Initial permutations (done fully)
                                    Qperms{m}{c}(p,:) = Q{m}{c};
                                    
                                elseif p == 2*plm.effrank{m}(c),
                                    
                                    % Generate a new basis
                                    plm.Uq{m}{c} = palm_lowrank(Qperms{m}{c});
                                    
                                    % Reconstruct past permutations in this
                                    % new basis, and save the moments
                                    [Qperms{m}{c},plm.momentsq{m}{c}] = ...
                                        palm_lowrank(Qperms{m}{c},plm.Uq{m}{c},plm.nsel(1));
                                    
                                    % Redefine the previous Q and counter, now
                                    % using the reconstructed stats
                                    plm.Q{m}{c} = Qperms{m}{c}(1,:);
                                    plm.Qpperm{m}{c} = sum(bsxfun(@ge,Qperms{m}{c},plm.Q{m}{c}),1);
                                    plm.Qmax  {m}{c}(1:p) = max(Qperms{m}{c},[],2);
                                    Qperms    {m}{c} = []; % free up a bit of memory
                                    
                                else
                                    
                                    % Once a basis is known, use it.
                                    Q{m}{c} = palm_lowrank(Q{m}{c},plm.Uq{m}{c},...
                                        ysel{1},plm.momentsq{m}{c});
                                end
                            end
                            
                            % If the user wants to save the statistic for each
                            % permutation, save it now. This isn't obviously allowed
                            % in negbin mode, as the images are not complete. Also,
                            % this is inside the loop to allow the two-tailed option
                            % not to use to much memory
                            if opts.saveperms,
                                palm_quicksave(Q{m}{c},0,opts,plm,[],m,c, ...
                                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}),sprintf('_perm%06d',p)));
                            end
                            if p == 1,
                                % In the first permutation, keep Q and start the counter.
                                plm.Q     {m}{c}  = Q{m}{c};
                                plm.Qpperm{m}{c}  = zeros(size(Q{m}{c}));
                            end
                            plm.Qpperm{m}{c}      = plm.Qpperm{m}{c} + mvrel(Q{m}{c},plm.Q{m}{c});
                            plm.Qmax  {m}{c}(p)   = mvextr(Q{m}{c},[],2);
                            
                            % Tail and gamma approximations
                            if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                if p == 1,
                                    plm.Qperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                                end
                                plm.Qperms{m}{c}(p,:) = Q{m}{c};
                            end
                            
                            % Now compute the spatial statistics, converting to z-score
                            % if not already.
                            if opts.spatial.mv,
                                if ~ opts.zstat,
                                    
                                    % Avoid infinities (this is ugly but needed).
                                    Qppara{m}{c}  = Qppara{m}{c}*0.999999999999999 + 1e-15;
                                    Q{m}{c}       = -erfinv(2*Qppara{m}{c}-1)*sqrt(2);
                                end
                            end
                            
                            % Cluster statistic is here
                            if opts.cluster.mv.do,
                                if p == 1,
                                    plm.Qclumax{m}{c} = zeros(plm.nP{m}(c),1);
                                    [plm.Qclumax{m}{c}(p),plm.Qclu{m}{c}] = ...
                                        clusterfunc(Q{m}{c},1,opts.cluster.mv.thr,opts,plm);
                                else
                                    plm.Qclumax{m}{c}(p) = ...
                                        clusterfunc(Q{m}{c},1,opts.cluster.mv.thr,opts,plm);
                                end
                            end
                            
                            % TFCE is here
                            if opts.tfce.mv.do,
                                Qtfce{m}{c} = tfcefunc(Q{m}{c},1,opts,plm);
                                if p == 1,
                                    plm.Qtfcemax  {m}{c} = zeros(plm.nP{m}(c),1);
                                    plm.Qtfce     {m}{c} = Qtfce{m}{c};
                                    plm.Qtfcepperm{m}{c} = zeros(size(Q{m}{c}));
                                end
                                plm.Qtfcepperm{m}{c} = plm.Qtfcepperm{m}{c} + ...
                                    (Qtfce{m}{c} >= plm.Qtfce{m}{c});
                                plm.Qtfcemax{m}{c}(p) = max(Qtfce{m}{c},[],2);
                                
                                % Tail and gamma approximations
                                if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                    if p == 1,
                                        plm.Qtfceperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                                    end
                                    plm.Qtfceperms{m}{c}(p,:) = Qtfce{m}{c};
                                end
                            end
                        end
                        if opts.approx.negbin && ~ any(yselq),
                            dotheMVorCCA = false;
                        end
                    end
                    
                elseif opts.CCA,
                    
                    % This if is for the negative binomial mode.
                    if dotheMVorCCA,
                        if opts.showprogress,
                            if opts.syncperms,
                                fprintf('\t [Shuffling %d/%d, Design %d/%d, Contrast %d/%d, CCA]\n', ...
                                    p,plm.nP{m}(c),m,plm.nM,c,plm.nC(m));
                            else
                                fprintf('\t [Design %d/%d, Contrast %d/%d, Shuffling %d/%d, CCA]\n', ...
                                    m,plm.nM,c,plm.nC(m),p,plm.nP{m}(c));
                            end
                        end
                        
                        % Initialise some vars
                        if p == 1,
                            yselq   = true(1,1,plm.Ysiz(1));
                            Q{m}{c} = zeros(1,plm.Ysiz(1));
                        end
                        
                        % Compute the CC coefficient
                        if opts.evperdat,
                            for t = find(yselq),
                                M(:,:,t)   = plm.Pset{p}*plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t);
                                Q{m}{c}(t) = cca(plm.Yq{m}{c}(:,:,t),M(:,:,t),opts.ccaparm);
                            end
                        else
                            M = plm.Pset{p}*plm.Rz{m}{c}*plm.X{m}{c};
                            for t = find(yselq),
                                Q{m}{c}(t) = cca(plm.Yq{m}{c}(:,:,t),M,opts.ccaparm);
                            end
                        end
                        
                        % Convert to zstat if that was asked
                        if opts.zstat,
                            Q{m}{c}(yselq) = atanh(Q{m}{c}(yselq));
                            if p == 1 && m == 1 && c == 1,
                                plm.mvstr = horzcat('_z',plm.mvstr(2:end));
                            end
                        end

                        % Save the CCA statistic
                        if p == 1,
                            palm_quicksave(Q{m}{c},0,opts,plm,[],m,c, ...
                                sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}));
                        end
                        
                        % In p = 1, there is no counter being incremented (stays at 0) and the number
                        % of permutations performed stays also at 0. In other words, in the negbin mode,
                        % the first permutation is entirely ignored, so that the Haldane equation can
                        % be applied directly.
                        % If the number of desired exceedances isn't found, this means the extra first
                        % permutation needs to be counted later, and the p-val is then computed as usual.
                        if opts.approx.negbin,
                            if p == 1,
                                % In the first permutation, keep Q and Qdf2,
                                % and start the counter.
                                % In the "negbin" mode, the plm.Qpperm variable isn't a counter,
                                % but the number of permutations until a statistic larger than
                                % the unpermuted was found.
                                plm.Q      {m}{c} = Q{m}{c};
                                plm.Qpperm {m}{c} = zeros(size(Q{m}{c}));
                                plm.Qppermp{m}{c} = zeros(size(Q{m}{c}));
                            else
                                % Otherwise, store the permutation in which a larger
                                % statistic happened, and remove this voxel/vertex/face
                                % from further runs.
                                % These increments only happen for p >= 2, such that the first permutation
                                % isn't counted (needs to be taken care of later).
                                plm.Qpperm{m}{c}(yselq) = plm.Qpperm{m}{c}(yselq) + ...
                                    (Q{m}{c}(yselq) >= plm.Q{m}{c}(yselq));
                                plm.Qppermp{m}{c}(yselq) = p - 1;  % the -1 is to ignore the 1st perm (unpermuted)
                                yselq = plm.Qpperm{m}{c} < opts.approx.negbin;
                            end
                        else
                            
                            % Low rank approximation
                            if opts.approx.lowrank,
                                if p == 1,
                                    Qperms{m}{c} = zeros(2*plm.effrank{m}(c),plm.Ysiz(1));
                                end
                                if p < 2*plm.effrank{m}(c),
                                    
                                    % Initial permutations (done fully)
                                    Qperms{m}{c}(p,:) = Q{m}{c};
                                    
                                elseif p == 2*plm.effrank{m}(c),
                                    
                                    % Generate a new basis
                                    plm.Uq{m}{c} = palm_lowrank(Qperms{m}{c});
                                    
                                    % Reconstruct past permutations in this
                                    % new basis, and save the moments
                                    [Qperms{m}{c},plm.momentsq{m}{c}] = ...
                                        palm_lowrank(Qperms{m}{c},plm.Uq{m}{c},plm.nsel(1));
                                    
                                    % Redefine the previous Q and counter, now
                                    % using the reconstructed stats
                                    plm.Q{m}{c} = Qperms{m}{c}(1,:);
                                    plm.Qpperm{m}{c} = sum(bsxfun(@ge,Qperms{m}{c},plm.Q{m}{c}),1);
                                    plm.Qmax  {m}{c}(1:p) = max(Qperms{m}{c},[],2);
                                    Qperms    {m}{c} = []; % free up a bit of memory
                                    
                                else
                                    
                                    % Once a basis is known, use it.
                                    Q{m}{c} = palm_lowrank(Q{m}{c},plm.Uq{m}{c},...
                                        ysel{1},plm.momentsq{m}{c});
                                end
                            end
                            
                            % If the user wants to save the statistic for each
                            % permutation, save it now. This isn't obviously allowed
                            % in negbin mode, as the images are not complete. Also,
                            % this is inside the loop to allow the two-tailed option
                            % not to use to much memory
                            if opts.saveperms,
                                palm_quicksave(Q{m}{c},0,opts,plm,[],m,c, ...
                                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{1},plm.mvstr,plm.Qname{m}{c},plm.mstr{m},plm.cstr{c}),sprintf('_perm%06d',p)));
                            end
                            
                            % In the first permutation, keep Q and start the counter.
                            if p == 1,
                                plm.Q     {m}{c} = Q{m}{c};
                                plm.Qpperm{m}{c} = zeros(size(Q{m}{c}));
                            end
                            plm.Qpperm{m}{c}     = plm.Qpperm{m}{c} + (Q{m}{c} >= plm.Q{m}{c});
                            plm.Qmax  {m}{c}(p)  = max(Q{m}{c},[],2);
                            
                            % Tail and gamma approximations
                            if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                if p == 1,
                                    plm.Qperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                                end
                                plm.Qperms{m}{c}(p,:) = Q{m}{c};
                            end
                            
                            % Now compute the MV spatial statistics, converting to z-score
                            % if not already.
                            if opts.spatial.mv,
                                if ~ opts.zstat,
                                    Q{m}{c} = atanh(Q{m}{c});
                                end
                            end
                            
                            % Cluster statistic is here
                            if opts.cluster.mv.do,
                                if p == 1,
                                    [plm.Qclumax{m}{c}(p),plm.Qclu{m}{c}] = ...
                                        clusterfunc(Q{m}{c},1,opts.cluster.mv.thr,opts,plm);
                                else
                                    plm.Qclumax{m}{c}(p) = ...
                                        clusterfunc(Q{m}{c},1,opts.cluster.mv.thr,opts,plm);
                                end
                            end
                            
                            % TFCE is here
                            if opts.tfce.mv.do,
                                Qtfce{m}{c} = tfcefunc(Q{m}{c},1,opts,plm);
                                if p == 1,
                                    plm.Qtfce     {m}{c} = Qtfce{m}{c};
                                    plm.Qtfcepperm{m}{c} = zeros(size(Q{m}{c}));
                                end
                                plm.Qtfcepperm{m}{c} = plm.Qtfcepperm{m}{c} + ...
                                    (Qtfce{m}{c} >= plm.Qtfce{m}{c});
                                plm.Qtfcemax{m}{c}(p) = max(Qtfce{m}{c},[],2);
                                
                                % Tail and gamma approximations
                                if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                                    if p == 1,
                                        plm.Qtfceperms{m}{c} = zeros(plm.nP{m}(c),plm.Ysiz(1));
                                    end
                                    plm.Qtfceperms{m}{c}(p,:) = Qtfce{m}{c};
                                end
                            end
                        end
                        if opts.approx.negbin && ~ any(yselq),
                            dotheMVorCCA = false;
                        end
                    end
                end
            end
            if ~ opts.syncperms,
                ProgressCon = ProgressCon + 1;
            end
        end
    end
    
    % NPC for modalities with design per input is here
    if opts.npcmod && opts.designperinput && ~ opts.npccon,
        
        % For each contrast (here we collapse designs and modalities, pairwise)
        for c = 1:plm.nC(1),
            if opts.showprogress,
                fprintf('\t [Combining modalities (pairwise designs), Contrast %d/%d]\n',c,plm.nC(1));
            end
            
            % Just a feedback message for some situations.
            if opts.showprogress && ...
                    ~ opts.zstat && ...
                    p == 1 && ...
                    opts.savepara && ...
                    ~ plm.nonpcppara && ...
                    ~ opts.spatial.npc && ...
                    any(strcmpi(opts.npcmethod,{ ...
                    'dudbridge-koeleman', ...
                    'dudbridge-koeleman2'})),
                fprintf('(1st perm is slower) ');
            end
            
            % Compute the combined statistic
            for y = 1:plm.nY,
                Gnpc  {1}(y,:) = G  {y}{y}{c};
                df2npc{1}(y,:) = df2{y}{y}{c};
            end
            T{1}{c} = fastnpc(Gnpc{1},0,df2npc{1});
            
            % Since computing the parametric p-value for some methods
            % can be quite slow, it's faster to run all these checks
            % to ensure that 'pparanpc' runs just once.
            if opts.zstat || ...
                    opts.spatial.npc || ...
                    (p == 1 && opts.savepara && ~ plm.nonpcppara),
                Tppara{1}{c} = pparanpc(T{1}{c},plm.nY);
                
                % Adjust the concordant signs for the concordant
                % test. This is the same as 1-(1-P)^2, but
                % preserves the precision.
                if opts.concordant,
                    Tppara{1}{c} = 2*Tppara{1}{c}-Tppara{1}{c}.^2;
                end
                
                % Reserve the p-parametric to save later.
                if p == 1,
                    plm.Tppara{1}{c} = Tppara{1}{c};
                end
            end
            
            % Convert T to zstat if that was asked (note that at this point,
            % G was already converted to z before making T).
            if opts.zstat,
                T{1}{c} = -erfinv(2*Tppara{1}{c}-1)*sqrt(2);
                if p == 1 && c == 1,
                    plm.npcstr = horzcat('_z',plm.npcstr(2:end));
                end
            end
            
            % Save the NPC Statistic (this is inside the loop because
            % of the two-tailed option)
            if p == 1,
                palm_quicksave(T{1}{c},0,opts,plm,[],1,c, ...
                    sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,plm.mstr{1},plm.cstr{c}));
            end
            
            % If the user wants to save the NPC statistic for each
            % permutation, save it now.
            if opts.saveperms,
                palm_quicksave(T{1}{c},0,opts,plm,[],1,c, ...
                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{1},plm.npcstr,plm.Tname,plm.mstr{1},plm.cstr{c}),sprintf('_perm%06d',p)));
            end
            
            % Increment counters
            if p == 1,
                plm.T{1}{c} = T{1}{c};
                plm.Tpperm{1}{c} = zeros(size(T{1}{c}));
            end
            plm.Tpperm{1}{c} = plm.Tpperm{1}{c} + ...
                bsxfun(npcrel,T{1}{c},plm.T{1}{c});
            plm.Tmax{1}{c}(p) = npcextr(T{1}{c},[],2);
            
            % Tail and gamma approximations
            if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                if p == 1,
                    plm.Tperms{1}{c} = zeros(plm.nP{1}(c),plm.Ysiz(1));
                end
                plm.Tperms{1}{c}(p,:) = T{1}{c};
            end
            
            % Be sure to use z-scores for the spatial statistics, converting
            % it if not already.
            if opts.spatial.npc,
                if ~ opts.zstat,
                    
                    % Avoid infinities
                    Tppara{1}{c} = Tppara{1}{c}*0.999999999999999 + 1e-15;
                    T{1}{c}      = -erfinv(2*Tppara{1}{c}-1)*sqrt(2);
                end
            end
            
            % Cluster statistic NPC is here
            if opts.cluster.npc.do,
                if p == 1,
                    [plm.Tclumax{1}{c}(p),plm.Tclu{1}{c}] = ...
                        clusterfunc(T{1}{c},1,opts.cluster.npc.thr,opts,plm);
                else
                    plm.Tclumax{1}{c}(p) = ...
                        clusterfunc(T{1}{c},1,opts.cluster.npc.thr,opts,plm);
                end
            end
            
            % TFCE NPC is here
            if opts.tfce.npc.do,
                Ttfce{1}{c} = tfcefunc(T{1}{c},1,opts,plm);
                if p == 1,
                    plm.Ttfcemax  {1}{c} = zeros(plm.nP{1}(c),1);
                    plm.Ttfce     {1}{c} = Ttfce{1}{c};
                    plm.Ttfcepperm{1}{c} = zeros(size(T{1}{c}));
                end
                plm.Ttfcepperm{1}{c} = plm.Ttfcepperm{1}{c} + ...
                    (Ttfce{1}{c} >= plm.Ttfce{1}{c});
                plm.Ttfcemax{1}{c}(p) = max(Ttfce{1}{c},[],2);
                
                % Tail and gamma approximations
                if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                    if p == 1,
                        plm.Ttfceperms{1}{c} = zeros(plm.nP{1}(c),plm.Ysiz(1));
                    end
                    plm.Ttfceperms{1}{c}(p,:) = Ttfce{1}{c};
                end
            end
        end
    end
    
    % NPC for contrasts is here
    if opts.npccon,
        
        % Just a feedback message for some situations.
        if opts.showprogress && ...
                ~ opts.zstat && ...
                po == 1 && ...
                opts.savepara && ...
                ~ plm.nonpcppara && ...
                ~ opts.spatial.npc && ...
                any(strcmpi(opts.npcmethod,{ ...
                'dudbridge-koeleman', ...
                'dudbridge-koeleman2'})),
            fprintf('(1st perm is slower) ');
        end
        
        % Assemble the stats for the combination
        if opts.npcmod,
            if opts.showprogress,
                fprintf('\t [Combining modalities and all contrasts]\n');
            end
            j = 1;
            for y = 1:plm.nY,
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        Gnpc  {1}(j,:) = G  {y}{m}{c};
                        df2npc{1}(j,:) = df2{y}{m}{c};
                        j = j + 1;
                    end
                end
            end
            if po == 1,
                plm.Tmax{1} = zeros(plm.nP{1}(1),1);
            end
        else
            if opts.showprogress,
                fprintf('\t [Combining all contrasts]\n');
            end
            for y = 1:plm.nY,
                j = 1;
                if opts.designperinput, loopM = y; else loopM = 1:plm.nM; end
                for m = loopM,
                    for c = 1:plm.nC(m),
                        Gnpc  {y}(j,:) = G  {y}{m}{c};
                        df2npc{y}(j,:) = df2{y}{m}{c};
                        j = j + 1;
                    end
                end
                if po == 1,
                    plm.Tmax{y} = zeros(plm.nP{1}(1),1);
                end
            end
        end
        
        % For each of set to be combined
        for j = 1:numel(Gnpc),
            
            % String with the counter
            if po == 1
                if j == 1,
                    plm.jstr = cell(numel(Gnpc),1);
                end
                if ~ opts.npcmod && (plm.nY > 1 || opts.verbosefilenames),
                    plm.jstr{j} = sprintf('_m%d',j);
                else
                    plm.jstr{j} = '';
                end
            end
            
            % Compute the combined statistic
            T{j} = fastnpc(Gnpc{j},0,df2npc{j});
            
            % Since computing the parametric p-value for some methods
            % can be quite slow, it's faster to run all these checks
            % to ensure that 'pparanpc' runs just once.
            if opts.zstat || ...
                    opts.spatial.npc || ...
                    (po == 1 && opts.savepara && ~ plm.nonpcppara),
                Tppara{j} = pparanpc(T{j});
                
                % Adjust the concordant signs for the concordant
                % test. This is the same as 1-(1-P)^2, but
                % preserves the precision.
                if opts.concordant,
                    Tppara{j} = 2*Tppara{j}-Tppara{j}.^2;
                end
                
                % Reserve the p-parametric to save later.
                if po == 1,
                    plm.Tppara{j} = Tppara{j};
                end
            end
            
            % Convert T to zstat if that was asked (note that at this point,
            % G was already converted to z before making T).
            if opts.zstat,
                T{j} = -erfinv(2*Tppara{j}-1)*sqrt(2);
                if po == 1 && j == 1,
                    plm.npcstr = horzcat('_z',plm.npcstr(2:end));
                end
            end
            
            % Save the NPC Statistic (this is inside the loop because
            % of the two-tailed option)
            if po == 1,
                palm_quicksave(T{j},0,opts,plm,j,[],[], ...
                    sprintf('%s',opts.o,plm.Ykindstr{j},plm.npcstr,plm.Tname,plm.jstr{j}));
            end
            
            % If the user wants to save the NPC statistic for each
            % permutation, save it now.
            if opts.saveperms,
                palm_quicksave(T{j},0,opts,plm,j,[],[], ...
                    horzcat(sprintf('%s',opts.o,plm.Ykindstr{j},plm.npcstr,plm.Tname,plm.jstr{j}),sprintf('_perm%06d',po)));
            end
            
            % Increment counters
            if po == 1,
                plm.T{j} = T{j};
                plm.Tpperm{j} = zeros(size(T{j}));
            end
            plm.Tpperm{j} = plm.Tpperm{j} + ...
                bsxfun(npcrel,T{j},plm.T{j});
            plm.Tmax{j}(po) = npcextr(T{j},[],2);
            
            % Tail and gamma approximations
            if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                if p == 1,
                    plm.Tperms{j} = zeros(plm.nP{1}(1),size(T{j},2));
                end
                plm.Tperms{j}(po,:) = T{j};
            end
            
            % Be sure to use z-scores for the spatial statistics, converting
            % it if not already.
            if opts.spatial.npc,
                if ~ opts.zstat,
                    
                    % Avoid infinities
                    Tppara{j} = Tppara{j}*0.999999999999999 + 1e-15;
                    T{j}      = -erfinv(2*Tppara{j}-1)*sqrt(2);
                end
            end
            
            % Cluster statistic NPC is here
            if opts.cluster.npc.do,
                if p == 1,
                    [plm.Tclumax{j}(po),plm.Tclu{j}] = ...
                        clusterfunc(T{j},1,opts.cluster.npc.thr,opts,plm);
                else
                    plm.Tclumax{j}(po) = ...
                        clusterfunc(T{j},1,opts.cluster.npc.thr,opts,plm);
                end
            end
            
            % TFCE NPC is here
            if opts.tfce.npc.do,
                Ttfce{j} = tfcefunc(T{j},1,opts,plm);
                if p == 1,
                    plm.Ttfcemax  {j} = zeros(plm.nP{1}(1),1);
                    plm.Ttfce     {j} = Ttfce{j};
                    plm.Ttfcepperm{j} = zeros(size(T{j}));
                end
                plm.Ttfcepperm{j} = plm.Ttfcepperm{j} + ...
                    (Ttfce{j} >= plm.Ttfce{j});
                plm.Ttfcemax{j}(po) = max(Ttfce{j},[],2);
                
                % Tail and gamma approximations
                if opts.saveuncorrected && (opts.approx.tail || opts.approx.gamma),
                    if p == 1,
                        plm.Ttfceperms{j} = zeros(plm.nP{1}(1),size(T{j},2));
                    end
                    plm.Ttfceperms{j}(po,:) = Ttfce{j};
                end
            end
        end
    end
end
tocP = toc(ticP);
fprintf('Elapsed time with permutations: ~ %g s\n',tocP);

% Free up a bit of memory after the loop.
clear M Y psi res G df2 T Q;

% Save everything, except the few bits saved above
ticS = tic;
palm_saveall(plm,opts);
tocS = toc(ticS);
fprintf('Elapsed time generating and saving results: ~ %g s\n',tocS);
fprintf('Overall elapsed time: ~ %g s\n',tocI+tocP+tocS);
csvwrite(sprintf('%s_elapsed.csv',opts.o),[tocI tocP tocS]);

% Finished.
fprintf('PALM finished at %s.\n',datestr(now));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  F U N C T I O N S  %%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================
% Below are the functions for each of the regression methods:
% ==============================================================
function [Mr,Y] = noz(P,Y,m,c,plm)
% This is equivalent to Draper-Stoneman, as when there is no Z
% Y remains unchanged.
Mr = P*plm.X{m}{c};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = noz3d(P,Y,m,c,plm)
Mr = zeros(size(plm.X{m}{c}));
for t = 1:size(Y,2),
    Mr(:,:,t) = P*plm.X{m}{c}(:,:,t);
end

% ==============================================================
function [Mr,Yr] = exact(P,Y,m,c,plm)
% The "exact" method, in which the coefficients for
% the nuisance are known.
Yr = Y - plm.Z{m}{c}*plm.g;
Mr = P*plm.X{m}{c};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = exact3d(P,Y,m,c,plm)
Mr = zeros(size(plm.X{m}{c}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t)   = Y(:,t) - plm.Z{m}{c}(:,:,t)*plm.g;
    Mr(:,:,t) = P*plm.X{m}{c}(:,:,t);
end

% ==============================================================
function [Mr,Y] = draperstoneman(P,Y,m,c,plm)
% Draper and Stoneman (1966) method.
% Y remains unchanged
Mr = horzcat(P*plm.X{m}{c},plm.Z{m}{c});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = draperstoneman3d(P,Y,m,c,plm)
Mr = zeros(size(plm.Mp{m}{c}));
for t = 1:size(Y,2),
    Mr(:,:,t) = horzcat(P*plm.X{m}{c}(:,:,t),plm.Z{m}{c}(:,:,t));
end

% ==============================================================
function [Mr,Yr] = stillwhite(P,Y,m,c,plm)
% A method following the same logic as the one
% proposed by Still and White (1981)
Yr = plm.Rz{m}{c}*Y;
Mr = P*plm.X{m}{c};
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = stillwhite3d(P,Y,m,c,plm)
Mr = zeros(size(plm.X{m}{c}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t)   = plm.Rz{m}{c}*Y(:,t);
    Mr(:,:,t) = P*plm.X{m}{c}(:,:,t);
end

% ==============================================================
function [Mr,Yr] = freedmanlane(P,Y,m,c,plm)
% The Freedman and Lane (1983) method.
Mr = plm.Mp{m}{c};
Yr = (P'*plm.Rz{m}{c} + plm.Hz{m}{c})*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = freedmanlane3d(P,Y,m,c,plm)
Mr = plm.Mp{m}{c};
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t) = (P'*plm.Rz{m}{c}(:,:,t) + plm.Hz{m}{c}(:,:,t))*Y(:,t);
end

% ==============================================================
function [Mr,Yr] = manly(P,Y,m,c,plm)
% The Manly (1986) method.
% There's no need for a 3D version of this method.
Mr = plm.Mp{m}{c};
Yr = P'*Y;

% ==============================================================
function [Mr,Yr] = terbraak(P,Y,m,c,plm)
% The ter Braak (1992) method.
Mr = plm.Mp{m}{c};
Yr = (P'*plm.Rm{m}{c} + plm.Hm{m}{c})*Y; % original method
% Yr = P'*plm.Rm{m}{c}*Y; % alternative (causes unpermuted stat to be 0)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = terbraak3d(P,Y,m,c,plm)
Mr = plm.Mp{m}{c};
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Yr(:,t) = (P'*plm.Rm{m}{c}(:,:,t) + plm.Hm{m}{c}(:,:,t))*Y(:,t); % original method
end

% ==============================================================
function [Mr,Yr] = kennedy(P,Y,m,c,plm)
% The Kennedy (1996) method. This method should NEVER be used.
Mr = plm.Rz{m}{c}*plm.X{m}{c};
Yr = P'*plm.Rz{m}{c}*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = kennedy3d(P,Y,m,c,plm)
Mr = zeros(size(plm.X{m}{c}));
Yr = zeros(size(Y));
for t = 1:size(Y,2),
    Mr(:,:,t) = plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t);
    Yr(:,t) = P'*plm.Rz{m}{c}(:,:,t)*Y(:,t);
end

% ==============================================================
function [Mr,Yr] = huhjhun(P,Y,m,c,plm)
% The Huh and Jhun (2001) method, that fixes the issues
% with Kennedy's, but doesn't allow block permutation.
Mr = plm.hj{m}{c}'*plm.Rz{m}{c}*plm.X{m}{c};
Yr = P'*plm.hj{m}{c}'*plm.Rz{m}{c}*Y;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Yr] = huhjhun3d(P,Y,m,c,plm)
Mr = zeros(size(plm.hj{m}{c},1),size(plm.X{m}{c},2));
Yr = zeros(size(plm.hj{m}{c},1),size(Y,2));
for t = 1:size(Y,2),
    Mr(:,:,t) = plm.hj{m}{c}(:,:,t)'*plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t);
    Yr(:,t)   = P'*plm.hj{m}{c}(:,:,t)'*plm.Rz{m}{c}(:,:,t)*Y(:,t);
end

% ==============================================================
function [Mr,Y] = smith(P,Y,m,c,plm)
% The Smith method, i.e., orthogonalization.
% Y remains unchanged
Mr = horzcat(P*plm.Rz{m}{c}*plm.X{m}{c},plm.Z{m}{c});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [Mr,Y] = smith3d(P,Y,m,c,plm)
Mr = zeros(size(plm.Mp{m}{c}));
for t = 1:size(Y,2),
    Mr(:,:,t) = horzcat(P*plm.Rz{m}{c}(:,:,t)*plm.X{m}{c}(:,:,t),plm.Z{m}{c}(:,:,t));
end

% ==============================================================
% Below are the functions to compute univariate statistics:
% ==============================================================
function G = fastr(M,psi,Y,m,c,plm)
% This only works if:
% - M and Y have zero mean.
% - rank(contrast) = 1
%
% Inputs:
% M   : design matrix (demeaned)
% psi : regression coefficients
% Y   : data (demeaned)
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Pearson's correlation coefficient (r).
G = fastrsq(M,psi,Y,m,c,plm);
G = sign(plm.eC{m}{c}'*psi).*G.^.5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function G = fastr3d(M,psi,Y,m,c,plm)
G = fastrsq3d(M,psi,Y,m,c,plm);
G = sign(plm.eC{m}{c}'*psi).*G.^.5;

% ==============================================================
function G = fastrsq(M,psi,Y,m,c,plm)
% This only works if:
% - M and Y have zero mean.
%
% Inputs:
% M   : design matrix (demeaned)
% psi : regression coefficients
% Y   : data (demeaned)
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : R^2, i.e., the coefficient of determination.
tmp = plm.eC{m}{c}/(plm.eC{m}{c}'/(M'*M)*plm.eC{m}{c})*plm.eC{m}{c}';
G   = sum((tmp'*psi).*psi,1);
den = sum(Y.^2,1);
G   = G./den;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function G = fastrsq3d(M,psi,Y,m,c,plm)
for t = 1:size(psi,2),
    tmp = plm.eC{m}{c}/(plm.eC{m}{c}'/(M(:,:,t)'*M(:,:,t))*plm.eC{m}{c})*plm.eC{m}{c}';
    G   = sum((tmp'*psi(:,t)).*psi(:,t),1);
end
den = sum(Y.^2,1);
G   = G./den;

% ==============================================================
function [G,df2] = fastt(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : t statistic.
% df2 : Degrees of freedom. df1 is 1 for the t statistic.
df2 = plm.N-plm.rM{m}(c);
G   = plm.eC{m}{c}'*psi;
den = sqrt(plm.eC{m}{c}'/(M'*M)*plm.eC{m}{c}*sum(res.^2)./df2);
G   = G./den;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastt3d(M,psi,res,m,c,plm)
df2 = plm.N-plm.rM{m}(c);
G   = plm.eC{m}{c}'*psi;
MtM = zeros(1,size(psi,2));
for t = 1:size(psi,2),
    MtM(t) = plm.eC{m}{c}'/(M(:,:,t)'*M(:,:,t))*plm.eC{m}{c};
end
den = sqrt(MtM.*sum(res.^2,1)./df2);
G   = G./den;

% ==============================================================
function [G,df2] = fastf(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups = 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : F-statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).
df2 = plm.N-plm.rM{m}(c);
cte = plm.eC{m}{c}/(plm.eC{m}{c}'/(M'*M)*plm.eC{m}{c})*plm.eC{m}{c}';
tmp = zeros(size(psi));
for j = 1:size(cte,2),
    tmp(j,:) = sum(bsxfun(@times,psi,cte(:,j)),1)';
end
G   = sum(tmp.*psi,1);
ete = sum(res.^2,1);
G   = G./ete*df2/plm.rC{m}(c);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastf3d(M,psi,res,m,c,plm)
df2 = plm.N-plm.rM{m}(c);
nT = size(res,2);
cte = zeros(size(psi,1),size(psi,1),nT);
for t = 1:nT,
    cte(:,:,t) = plm.eC{m}{c}/(plm.eC{m}{c}' ...
        /(M(:,:,t)'*M(:,:,t))*plm.eC{m}{c})*plm.eC{m}{c}';
end
ppsi = permute(psi,[1 3 2]);
ppsi = sum(bsxfun(@times,ppsi,cte),1);
ppsi = permute(ppsi,[2 3 1]);
G    = sum(ppsi.*psi,1);
ete  = sum(res.^2,1);
G    = G./ete*df2/plm.rC{m}(c);

% ==============================================================
function [G,df2] = fastv(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups > 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Aspin-Welch v statistic.
% df2 : Degrees of freedom 2. df1 is 1.
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
den  = zeros(1,nT);
r    = size(M,2);
dRmb = zeros(plm.nVG,1);
cte  = zeros(r^2,nT);
for b = 1:plm.nVG,
    bidx    = plm.VG == b;
    dRmb(b) = sum(plm.dRm{m}{c}(bidx),1);
    W(b,:)  = dRmb(b)./sum(res(bidx,:).^2,1);
    Mb      = M(bidx,:)'*M(bidx,:);
    cte     = cte + Mb(:)*W(b,:);
    W(b,:)  = W(b,:)*sum(bidx);
end
for t = 1:nT,
    den(t) = plm.eC{m}{c}'/reshape(cte(:,t),[r r])*plm.eC{m}{c};
end
G = plm.eC{m}{c}'*psi./sqrt(den);
sW1 = sum(W,1);
bsum = sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,W,sW1)).^2,dRmb),1);
df2 = 1/3./bsum;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastv3d(M,psi,res,m,c,plm)
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
den  = zeros(1,nT);
r    = size(M,2);
dRmb = zeros(plm.nVG,nT);
cte  = zeros(r,r,nT);
for b = 1:plm.nVG,
    bidx = plm.VG == b;
    dRmb(b,:) = sum(plm.dRm{m}{c}(bidx,:),1);
    W(b,:)    = dRmb(b,:)./sum(res(bidx,:).^2,1);
    for t = 1:nT,
        cte(:,:,t) = cte(:,:,t) + (M(bidx,:,t)'*M(bidx,:,t)).*W(b,t);
    end
    W(b,:) = W(b,:)*sum(bidx);
end
for t = 1:nT,
    den(t) = plm.eC{m}{c}'/cte(:,:,t)*plm.eC{m}{c};
end
G    = plm.eC{m}{c}'*psi./sqrt(den);
sW1  = sum(W,1);
bsum = sum((1-bsxfun(@rdivide,W,sW1)).^2./dRmb,1);
df2  = 1/3./bsum;

% ==============================================================
function [G,df2] = fastg(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups > 1
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% G   : Welch v^2 statistic.
% df2 : Degrees of freedom 2. df1 is rank(C).
r    = size(M,2);
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
dRmb = zeros(plm.nVG,1);
cte  = zeros(r^2,nT);
for b = 1:plm.nVG,
    bidx    = plm.VG == b;
    dRmb(b) = sum(plm.dRm{m}{c}(bidx));
    W(b,:)  = dRmb(b)./sum(res(bidx,:).^2);
    Mb      = M(bidx,:)'*M(bidx,:);
    cte     = cte + Mb(:)*W(b,:);
    W(b,:)  = W(b,:)*sum(bidx);
end
G = zeros(1,nT);
for t = 1:nT,
    A = psi(:,t)'*plm.eC{m}{c};
    G(t) = A/(plm.eC{m}{c}'/(reshape(cte(:,t),[r r]))* ...
        plm.eC{m}{c})*A'/plm.rC{m}(c);
end
sW1  = sum(W,1);
bsum = sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,W,sW1)).^2,dRmb),1);
bsum = bsum/plm.rC{m}(c)/(plm.rC{m}(c)+2);
df2  = 1/3./bsum;
G    = G./(1 + 2*(plm.rC{m}(c)-1).*bsum);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [G,df2] = fastg3d(M,psi,res,m,c,plm)
r    = size(M,2);
nT   = size(res,2);
W    = zeros(plm.nVG,nT);
dRmb = zeros(plm.nVG,nT);
cte  = zeros(r,r,nT);
for b = 1:plm.nVG,
    bidx = plm.VG == b;
    dRmb(b,:) = sum(plm.dRm{m}{c}(bidx,:),1);
    W(b,:) = dRmb(b,:)./sum(res(bidx,:).^2,1);
    for t = 1:nT,
        cte(:,:,t) = cte(:,:,t) + (M(bidx,:,t)'*M(bidx,:,t))*W(b,t);
    end
    W(b,:) = W(b,:)*sum(bidx);
end
G = zeros(1,nT);
for t = 1:nT,
    A = psi(:,t)'*plm.eC{m}{c};
    G(t) = A/(plm.eC{m}{c}'/cte(:,:,t) * ...
        plm.eC{m}{c})*A'/plm.rC{m}(c);
end
sW1  = sum(W,1);
bsum = sum((1-bsxfun(@rdivide,W,sW1)).^2./dRmb,1);
bsum = bsum/plm.rC{m}(c)/(plm.rC{m}(c)+2);
df2  = 1/3./bsum;
G    = G./(1 + 2*(plm.rC{m}(c)-1).*bsum);

% ==============================================================
% Below are the functions to compute multivariate statistics:
% ==============================================================
function Q = fasttsq(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) = 1
% - number of variance groups = 1
% - psi and res are 3D
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% Q    : Hotelling's T^2 statistic.

% Swap dimensions so that dims 1 and 2 are subjects and variables
% leaving the voxels/tests as the 3rd.
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT  = size(res,3);
df0 = plm.N-plm.rM{m}(c);
S = spr(res)/df0;
Q = zeros(1,nT);
cte2 = plm.eC{m}{c}'/(M'*M)*plm.eC{m}{c};
for t = 1:nT,
    cte1 = plm.eC{m}{c}'*psi(:,:,t)*plm.Dset{m}{c};
    Q(1,t) = cte1/(plm.Dset{m}{c}'*S(:,:,t)*plm.Dset{m}{c})/cte2*cte1';
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Q = fasttsq3d(M,psi,res,m,c,plm)
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT  = size(res,3);
df0 = plm.N-plm.rM{m}(c);
S = spr(res)/df0;
Q = zeros(1,nT);
for t = 1:nT,
    cte1 = plm.eC{m}{c}'*psi(:,:,t)*plm.Dset{m}{c};
    cte2 = plm.eC{m}{c}'/(M(:,:,t)'*M(:,:,t))*plm.eC{m}{c};
    Q(1,t) = cte1/(plm.Dset{m}{c}'*S(:,:,t)*plm.Dset{m}{c})/cte2*cte1';
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fasttsqp(Q,df2,p)
% P-value for Hotelling's T^2
P = palm_gpval(Q*(df2-p+1)/p/df2,p,df2-p+1);

% ==============================================================
function Q = fastq(M,psi,res,m,c,plm)
% This works only if:
% - rank(contrast) > 1
% - number of variance groups = 1
% - psi and res are 3D
%
% Inputs:
% M   : design matrix
% psi : regression coefficients
% res : residuals
% plm : a struct with many things as generated by
%       'palm_backend.m' and 'palm_takeargs.m'
%
% Outputs:
% Q    : Multivariate (yet scalar) statistic.

% Swap dimensions so that dims 1 and 2 are subjects and variables
% leaving the voxels/tests as the 3rd.
res = permute(res,[1 3 2]);
psi = permute(psi,[1 3 2]);
nT   = size(res,3);
cte2 = plm.eC{m}{c}'/(M'*M)*plm.eC{m}{c};
E    = spr(res);
Q    = zeros(1,nT);
for t = 1:nT,
    cte1   = plm.Dset{m}{c}'*psi(:,:,t)'*plm.eC{m}{c};
    H      = cte1/cte2*cte1';
    Q(1,t) = plm.qfun(plm.Dset{m}{c}'*E(:,:,t)*plm.Dset{m}{c},H);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Q = fastq3d(M,psi,res,m,c,plm)
res  = permute(res,[1 3 2]);
psi  = permute(psi,[1 3 2]);
nT   = size(res,3);
E    = spr(res);
Q    = zeros(1,nT);
for t = 1:nT,
    cte1   = plm.Dset{m}{c}'*psi(:,:,t)'*plm.eC{m}{c};
    cte2   = plm.eC{m}{c}'/(M(:,:,t)'*M(:,:,t))*plm.eC{m}{c};
    H      = cte1/cte2*cte1';
    Q(1,t) = plm.qfun(plm.Dset{m}{c}'*E(:,:,t)*plm.Dset{m}{c},H);
end

% ==============================================================
function Q = wilks(E,H)
% Wilks' lambda.
Q = det(E)/det(E+H);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = wilksp(Q,df1,df2,p)
r = df2-(p-df1+1)/2;
u = (p*df1-2)/4;
cden = (p^2+df1^2-5);
if cden > 0,
    t = sqrt((p^2*df1^2-4)/cden);
else
    t = 1;
end
F = (r*t-2*u)*(1-Q.^(1/t))./(Q.^(1/t)*p*df1);
Fdf1 = p*df1;
Fdf2 = r*t-2*u;
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = lawley(E,H)
% Lawley-Hotelling's trace.
Q = trace(H/E);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = lawleyp(Q,df1,df2,p)
m = (abs(p-df1)-1)/2;
n = (df2-p-1)/2;
s = min(p,df1);
if n > 0,
    b = (p+2*n)*(df1+2*n)/(2*(2*n+1)*(n-1));
    c = (2+(p*df1+2)/(b-1))/(2*n);
    Fdf1 = p*df1;
    Fdf2 = 4+(p*df1+2)/(b-1);
    F = (Q/c)*Fdf2/Fdf1;
else
    Fdf1 = s*(2*m+s+1);
    Fdf2 = 2*(s*n+1);
    F = (Q/s)*Fdf2/Fdf1;
end
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = pillai(E,H)
% Pillai's trace.
Q = trace(H/(E+H));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = pillaip(Q,df1,df2,p)
m = (abs(p-df1)-1)/2;
n = (df2-p-1)/2;
s = min(p,df1);
F = (2*n+s+1)/(2*m+s+1)*(Q./(s-Q));
Fdf1 = s*(2*m+s+1);
Fdf2 = s*(2*n+s+1);
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = roy_ii(E,H)
% Roy's (ii) largest root (analogous to F).
Q = max(eig(H/E));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = roy_iip(Q,df1,df2,p)
Fdf1 = max(p,df1);
Fdf2 = df2-Fdf1+df1;
F = Q*Fdf2/Fdf1;
P = palm_gpval(F,Fdf1,Fdf2);

% ==============================================================
function Q = roy_iii(E,H)
% Roy's (iii) largest root (analogous to R^2).
% No p-vals for this (not even approximate or bound).
Q = max(eig(H/(E+H)));

% ==============================================================
function cc = cca(Y,X,k)
% Do CCA via QR & SVD.
% The ranks of X and Y aren't checked for speed.
% Inputs are assumed to have been mean-centered and be free
% of nuisance (partial CCA) via Y=Rz*Y and X=Rz*X.
% k is the k-th CC (typically we want the 1st).
% Based on the algorithm proposed by:
% Bjorck A, Golub GH. Numerical methods for
% computing angles between linear subspaces.
% Math Comput. 1973;27(123):579-579.
[Qy,~]  = qr(Y,0);
[Qx,~]  = qr(X,0);
[~,D,~] = svd(Qy'*Qx,0);
cc      = max(min(D(k,k),1),0);

% ==============================================================
% Below are the functions to combine statistics:
% ==============================================================
function T = tippett(G,df1,df2)
T = min(palm_gpval(G,df1,df2),[],1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = tippettp(T,nG)
%P = T.^plm.nY;
% Note it can't be simply P = 1-(1-T)^K when implementing
% because precision is lost if the original T is smaller than eps,
% something quite common. Hence the need for the Pascal
% triangle, etc, as done below.
pw  = nG:-1:1;
cf  = pascaltri(nG);
sgn = (-1)*(-1).^pw;
P   = sgn.*cf*bsxfun(@power,T,pw');

% ==============================================================
function T = fisher(G,df1,df2)
T = -2*sum(log(palm_gpval(G,df1,df2)),1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fisherp(T,nG)
P = palm_gpval(T,-1,2*nG);

% ==============================================================
function T = stouffer(G,df1,df2)
T = sum(palm_gtoz(G,df1,df2),1)/sqrt(size(G,1));

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = stoufferp(T,~)
P = palm_gpval(T,0);

% ==============================================================
function T = wilkinson(G,df1,df2,parma)
T = sum(palm_gpval(G,df1,df2) <= parma);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = wilkinsonp(T,nG,parma)
lfac    = palm_factorial(nG);
lalpha  = log(parma);
l1alpha = log(1-parma);
P = zeros(size(T));
for k = 1:nG,
    lp1 = lfac(nG+1) - lfac(k+1) - lfac(nG-k+1);
    lp2 = k*lalpha;
    lp3 = (nG-k)*l1alpha;
    P = P + (k>=T).*exp(lp1+lp2+lp3);
end

% ==============================================================
function T = winer(G,df1,df2)
df2 = bsxfun(@times,ones(size(G)),df2);
cte = sqrt(sum(df2./(df2-2),1));
gp  = palm_gpval(G,df1,df2);
gt  = sign(gp-.5).*sqrt(df2./betainv(2*min(gp,1-gp),df2/2,.5)-df2); % =tinv(gp,df2)
T   = -sum(gt)./cte;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = winerp(T,~)
P = palm_gcdf(-T,0);

% ==============================================================
function T = edgington(G,df1,df2)
T = sum(palm_gpval(G,df1,df2),1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = edgingtonp(T,nG)
lfac = palm_factorial(nG);
fT   = floor(T);
mxfT = max(fT(:));
P = zeros(size(T));
for j = 0:mxfT,
    p1  = (-1)^j;
    lp2 = - lfac(j+1) - lfac(nG-j+1);
    lp3 = nG*log(T-j);
    P = P + (j<=fT).*p1.*exp(lp2+lp3);
end

% ==============================================================
function T = mudholkargeorge(G,df1,df2)
nG = size(G,1);
mhcte = sqrt(3*(5*nG+4)/nG/(5*nG+2))/pi;
T = mhcte*sum(log(...
    palm_gcdf(G,df1,df2)./...
    palm_gpval(G,df1,df2)),1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = mudholkargeorgep(T,nG)
P = palm_gpval(T,1,5*nG+4);

% ==============================================================
function [T,Gpval] = friston(G,df1,df2)
Gpval = palm_gpval(G,df1,df2);
T = max(Gpval,[],1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = fristonp(T,nG,parmu)
P = T.^(nG - parmu + 1);

% ==============================================================
function T = darlingtonhayes(G,df1,df2,parmr)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= parmr;
G       = reshape(G(idx),horzcat(parmr,size(G,2)));
df2     = reshape(df2(idx),horzcat(parmr,size(df2,2)));
P       = palm_gcdf(G,df1,df2);
Z       = erfinv(2*P-1)*sqrt(2);
T       = mean(Z,1);

% ==============================================================
function T = zaykin(G,df1,df2,parma)
P = -log10(palm_gpval(G,df1,df2));
P(P < -log10(parma)) = 0;
T = sum(P,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = zaykinp(T,nG,parma)
lT      = -T;
lfac    = palm_factorial(nG);
lalpha  = log10(parma);
l1alpha = log10(1-parma);
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

% ==============================================================
function T = dudbridgekoeleman(G,df1,df2,parmr)
df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= parmr;
G       = reshape(G(idx),horzcat(parmr,size(G,2)));
df2     = reshape(df2(idx),horzcat(parmr,size(df2,2)));
P       = -log10(palm_gpval(G,df1,df2));
T       = sum(P,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = dudbridgekoelemanp(T,nG,parmr)
lT = -T;
lfac = palm_factorial(nG);
P    = zeros(size(lT));
lp1  = lfac(nG+1)  - ...
    lfac(parmr+2)  - ...
    lfac(nG-parmr) + ...
    log10(parmr+2);
for v = 1:numel(lT);
    P(v) = quad(@(t)dkint(t,lp1,lT(v),nG,...
        parmr,lfac(1:parmr)),eps,1);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function T = dudbridgekoeleman2(G,df1,df2,parmr,parma)
df2 = bsxfun(@times,ones(size(G)),df2);
P = -log10(palm_gpval(G,df1,df2));
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
P(tmp > parmr) = 0;
P(P < -log10(parma)) = 0;
T = sum(P,1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function P = dudbridgekoeleman2p(T,nG,parmr,parma)
lT = -T;
lfac = palm_factorial(nG);
P    = zeros(1,size(T,2));
for k = 1:parmr,
    kk = (nG-k)*log(1-parma);
    if isnan(kk), kk = 0; end
    p1 = exp(lfac(nG+1) - lfac(k+1) - lfac(nG-k+1) + kk);
    p2 = awtk(lT,parma,k,lfac(1:k));
    P = P + p1.*p2;
end
if k < nG,
    lp1 = lfac(nG+1)   - ...
        lfac(parmr+2)  - ...
        lfac(nG-parmr) + ...
        log(parmr+2);
    for v = 1:numel(lT);
        P(v) = P(v) + ...
            quad(@(t)dkint(t,lp1,lT(v),nG,parmr, ...
            lfac(1:parmr)),eps,parma);
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function A = awtk(lw,t,k,lfac)
ltk = k.*log(t);
tk = real(exp(ltk));
s = (1:k)';
L = bsxfun(@plus,lw,...
    bsxfun(@minus,(s-1)*log(k*log(t)-lw),lfac(s)));
S = sum(real(exp(L)),1);
A = (lw <= ltk).*S + (lw > ltk).*tk;

% ==============================================================
function T = taylortibshirani(G,df1,df2)
nG = size(G,1);
P = palm_gpval(G,df1,df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum(1-P.*(nG+1)./prank)/nG;

function P = taylortibshiranip(T,nG)
P = palm_gcdf(-T./sqrt(nG),0);

% ==============================================================
function T = jiang(G,df1,df2,parma)
nG = size(G,1);
P = palm_gpval(G,df1,df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum((P<=parma).*(1-P.*(nG+1)./prank))/nG;

% ==============================================================
% Other useful functions:
% ==============================================================
function savedof(df1,df2,fname)
% Save the degrees of freedom.
% This is faster than dlmwrite.

fdof = fopen(fname,'w');
fprintf(fdof,'%g\n',df1);
fprintf(fdof,'%g,',df2);
fseek(fdof,-1,'cof');
fprintf(fdof,'\n');
fclose(fdof);

% ==============================================================
function S = spr(X)
% Compute the matrix with the sum of products.
% X is a 3D array, with the resilduals of the GLM.
% - 1st dimension are the subjects
% - 2nd dimension the modalities.
% - 3rd dimension would tipically be voxels
%
% S is the sum of products that make up the covariance
% matrix:
% - 1st and 3rd dimension have the same size as the number of
%   modalities and the 2nd dimension are typically the voxels.

% To make it faster, the check should be made just once, and
% the result kept throughout runs.
persistent useway1;
if isempty(useway1),
    
    % Test both ways and compute the timings.
    tic; S1 = way1(X); w1 = toc;
    tic; S2 = way2(X); w2 = toc;
    
    % The variables sp1 and sp2 should be absolutely
    % identical but they may have sightly different numerical
    % precisions so to be consistent, choose the same that will
    % be used for all permutations later
    if w1 < w2,
        useway1 = true;
        S = S1;
    else
        useway1 = false;
        S = S2;
    end
else
    if useway1,
        S = way1(X);
    else
        S = way2(X);
    end
end

function sp = way1(X)
% Way 1: this tends to be faster for Octave and if
% the number of levels in X is smaller than about 5.
[~,nY,nT] = size(X);
sp = zeros(nY,nY,nT);
for y1 = 1:nY,
    for y2 = 1:y1,
        sp(y1,y2,:) = sum(X(:,y1,:).*X(:,y2,:),1);
        if y1 ~= y2,
            sp(y2,y1,:) = sp(y1,y2,:);
        end
    end
end

function sp = way2(X)
% Way 2: This tends to be faster in Matlab or if
% there are many levels in X, e.g., more than about 7.
[~,nY,nT] = size(X);
sp = zeros(nY,nY,nT);
for t = 1:nT,
    sp(:,:,t) = (X(:,:,t)'*X(:,:,t));
end

% ==============================================================
function C = pascaltri(K)
% Returns the coefficients for a binomial expansion of
% power K, except the last term. This is used by the Tippett
% method to avoid issues with numerical precision.

persistent Cp;
if isempty(Cp),
    K = K + 1;
    if K <= 2,
        Cp = horzcat(ones(1,K),0);
    elseif K >= 3,
        Rprev = [1 1 0];
        for r = 3:K,
            Cp = horzcat(Rprev + fliplr(Rprev),0);
            Rprev = Cp;
        end
    end
end
C = Cp(1:end-2);

% Finished! :-)
