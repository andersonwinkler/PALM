% #!/usr/bin/octave -q
function palm(varargin)
% ===========================================
% PALM: Permutation Analysis of Linear Models
% ===========================================
%
% PALM performs permutation inference for the general linear
% models (GLMs) of arbitrary complexity, taking as inputs data
% in various formats, and being able to take into account certain
% cases of well-structured non-independence.
%
% The options are:
%
% -i <file>         : Input(s). More than one can be specified.
% -m <file>         : Mask(s). Either one for all inputs, or one per input,
%                     supplied in the same order as the respective -i appear.
% -s <file>         : Surface file(s), for the inputs that are scalars over
%                     a surface geometry. Each -s should be supplied in the
%                     same order as the respective -i.
% -d <file>         : Design matrix. It can be in .csv format, or in FSL's
%                     "VEST" format.
% -t <file>         : t-contrasts file, in VEST format.
% -f <file>         : F-contrasts file, in VEST format, just as in FSL.
% -eb <file>        : Exchangeability blocks file, in VEST format.
% -vg <file>        : Variance groups file, in VEST format.
% -o <string>       : Output prefix. It may itself be prefixed by a path.
% -n <integer>      : Number of permutations. Use -n 0 for exhaustive.
% -c <real>         : Theshold for cluster extent for t contrasts.
% -C <real>         : Theshold for cluster mass for t contrasts.
% -F <positive>     : Theshold for cluster extent for F contrasts.
% -S <positive>     : Theshold for cluster mass for F contrasts.
% -T                : Enables TFCE inference for 3D (volume) data, i.e.,
%                     with H=2, E=0.5, C=6.
% -T2               : Enables TFCE inference for 2D (surface, TBSS) data,
%                     i.e., H=2, E=1, C=26.
% -cnpc <real>      : Threshold for NPC cluster extent (z-stat).
% -Cnpc <real>      : Threshold for NPC cluster mass (z-stat).
% -Tnpc             : Enables TFCE inference for NPC on 3D data,
%                     with H=2, E=0.5, C=6.
% -T2npc            : Enables TFCE inference for NPC on 2D data,
%                     i.e., H=2, E=1, C=26.
% -tfce_H <real>    : Set the TFCE H (height power) parameter.
% -tfce_E <real>    : Set the TFCE E (extent power) parameter.
% -tfce_C <integer> : Set the TFCE C (connectivity) parameter.
% -pb               : Permute blocks as a whole.
% -ise              : Assume independent and symmetric errors, to
%                     allow sign-flipping.
% -ee               : Assume exchangeable errors, to allow
%                     permutations.
% -pmethod <string> : Method to partition the model. It can be one of:
%                     'Guttman', 'Beckmann' or 'Ridgway'.
%                     Default is 'Beckmann'.
% -corrmod          : Applies FWER-correction of p-values over multiple
%                     modalities.
% -corrcon          : Applies FWER-correction of p-values over multiple
%                     contrasts.
% -saveparametric   : Saves also the parametric p-values.
% -savemask         : Save the effective masks used for each modality.
% -rmethod <string> : Method for regression/permutation. It can be one of:
%                     'Draper-Stoneman','Still-White', 'Freedman-Lane',
%                     'terBraak', 'Kennedy', 'Manly', 'Huh-Jhun' or 'Smith'.
%                     Default is 'Freedman-Lane'.
% -npc              : Use non-parametric combination (NPC).
% -cmethod <string> : Method for combination in the NPC. It can be one
%                     of: 'Tippett', 'Fisher', 'Pearson-David', 'Stouffer',
%                     'Wilkinson', 'Winer', 'Edgington', 'Mudholkar-George',
%                     'Friston', 'Darlington-Hayes', 'Zaykin',
%                     'Dudbridge-Koeleman', 'Dudbridge-Koeleman2',
%                     'Nichols', 'Taylor-Tibshirani', 'Jiang'.
%                     Default is 'Tippett'.
% -draft            : Run in the "draft mode". No FWER correction is
%                     possible, only FDR-adjustment.
% -fdr              : Produces FDR-adjusted p-values.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Sep/2013
% http://brainder.org

% If Octave
if palm_isoctave,
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % If running as a script, take the input arguments
    cmdname = program_invocation_name();
    if ~strcmpi(cmdname(end-5:end),'octave'),
        varargin = argv();
    end
end

% This is probably redundant but fix a bug in an old version
nargin = numel(varargin);

% Print usage if no inputs are given
if nargin == 0 || strcmp(varargin{1},'-q'),
    fprintf('===========================================\n');
    fprintf('PALM: Permutation Analysis of Linear Models\n');
    fprintf('===========================================\n');
    fprintf('\n');
    fprintf('PALM performs permutation inference for the general linear\n');
    fprintf('models (GLMs) of arbitrary complexity, taking as inputs data\n');
    fprintf('in various formats, and being able to take into account certain\n');
    fprintf('cases of well-structured non-independence.\n');
    fprintf('\n');
    fprintf('The options are:\n');
    fprintf('\n');
    fprintf('-i <file>         : Input(s). More than one can be specified.\n');
    fprintf('-m <file>         : Mask(s). Either one for all inputs, or one per input,\n');
    fprintf('                    supplied in the same order as the respective -i appear.\n');
    fprintf('-s <file>         : Surface file(s), for the inputs that are scalars over\n');
    fprintf('                    a surface geometry. Each -s should be supplied in the\n');
    fprintf('                    same order as the respective -i.\n');
    fprintf('-d <file>         : Design matrix. It can be in .csv format, or in FSL''s\n');
    fprintf('                    "VEST" format.\n');
    fprintf('-t <file>         : t-contrasts file, in VEST format.\n');
    fprintf('-f <file>         : F-contrasts file, in VEST format, just as in FSL.\n');
    fprintf('-eb <file>        : Exchangeability blocks file, in VEST format.\n');
    fprintf('-vg <file>        : Variance groups file, in VEST format.\n');
    fprintf('-o <string>       : Output prefix. It may itself be prefixed by a path.\n');
    fprintf('-n <integer>      : Number of permutations. Use -n 0 for exhaustive.\n');
    fprintf('-c <real>         : Theshold for cluster extent for t contrasts.\n');
    fprintf('-C <real>         : Theshold for cluster mass for t contrasts.\n');
    fprintf('-F <positive>     : Theshold for cluster extent for F contrasts.\n');
    fprintf('-S <positive>     : Theshold for cluster mass for F contrasts.\n');
    fprintf('-T                : Enables TFCE inference for 3D (volume) data, i.e.,\n');
    fprintf('                    with H=2, E=0.5, C=6.\n');
    fprintf('-T2               : Enables TFCE inference for 2D (surface, TBSS) data,\n');
    fprintf('                    i.e., H=2, E=1, C=26.\n');
    fprintf('-cnpc <real>      : Threshold for NPC cluster extent (z-stat).\n');
    fprintf('-Cnpc <real>      : Threshold for NPC cluster mass (z-stat).\n');
    fprintf('-Tnpc             : Enables TFCE inference for NPC on 3D data,\n');
    fprintf('                    with H=2, E=0.5, C=6.\n');
    fprintf('-T2npc            : Enables TFCE inference for NPC on 2D data,\n');
    fprintf('                    i.e., H=2, E=1, C=26.\n');
    fprintf('-tfce_H <real>    : Set the TFCE H (height power) parameter.\n');
    fprintf('-tfce_E <real>    : Set the TFCE E (extent power) parameter.\n');
    fprintf('-tfce_C <integer> : Set the TFCE C (connectivity) parameter.\n');
    fprintf('-pb               : Permute blocks as a whole.\n');
    fprintf('-ise              : Assume independent and symmetric errors, to\n');
    fprintf('                    allow sign-flipping.\n');
    fprintf('-ee               : Assume exchangeable errors, to allow\n');
    fprintf('                    permutations.\n');
    fprintf('-pmethod <string> : Method to partition the model. It can be one of:\n');
    fprintf('                    ''Guttman'', ''Beckmann'' or ''Ridgway''.\n');
    fprintf('                    Default is ''Beckmann''.\n');
    fprintf('-corrmod          : Applies FWER-correction of p-values over multiple\n');
    fprintf('                    modalities.\n');
    fprintf('-corrcon          : Applies FWER-correction of p-values over multiple\n');
    fprintf('                    contrasts.\n');
    fprintf('-saveparametric   : Saves also the parametric p-values.\n');
    fprintf('-savemask         : Save the effective masks used for each modality.\n');
    fprintf('-rmethod <string> : Method for regression/permutation. It can be one of:\n');
    fprintf('                    ''Draper-Stoneman'',''Still-White'', ''Freedman-Lane'',\n');
    fprintf('                    ''terBraak'', ''Kennedy'', ''Manly'', ''Huh-Jhun'' or ''Smith''.\n');
    fprintf('                    Default is ''Freedman-Lane''.\n');
    fprintf('-npc              : Use non-parametric combination (NPC).\n');
    fprintf('-cmethod <string> : Method for combination in the NPC. It can be one\n');
    fprintf('                    of: ''Tippett'', ''Fisher'', ''Pearson-David'', ''Stouffer'',\n');
    fprintf('                    ''Wilkinson'', ''Winer'', ''Edgington'', ''Mudholkar-George'',\n');
    fprintf('                    ''Friston'', ''Darlington-Hayes'', ''Zaykin'',\n');
    fprintf('                    ''Dudbridge-Koeleman'', ''Dudbridge-Koeleman2'',\n');
    fprintf('                    ''Nichols'', ''Taylor-Tibshirani'', ''Jiang''.\n');
    fprintf('                    Default is ''Tippett''.\n');
    fprintf('-draft            : Run in the "draft mode". No FWER correction is\n');
    fprintf('                    possible, only FDR-adjustment.\n');
    fprintf('-fdr              : Produces FDR-adjusted p-values.\n');
    fprintf('\n');
    fprintf('_____________________________________\n');
    fprintf('Anderson M. Winkler\n');
    fprintf('FMRIB / Univ. of Oxford\n');
    fprintf('Sep/2013\n');
    fprintf('http://brainder.org\n');
    return;
end

% Now run what matters
palm_backend(varargin);