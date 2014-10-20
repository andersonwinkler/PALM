function palm_help(varargin)
% Shows a help text. Call palm_help('logo') to show just the logo.

if nargin == 0 || strcmpi(varargin{1},'basic'),
    showlogo;
    basic_help;
elseif strcmpi(varargin{1},'advanced'),
    showlogo;
    advanced_help;
elseif strcmpi(varargin{1},'logo'),
    showlogo;
end

% ==============================================================
function basic_help
% Show the most common options.

fprintf('\nThe main options are:\n\n');

fprintf('-i <file> : Input(s). More than one can be specified, each one preceded\n');
fprintf('	by its own -i. All input files must contain the same number of\n');
fprintf('	observations (e.g., the same number of subjects). Except for NPC\n');
fprintf('	and MV, mixing is allowed (e.g., voxelwise, vertexwise and\n');
fprintf('	non-imaging data can be all loaded at once, and later will be all\n');
fprintf('	corrected across).\n\n');

fprintf('-m <file> : Mask(s). Either one for all inputs, or one per input,\n');
fprintf('	supplied in the same order as the respective -i appear.\n\n');

fprintf('-s <file> : Surface file(s). When more than one is supplied, each -s\n');
fprintf('	should be entered in the same order as the respective -i. This\n');
fprintf('	option is needed when the input data is a scalar field over a\n');
fprintf('	surface and spatial statistics (cluster extent, cluster mass or\n');
fprintf('	TFCE) have been enabled.\n\n');

fprintf('-d <file> : Design matrix. It can be in csv format, or in fsl''s vest\n');
fprintf('	format. For information on how to construct the design matrix,\n');
fprintf('	see the FSL GLM manual.\n\n');

fprintf('-t <file> : t-contrasts file, in csv or vest format (the format used by\n');
fprintf('	FSL). The option -t can be used more than once, so that more than\n');
fprintf('	one t-contrasts file can be loaded.\n\n');

fprintf('-f <file> : F-contrasts file, in csv or vest format. The option -f can\n');
fprintf('	be used more than once, so that more than one F-contrasts file\n');
fprintf('	can be loaded. Each file supplied with a -f corresponds to the\n');
fprintf('	file supplied with the option -t in the same order. The option -f\n');
fprintf('	cannot be used more than the number of times the option -t was\n');
fprintf('	used.\n\n');

fprintf('-fonly : Run only the F-contrasts, not the t-contrasts.\n\n');

fprintf('-rmethod <string> : Method for regression/permutation. It can be one of\n');
fprintf('	Freedman-Lane, Smith, terBraak, Manly, Draper-Stoneman,\n');
fprintf('	Still-White and Huh-Jhun. Default, and recommended, is\n');
fprintf('	Freedman-Lane.\n\n');

fprintf('-n <integer> : Number of permutations. Use -n 0 to run all permutations\n');
fprintf('	and/or sign-flips exhaustively. Default is 10000.\n\n');

fprintf('-eb <file> : Exchangeability blocks file, in csv or vest format. If\n');
fprintf('	omitted, all observations are treated as exchangeable and\n');
fprintf('	belonging to a single large exchangeability block.\n\n');

fprintf('-within : If the file supplied with -eb has a single column, this option\n');
fprintf('	runs within-block permutation (default). Can be used with ''-whole''.\n\n');

fprintf('-whole : If the file supplied with -eb has a single column, this option\n');
fprintf('	runs whole-block permutation. Can be used with ''-within''.\n\n');

fprintf('-vg <file> : Variance groups file, in csv or vest format. If omitted,\n');
fprintf('	all observations are assumed to belong to the same variance group (i.e.\n');
fprintf('	the data is treated as homoscedastic. Use ''-vg auto'' to define the\n');
fprintf('	automatically using a structure that is compatible with the\n');
fprintf('	exchangeability blocks (option -eb).\n\n');

fprintf('-ee : Assume exchangeable errors (EE), to allow permutations.\n\n');

fprintf('-ise : Assume independent and symmetric errors (ISE), to allow\n');
fprintf('	sign-flipping.\n\n');

fprintf('-npc <method> : Do non-parametric combination (NPC), using the the\n');
fprintf('	specified method, which can be one of: Tippett, Fisher,\n');
fprintf('	Pearson-David, Stouffer, Wilkinson <alpha>, Winer, Edgington,\n');
fprintf('	Mudholkar-George, Friston, Darlington-Hayes <r>, Zaykin <alpha>,\n');
fprintf('	Dudbridge-Koeleman <r>, Dudbridge-Koeleman2 <r> <alpha>, Nichols,\n');
fprintf('	Taylor-Tibshirani or Jiang <alpha>. Default is Tippett. Note that\n');
fprintf('	some methods require 1 or 2 additional parameters to be provided.\n');
fprintf('	All methods except Darlington-Hayes and Jiang can also be used to\n');
fprintf('	produce parametric p-values and spatial statistics.\n\n');

fprintf('-npcmod : Enable NPC over modalities.\n\n');

fprintf('-npccon : Enable NPC over contrasts.\n\n');

fprintf('-mv <statistic> : Do classical multivariate analysis (MV), such as\n');
fprintf('	MANOVA and MANCOVA, using the the specified statistic, which can\n');
fprintf('	be one of: Wilks, Hotelling, Pillai, Roy_ii or Roy_iii. All but\n');
fprintf('	the last can also be used to calculate parametric p-values and\n');
fprintf('	spatial statistics.\n\n');

fprintf('-o <string> : Output prefix. It may itself be prefixed by a path.\n');
fprintf('	Default is palm.\n\n');

fprintf('-pearson : Instead of t, F, v or G, compute either the Pearson''s\n');
fprintf('	correlation coefficient, r (if the constrast has rank=1), or the\n');
fprintf('	coefficient of determination R2 (if the constrast has rank>1).\n');
fprintf('	For the contrasts in which some EVs are zeroed out, this option\n');
fprintf('	computes the partial correlation (or partial R2).\n\n');

fprintf('-T : Enable TFCE inference for univariate (partial) tests, as well as\n');
fprintf('	for NPC and/or MV if these options have been enabled.\n\n');

fprintf('-c <real> : Enable cluster extent for t contrasts for univariate\n');
fprintf('	(partial) tests, with the supplied cluster-forming threshold\n');
fprintf('	(supplied as the equivalent z-score), as well as for NPC and/or\n');
fprintf('	MV if these options have been enabled.\n\n');

fprintf('-C <real> : Enable cluster mass for t contrasts for univariate\n');
fprintf('	(partial) tests, with the supplied cluster-forming threshold\n');
fprintf('	(supplied as the equivalent z-score), as well as for NPC and/or\n');
fprintf('	MV if these options have been enabled.\n\n');

fprintf('-tfce1D : Set TFCE parameters for 1D data (synchronised timeseries) i.e.,\n');
fprintf('	H = 2, E = 2, C = 6. Use this option together with -T.\n\n');

fprintf('-tfce2D : Set TFCE parameters for 2D data (surface, TBSS)  i.e.,\n');
fprintf('	H = 2, E = 1, C = 26. Use this option together with -T.\n\n');

fprintf('-corrmod : Apply FWER-correction of p-values over multiple modalities.\n\n');

fprintf('-corrcon : Apply FWER-correction of p-values over multiple contrasts.\n\n');

fprintf('-fdr : Produce FDR-adjusted p-values.\n\n');

fprintf('-save1-p : Save (1-p) instead of the actual p-values.\n\n');

fprintf('-logp : Save the output p-values as -log(p) (or -log(1-p) if the option\n');
fprintf('	-save1-p is used).\n\n');

fprintf('-draft : Run in the "draft mode". No NPC, nor FWER correction are\n');
fprintf('	possible. MV and FDR-adjustment are possible.\n\n');

fprintf('-demean : Mean center the data, as well as all columns of the design\n');
fprintf('	matrix. If the original design had an intercept, the intercept is\n');
fprintf('	removed.\n\n');

fprintf('-twotail : Run two-tailed tests for all the t-contrasts instead of\n');
fprintf('	one-tailed. If NPC is used, it naturally becomes two-tailed\n\n');

fprintf('-advanced : Show advanced options.\n\n');

fprintf('_____________________________________\n');
fprintf('Anderson M. Winkler\n');
fprintf('FMRIB / University of Oxford\n');
fprintf('Oct/2014\n');
fprintf('http://brainder.org\n');

% ==============================================================
function advanced_help
% Show advanced options.

fprintf('\nThe advanced or less commonly used options are:\n\n');

fprintf('-inormal : Apply an inverse-normal transformation to the data.\n');
fprintf('	This procedure can go faster if the data is known to be quantitative\n');
fprintf('	(in which case, use ''-inormal quanti''). There are four different\n');
fprintf('	methods available, which can be specified as ''-inormal <method>'' or\n');
fprintf('	''-inormal quanti <method>''. The methods are ''Waerden'' (default),\n');
fprintf('	''Blom'', ''Tukey'' and ''Bliss''.\n\n');

fprintf('-syncperms : If possible, use synchronized permutations even they wouldn''t\n');
fprintf('	be used by default.\n\n');

fprintf('-conskipcount <integer> : Normally the contrasts are numbered from 1, but\n');
fprintf('	this option allows staring the counter from the specified number.\n');
fprintf('	This option doesn''t affect which contrasts are performed.\n\n');

fprintf('-tfce_H <real> : Set the TFCE H parameter (height power).\n\n');

fprintf('-tfce_E <real> : Set the TFCE E parameter (extent power).\n\n');

fprintf('-tfce_C <integer> : Set the TFCE C parameter (connectivity).\n\n');

fprintf('-Tuni : Enable TFCE inference for univariate (partial) tests.\n\n');

fprintf('-cuni <real> : Enable cluster extent for t contrasts for univariate\n');
fprintf('	(partial) tests, with the supplied cluster-forming threshold (as\n');
fprintf('	a z-score).\n\n');

fprintf('-Cuni <real> : Enable cluster mass for t contrasts for univariate\n');
fprintf('	(partial) tests, with the supplied cluster-forming threshold (as\n');
fprintf('	a z-score).\n\n');

fprintf('-Tnpc : Enable TFCE inference for NPC.\n\n');

fprintf('-cnpc <real> : Enable cluster extent for t contrasts for NPC, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-Cnpc <real> : Enable cluster mass for t contrasts for NPC, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-Tmv : Enable TFCE inference for MV.\n\n');

fprintf('-cmv <real> : Enable cluster extent for t contrasts for MV, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-Cmv <real> : Enable cluster mass for t contrasts for MV, with the\n');
fprintf('	supplied cluster-forming threshold (as a z-score).\n\n');

fprintf('-saveparametric : Save also uncorrected parametric p-values.\n\n');

fprintf('-savemask : Save the effective masks used for each modality, as well as\n');
fprintf('	an intersection mask used for NPC and/or MV of these have been\n');
fprintf('	selected.\n\n');

fprintf('-noniiclass : Do not use the NIFTI class (use this option only if you\n');
fprintf('	have problems and even so, only for small datasets).\n\n');

fprintf('-saveperms : Save one statistic image per permutation, as well as a csv\n');
fprintf('	file with the permutation indices (one permutation per row, one\n');
fprintf('	index per column; sign-flips are represented by the negative\n');
fprintf('	indices). Use cautiously, as this may consume large amounts of\n');
fprintf('	disk space.\n\n');

fprintf('-savemetrics : Save a csv file with various permutation metrics.\n\n');

fprintf('-seed <positive> : Seed used for the random number generator. Use a\n');
fprintf('	positive integer up to 2^32. Default is 0. To start with a random\n');
fprintf('	number, use the word ''twist'' instead of the integer. Note that a\n');
fprintf('	given seed used in Matlab isn''t equivalent to the same seed used\n');
fprintf('	in Octave.\n\n');

fprintf('-designperinput : Use one design file for each input modality.\n\n');

fprintf('-ev4vg : Add to the design one EV modelling the mean of each VG.\n\n');

fprintf('-vgdemean : Mean center the data, as well as all columns of the design\n');
fprintf('	matrix, within each VG. Intercepts are removed.\n\n');

fprintf('-pmethodp : Partition method used when defining the set of permutations.\n');
fprintf('	Cab be ''Guttman'', ''Beckmann'', ''Ridgway'' or ''None''.\n');
fprintf('	Default is ''Beckmann''\n\n');

fprintf('-pmethodr : Partition method used during the regression. Valid values\n');
fprintf('	are ''Guttman'', ''Beckmann'', ''Ridgway'' or ''None''.\n\n');
fprintf('	Default is ''Beckmann''\n\n');

fprintf('-cmcp : Ignore the possibility of repeated permutations. This option\n');
fprintf('	will be ignored if the number of shufflings chosen is larger than the\n');
fprintf('	maximum number of possible shufflings, in which case all possible\n');
fprintf('	shufflings will be performed.\n\n');

fprintf('-cmcx : Ignore the possibility of repeated rows in X when\n');
fprintf('	constructing the set of permutations, such that each row is\n');
fprintf('	treated as unique, regardless of potential repetitions (ties).\n\n');

fprintf('-removevgbysize <integer> : Remove from the data and design those\n');
fprintf('	observations that are in VGs of size smaller or equal than specified.\n');
fprintf('	This is specially useful with the option ''-vg auto''.\n\n')

fprintf('-zstat : Convert the original statistic (t, F, v, G, r, R2, or any of\n');
fprintf('	the MV statistics) to a normally distributed, z-statistic.\n\n');

fprintf('-noranktest : For MV, don''t check the rank of the data before trying to\n');
fprintf('	compute the multivariate statistics.\n\n');

fprintf('-evperdat : Treat the design matrix supplied with -d as having one column\n');
fprintf('	for each column in the observed data (entered with -i). This\n');
fprintf('	enables voxelwise/facewise/vertexwise regressors.\n\n');

fprintf('-transposedata : For input data (-i) that are csv tables (2D), transpose\n');
fprintf('	rows and columns.\n\n');

fprintf('-verbosefilenames : Use lengthy filenames, even if the indices go up\n');
fprintf('	to 1 only.\n\n');

fprintf('_____________________________________\n');
fprintf('Anderson M. Winkler\n');
fprintf('FMRIB / University of Oxford\n');
fprintf('Oct/2014\n');
fprintf('http://brainder.org\n');

% ==============================================================
function showlogo
% Show just the logo
fprintf('=======================================================================\n');
fprintf('             ___         ___                         ___\n');
fprintf('            /  /\\       /  /\\                       /__/\\\n');
fprintf('           /  /::\\     /  /::\\                     |  |::\\\n');
fprintf('          /  /:/\\:\\   /  /:/\\:\\    ___     ___     |  |:|:\\\n');
fprintf('         /  /:/~/:/  /  /:/~/::\\  /__/\\   /  /\\  __|__|:|\\:\\\n');
fprintf('        /__/:/ /:/  /__/:/ /:/\\:\\ \\  \\:\\ /  /:/ /__/::::| \\:\\\n');
fprintf('        \\  \\:\\/:/   \\  \\:\\/:/__\\/  \\  \\:\\  /:/  \\  \\:\\~~\\__\\/\n');
fprintf('         \\  \\::/     \\  \\::/        \\  \\:\\/:/    \\  \\:\\\n');
fprintf('          \\  \\:\\      \\  \\:\\         \\  \\::/      \\  \\:\\\n');
fprintf('           \\  \\:\\      \\  \\:\\         \\__\\/        \\  \\:\\\n');
fprintf('            \\__\\/       \\__\\/                       \\__\\/\n');
fprintf('\n');
fprintf('=======================================================================\n');
fprintf('                 Permutation Analysis of Linear Models\n');
fprintf('=======================================================================\n');
