# PALM User Guide

The overall syntax and interface is purposefully similar to that of [randomise](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/randomise.html), so moving from one to another should be straightforward.

## Main options

A description of the main options is below. The same information can be seen if the command `palm` is executed without arguments.

| Option | Description |
| --- | --- |
| `-i <file>` | Input(s). More than one can be specified, each one preceded by its own `-i`. All input files must contain the same number of observations (e.g., the same number of subjects). Except for NPC and MV, mixing is allowed (e.g., voxelwise, vertexwise and non-imaging data can be all loaded at once, and later will be all corrected across). |
| `-m <file>` | Mask(s). Either one for all inputs, or one per input supplied in the same order as the respective `-i` appear. |
| `-s <filesurf> [filearea]` | Surface file(s). When more than one is supplied, each `-s` should be entered in the same order as the respective `-i`. This option is needed when the input data is a scalar field over a surface and cluster extent or TFCE have been enabled. The first argument is the surface file itself. The second is an optional area-per-vertex or area-per-face file, or simply a number. If only the surface file is provided, its area is calculated and used for the computation of spatial statistics (cluster extent and TFCE). If the second argument is given, it should contain the areas, which are then used (e.g., average areas from native geometry after areal interpolation). Alternatively, if the areas are not meaningful for cluster extent or TFCE, this argument can be simply a number, such as "1", which is then used as the area of all vertices or faces. |
| `-d <file>` | Design matrix. It can be in csv format, or in fsl's `vest` format. For information on how to construct the design matrix, see the [FSL GLM manual](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/glm.html). |
| `-t <file>` | t-contrasts file, in `csv` or `vest` format (the format used by FSL). The option `-t` can be used more than once, so that more than one t-contrasts file can be loaded. |
| `-f <file>` | F-contrasts file, in `csv` or `vest` format. The option `-f` can be used more than once, so that more than one F-contrasts file can be loaded. Each file supplied with a `-f` corresponds to the file supplied with the option `-t` immediately before. The option `-f` cannot be used more than the number of times the option `-t` is used. |
| `-fonly` | Run only the F-contrasts, not the t-contrasts. |
| `-n <integer>` | Number of permutations. Use `-n 0` to run all permutations and/or sign-flips exhaustively. Default is 10000. |
| `-eb <file>` | Exchangeability blocks file, in `csv` or `vest` format. If omitted, all observations are treated as exchangeable and belonging to a single large exchangeability block. |
| `-within` | If the file supplied with `-eb` has a single column, this option runs within-block permutation (default). Can be used with `-whole`. |
| `-whole` | If the file supplied with `-eb` has a single column, this option runs whole-block permutation. Can be used with `-within`. |
| `-ee` | Assume exchangeable errors (EE), to allow permutations (more details below). |
| `-ise` | Assume independent and symmetric errors (ISE), to allow sign-flipping (more details below). |
| `-vg <file>` | Variance groups file, in `csv` or `vest` format. If omitted, all observations are assumed to belong to the same variance group (i.e. the data is treated as homoscedastic. Use '-vg auto' to define the automatically using a structure that is compatible with the exchangeability blocks (option `-eb`). |
| `-npcmethod <method>` | Do a modified non-parametric combination (NPC), using the the specified method (combining function), which can be one of: `Tippett`, `Fisher`, `Stouffer`, `Wilkinson <alpha>`, `Winer`, `Edgington`, `Mudholkar-George`, `Friston <u>`, `Darlington-Hayes <r>`, `Zaykin <alpha>`, `Dudbridge-Koeleman <r>`, `Dudbridge-Koeleman2 <r> <alpha>`, `Taylor-Tibshirani` or `Jiang <alpha>`. Default is `Fisher`. Note that some methods require 1 or 2 additional parameters to be provided. All methods except `Darlington-Hayes` and `Jiang` can also be used to produce parametric p-values (under certain assumptions) and spatial statistics. |
| `-npcmod` | Enable NPC over modalities. |
| `-npccon` | Enable NPC over contrasts. |
| `-npc` | Shortcut to `-npcmethod <default method> -npcmod`. The default method is `Fisher` (this can be changed in the file `palm_defaults.m`). |
| `-mv <statistic> <k>` | Do classical multivariate analysis (MV), such as MANOVA and MANCOVA, using the the specified statistic, which can be one of: `Wilks`, `HotellingTsq`, `Lawley`, `Pillai`, `Roy_ii`, `Roy_iii`. All but `Roy_iii` can be used with spatial statistics. |
| `-pearson` | Instead of t, F, v or G, compute the Pearson's correlation coefficient, r (if the constrast has rank=1), or the coefficient of determination, R2 (if the constrast has rank>1). For the contrasts in which some EVs are zeroed out, this option computes the multiple correlation coefficient (or R2) corresponding to the EVs of interest. |
| `-T` | Enable TFCE inference for univariate (partial) tests, as well as for NPC and/or MV if these options have been enabled. |
| `-C <real>` | Enable cluster inference for univariate (partial) tests, with the supplied cluster-forming threshold (supplied as the equivalent z-score), as well as for NPC and/or MV if these options have been enabled. |
| `-Cstat <name>` | Choose which cluster statistic should be used. Accepted statistics are `extent` and `mass`. |
| `-tfce1D` | Set TFCE parameters for 1D data (timeseries), i.e., H = 2, E = 2, C = 6. Use this option together with `-T`. |
| `-tfce2D` | Set TFCE parameters for 2D data (surface, TBSS), i.e., H = 2, E = 1, C = 26. Use this option together with `-T`. |
| `-corrmod` | Apply FWER-correction of p-values over multiple modalities. |
| `-corrcon` | Apply FWER-correction of p-values over multiple contrasts. |
| `-fdr` | Produce FDR-adjusted p-values. |
| `-o <string>` | Output prefix. It may itself be prefixed by a path. Default is `palm`. |
| `-save1-p` | Save (1-p) instead of the actual p-values. Instead of `-save1-p`, consider using `-logp`. |
| `-logp` | Save the output p-values as \(-log10(p)\) (or \(-log10(1-p)\) if the option `-save1-p` is also used; using both together is not recommended). The `-logp` is not default but it is strongly recommended. |
| `-demean` | Mean center the data, as well as all columns of the design matrix. If the original design had an intercept, the intercept is removed. |
| `-twotail` | Run two-tailed tests for all the t-contrasts instead of one-tailed. If NPC is used, it also becomes two-tailed for the methods which statistic are symmetric around zero under the null. |
| `-concordant` | For the NPC, favour alternative hypotheses with concordant signs. Cannot be used with `-twotail`. |
| `-accel <method>` | Run one of various acceleration methods (`negbin`, `tail`, `noperm`, `gamma`, `lowrank`). Click [here](faster_inference.md) for details. |
| `-reversemasks` | Reverse 0/1 in the masks, so that the zero values are then used to select the voxels/vertices/faces. |
| `-quiet` | Don't show progress as the shufflings are performed. |
| `-advanced` | Show advanced options (see below). |

Advanced options are:

| Option | Description |
| --- | --- |
| `-con <file1> <file2>` | Contrast file(s) in `.mset` format. For hypotheses of the form `H0: C'*Beta*D`, `file1` contains a set of C contrasts, and `file2` (optional) contains a set of D contrasts, that are performed in pairs (C,D). |
| `-tonly` | Run only the t-contrasts, not the F-contrasts. |
| `-cmcp` | Ignore the possibility of repeated permutations. This option will be ignored if the number of shufflings chosen is larger than the maximum number of possible shufflings, in which case all possible shufflings will be performed. |
| `-cmcx` | Ignore the possibility of repeated rows in X when constructing the set of permutations, such that each row is treated as unique, regardless of potential repetitions (ties). |
| `-conskipcount <integer>` | Normally the contrasts are numbered from 1, but this option allows staring the counter from the specified number. This option doesn't affect which contrasts are performed. |
| `-Cuni <real>` | Enable cluster inference for univariate (partial) tests, with the supplied cluster-forming threshold (as a z-score). |
| `-Cnpc <real>` | Enable cluster inference for NPC, with the supplied cluster-forming threshold (as a z-score). |
| `-Cmv <real>` | Enable cluster inference for MV, with the supplied cluster-forming threshold (as a z-score). |
| `-designperinput` | Use one design file for each input modality. |
| `-ev4vg` | Add to the design one EV modelling the mean of each VG. |
| `-evperdat <string> [integer] [integer]` | Treat the design matrix supplied with `-d` as having one column for each column in the observed data (entered with `-i`). This enables voxelwise/facewise/vertexwise regressors. For details, click [here](masks_and_voxelwise_regressors.md). |
| `-inormal` | Apply an inverse-normal transformation to the data. This procedure can go faster if the data is known to be quantitative (in which case, use `-inormal quanti`). There are four different methods available, which can be specified as `-inormal <method>` or `-inormal quanti <method>`. The methods are `Waerden` (default), `Blom`, `Tukey` and `Bliss`. |
| `-probit` | Apply a probit transformation to the data. |
| `-inputmv` | Treat the (sole) input as multivariate, that is, each column is a variable in a multivariate model, as opposed to independent univariate tests. Useful with non-imaging data. |
| `-noranktest` | For MV, don't check the rank of the data before trying to compute the multivariate statistics. |
| `-noniiclass` | Do not use the NIFTI class (use this option only if you have problems and even so, only for small datasets). |
| `-precision <string>` | Precision (`single` or `double`) for input files. |
| `-nounivariate` | Don't save univariate results. |
| `-nouncorrected` | Don't save uncorrected results. |
| `-pmethodp` | Partition method used when defining the set of permutations. Can be `Guttman`, `Beckmann`, `Ridgway` or `None`. Default is `Beckmann`. |
| `-pmethodr` | Partition method used during the regression. Valid values are `Guttman`, `Beckmann`, `Ridgway` or `None`. Default is `Beckmann`. |
| `-removevgbysize <integer>` | Remove from the data and design those observations that are in VGs of size smaller or equal than specified. This is especially useful with the option `-vg auto`. |
| `-rmethod <string>` | Method for regression/permutation. It can be one of `Freedman-Lane`, `Dekker`, `terBraak`, `Manly`, `Draper-Stoneman`, `Still-White` and `Huh-Jhun`. Default, and recommended, is `Freedman-Lane`. |
| `-savedof` | Save file with the degrees of freedom. |
| `-savemask` | Save the effective masks used for each modality, as well as an intersection mask used for NPC and/or MV of these have been selected. |
| `-savemetrics` | Save a `csv` file with various permutation metrics. |
| `-saveparametric` | Save also uncorrected parametric p-values. These are only valid if all assumptions are met, including iid and normality. |
| `-saveglm` | Save COPEs and VARCOPEs in the first permutation. |
| `-saveperms` | Save one statistic image per permutation, as well as a `csv` file with the permutation indices (one permutation per row, one index per column; sign-flips are represented by the negative indices). Use cautiously, as this may consume large amounts of disk space. |
| `-savemax` | Save the distribution of the maximum statistic across (taken within each input file). |
| `-nounivariate` | For classical multivariate and NPC, don't save the univariate results. |
| `-seed <positive>` | Seed used for the random number generator. Use a positive integer up to \(2^32\). Default is 0. To start with a random number, use the word `twist` instead of the integer. Note that a given seed used in Matlab isn't equivalent to the same seed used in Octave. |
| `-syncperms` | If possible, use synchronized permutations even they wouldn't be used by default. |
| `-subjidx <file>` | Indices of input data to keep in the design. |
| `-tfce_H <real>` | Set the TFCE H parameter (height power). |
| `-tfce_E <real>` | Set the TFCE E parameter (extent power). |
| `-tfce_C <integer>` | Set the TFCE C parameter (connectivity). |
| `-tfce_dh <real>` | Set the "delta h" for the calculation of TFCE. Default is `auto` (or 0), i.e., defined on a per-case basis, with `max(:)/100` for each map. |
| `-Tuni` | Enable TFCE inference for univariate (partial) tests. |
| `-Tnpc` | Enable TFCE inference for NPC. |
| `-Tmv` | Enable TFCE inference for MV. |
| `-transposedata` | For input data (`-i`) that are csv tables, transpose rows and columns. If `-evperdat` is used, these are also transposed. |
| `-verbosefilenames` | Use lengthy filenames, even if the numbering goes up to 1 only. |
| `-vgdemean` | Mean center the data, as well as all columns of the design matrix, within each VG. Intercepts are removed. |
| `-zstat` | Convert the original statistic (t, F, v, G, r, R2, or any of the NPC or MV statistics) to a normally distributed, z-statistic. |

### Input files

Input files are recognised by the file extension:

| Extension | Read as |
| --- | --- |
| `.nii`, `.hdr`, `.img` | NIFTI file. Can be used to input data (option `-i`) and mask (`-m`). Note that the old ANALYZE format is not supported. |
| `.nii.gz` | Compressed NIFTI file. Should be uncompressed manually (with `gunzip`) and loaded as `nii` (it uses more disk space, but far less memory). Or can be read directly with the option `-noniiclass`, but beware that for large datasets, the data may easily occupy 'all' the computer memory. |
| `.mgh`, `.mgz` | FreeSurfer data, either volumetric or surface-based. Can be used to input data (option `-i`) or mask for the same kind of data (option `-m`). |
| `.dpv`, `.dpf`, `.dpx` | Data-per-vertex (vertexwise), data-per-face (facewise), or unspecified surface-based data. Can be used to specify a mask (option `-m`) for surface-based data. Multiple such files can be merged into a `csv` table, and then input with the option `-i`; see details below. |
| `.csv` | Table in `csv` format. Can be used to specify input data (option `-i`), mask for the same kind of data (option `-m`), design matrix (option `-d`), contrasts (options `-t` and `-f`), exchangeability blocks (`-eb`) and variance groups (`-vg`). |
| `.mat`, `.con`, `.fts`, `.grp` | Table in the FSL `vest` format. Can be used essentially in the same way as `csv` files, i.e., to specify input data (option `-i`), mask for the same kind of data (`-m`) design matrix (`-d`), contrasts (options `-t` and `-f`), exchangeability blocks (`-eb`) and variance groups (`-vg`). |
| `.mset` | Multiple tables (arrays) in a single ASCII file. This format can be used with the option `-con`. |
| `.srf` | Surface in ASCII format (option `-s`). |
| `.inflated`, `.nofix`, `.orig`, `.pial`, `.smoothwm`, `.sphere`, `.reg`, `.white` | Files with these extensions are read as FreeSurfer surface (option `-s`). |
| `.gii` | GIFTI file. |
| `.dtseries.nii`, `.ptseries.nii`, `.dscalar.nii`, `.pscalar.nii` | CIFTI file. Support is currently available for `.dscalar.nii`, `.dtseries.nii`, `.pscalar.nii` and `.ptseries.nii.` It is expected that in the future there will be complete support for CIFTI files. |

Files with extensions not listed above won't be read.

#### Support for NIFTI files

Support for NIFTI files is provided, internally, by the publicly available NIFTI class. This allows reading and writing even huge files without using too much computer memory. However, the NIFTI class does not operate on compressed files, i.e., with extension `.nii.gz.` To read these files, it is recommended that they are uncompressed (with `gunzip`).

Alternatively, if the datasets are small, and if either FSL or FreeSurfer are installed, the NIFTI class can be disabled with the option `-noniiclass`. This allows to read/write these `.nii.gz` files directly. However, if the files are too large, this can easily use all the computer memory and the system may become unstable/unusable. The option `-noniiclass` should be used with great caution for large datasets.

The NIFTI class, that is used by default, already comes with precompiled binaries for Matlab in various platforms, and for Octave for most 64-bit Linux distributions. Nonetheless, if compilation is needed, try using:

```
cd /full/path/to/palm/fileio/@file_array/private
./compile.sh
```

#### Support for FreeSurfer files

To be able to read/write FreeSurfer binary files (surfaces and curvatures), FreeSurfer needs to be installed and configured correctly. Surface files in ASCII format are read directly by PALM as long as their file extension is srf, so FreeSurfer isn't necessary for these files.

Freesurfer "curvature" files converted to ASCII (`asc`/`dpv`/`dpf`/`dpx`) need to be merged and converted to `csv` tables. To generate a valid `csv` file, use the command `dpx2csv` (available [here](https://raw.githubusercontent.com/andersonwinkler/toolbox/master/bin/dpx2csv)) to make the `csv` file, initially with one column per subject and one row per vertex or face, then transpose the rows and columns of this file with the command `transpose` (available [here](https://raw.githubusercontent.com/andersonwinkler/toolbox/master/bin/transpose)).

#### Support for MZ3 files

MZ3 is a highly efficient, fast and compressed format to store triangular surfaces developed by Chris Rorden and colleagues. A full description of the format is available [here](https://github.com/neurolabusc/surf-ice/tree/master/mz3). Support is provided natively.

#### Support for CSV files

Although PALM was created with imaging in mind, it can be used for non-imaging data. Any kind of data that can be arranged in a table and saved as comma-separated values file (`csv`) (i.e., any data) can be used as input. Each row constitutes a subject or observation, and each column represents a measurement. The `csv` files must contain numeric fields only, i.e., title labels are not accepted and will cause errors.

If the data is arranged as space- or tab-separated values instead of comma-separated values, it can be converted quickly using `awk` (or `gawk`):

```
awk 'BEGIN { OFS="," } { $1=$1; print $0 }' oldfile.txt > newfile.csv
```

If your data has any other separator, e.g., a semicolon, the same applies with using a small modification:

```
awk 'BEGIN { FS="x"; OFS="," } { $1=$1; print $0 }' oldfile.txt > newfile.csv
```

where the `"x"` should be replaced by the separator in the original table (e.g., a semicolon).

#### Support for MSET files

To allow multivariate contrasts, PALM uses a simple ASCII format that can contain multiple matrices. An example of such file is:

```
Matrix 1 3
1 -1 0

Matrix 2 3
1 -1 0
1 0 -1
```

In the example, the first matrix is defined as having 1 row and 3 columns (the numbers after the keyword `Matrix`). The second matrix is defined as having 2 rows and 3 columns. Files as these are meant to be used to test multivariate hypotheses as `H: C'*Beta*D` using the option `-con <file1> <file2>`, where `<file1>` is a file with multiple contrasts C and `<file2>` with multiple contrasts D. Each contrast in C pairs with a contrast in D.

#### Support for CIFTI files

Both the surface and volume components in CIFTI files that are of the type `dtseries`, `ptseries`, `dscalar` and `pscalar` can be read directly. If no spatial statistics are requested (i.e., no TFCE or cluster-level inference), all “grayordinates” (whether surface or volume) are processed without any complications. However, when spatial statistics are requested, it is currently necessary to first manually split the CIFTI files into separate surface (GIFTI) and volume (NIFTI) components, which then can be loaded and processed (see examples [in this page](examples.md)). The [Connectome Workbench](http://www.humanconnectome.org/software/connectome-workbench.html) must be installed. CIFTI can be used with PALM in Octave (both Mac and Linux) and in Matlab (Linux only).

#### Support for GIFTI files

GIFTI files are supported. However, note that the access to the actual data depends on parsing a potentially large XML tree, which can be slow in both Matlab and Octave. Other equivalent formats are usually faster to load.

#### Support for other file formats

If there is a file format that cannot be easily converted to `csv`, and which you'd like to be able to read directly in PALM, feel free to contact us. New formats may be added for future releases.

### Permutations or sign-flips?

Exchangeable errors (EE) are assumed by default, such that permutations are performed. If the option `-ise` is supplied, the errors are instead assumed to be independent and symmetric (ISE) and only sign-flips are performed. If both `-ee` and `-ise` options are given, the errors are treated as being exchangeable, independent and symmetric, and both permutations and sign-flips are performed.

### Defining exchangeability blocks

A basic definition of exchangeability blocks (EBs) can be made using a text file containing a single column of positive integers representing the blocks, i.e., one value per line. This file is supplied with the option `-eb`, and indicates that the observations corresponding to each line can only be shuffled within block of the same index. If a file with this simple structure is given with `-eb` and the option `-whole` is used, then the blocks are permuted as a whole, with the order of observations inside each block being kept unchanged. If the option `-ise` is used, sign-flips ignore block definitions unless the option `-whole` is given, in which case the signs of the blocks as a whole are flipped. If the option `-whole` is used together with the option `-within`, permutations happen for the blocks as a whole, and further within each block.

While the simple block specification as above, with just one column of indices, can cover various designs with structured dependence between observations, cases of more complex dependence can be accommodated using multiple, nested block definitions. In this case, the file containing the block definitions needs not just one, but multiple columns, each with a sequentially deeper level of dependence. Indices on one level indicate how the indices of the next level should be permuted. Positive indices at one level indicate that the corresponding indices of the next level should be permuted as a whole, akin to whole-block permutation with the `-whole` option described above. Negative indices at one level indicate that the corresponding indices in the next level should not be permutted, and their own subindices should be permuted within block, akin to within-block permutation described above. With this more complex structure, the option `-whole`, even if supplied, is ignored, as this information is embedded, at multiple levels, in the file with the block definitions.

For more details and examples, [click here](exchangeability_blocks.md).

### Output files

All output file names are prefixed by the string supplied with the option `-o`. If this option is not supplied, the outputs will be saved to the directory from which palm is executed, all named as `palm_*`. The file names follow a consistent convention described below.

The file format of the outputs is the same as the input data (supplied with the option `-i`, collapsed along the dimension of the observations, i.e., the 4th dimension for 4D datasets (NIFTI or FreeSurfer surface formats), or the 1st dimension (rows) of `csv` files or 2D NIFTI files.

There are 2 outputs of general interest: the statistic image and its associated p-value. If the option `-saveperms` is used, the statistic computed for each permutation is also saved, as well as a table with the permutation indices for these images. These files are named according to the following convention:

| File name | File contents |
| --- | --- |
| `{prefix}_{unit}_{stat}_m{#}_d{#}_c{#}.{ext}` | Statistic for the modality `m{#}`, for the design `d{#}` and for contrast `c{#}`. |
| `{prefix}_{unit}_{stat}_{pval}_m{#}_d{#}_c{#}.{ext}` | P-value for the modality `m{#}`, for the design `d{#}` and for contrast `c{#}`. |
| `{prefix}_{unit}_{stat}_m{#}_d{#}_c{#}.dof` | A text file containing the degrees of freedom for the statistic computed for the modality `m{#}`, for the design `d{#}` and for contrast `c{#}`. This is only useful to refer to the parametric distribution to find the p-values (which can, nonetheless, be produced directly with the option `-saveparametric`. For the v and G statistics, the dof of the errors are saved in the same format as the inputs. |
| `{prefix}_{unit}_{stat}_m{#}_d{#}_c{#}_perm{#}.{ext}` | Statistic for the permutation perm{#} for the modality `m{#}`, for the design `d{#}` and for contrast `c{#}`, if the option `-saveperms` is used. |
| `{prefix}_{unit}_{stat}_m{#}_d{#}_c{#}_permidx.csv` | Permutation indices, with negatives for the sign flips, if the option `-saveperms` is used. This `csv` file has one row per permutation and one column per permutation index. |

Each of these fields are:

| Field | Description |
| --- | --- |
| `{prefix}` | The prefix supplied with the option `-o`, which can contain full path information. If `-o` isn't supplied, the default prefix is `palm`. |
| `{unit}` | The unit in space in which the statistic is calculated. This can be `vox` for voxelwise, `vtx` for vertexwise, `fac` for facewise, `dpx` for surface-data without specification (it can be vertexwise or facewise if the actual surface geometry wasn't supplied with the `-s` option), `dat` for any data that doesn't necessarily have a spatial meaning (e.g., data entered in csv format), `clustere` for cluster extent, `clusterm` for cluster mass, or `tfce` for TFCE. |
| `{stat}` | The name of the statistic, which can be `tstat` for the t-statistic, `fstat` for the F-statistic, `vstat` for the Aspin–Welch's v-statistic, `gstat` for the generalised G-statistic, `rstat` for the Person's r correlation coefficient, `rsqstat` for the R2 coefficient of determination. For NPC, it is `npc`, appended by the name of the combining function. For MV, it is `mv` appended by the name of the statistic, or by `tsqstat` for Hotelling's T2. If the option `-zstat` was enabled, all these names are prefixed by the letter `z`. |

|`{pval}`|The type of p-value stored in the image. It can be `uncp` for uncorrected permutation p-values, `fwep` for the p-values FWER-corrected within modality and contrast, `cfwep` for the p-values FWER-corrected across contrasts if the option `-corrcon` was used (so, within modality), `mfwep` for the p-values FWER-corrected across modalities if the option `-corrmod` was used (so, within contrast), `mcfwep` for the p-values FWER-corrected across modalities and contrasts, `fdrp` for p-values FDR-adjusted within modality and contrast, `uncparap` for uncorrected parametric p-value and `fdrparap` for parametric FDR-corrected p-values.|

|`m{#}`|Modality index, in the same order as supplied as input with the option `-i`. This is omitted for the NPC over modalities and for MV statistics, or if there is just one modality.|

|`d{#}`|Design index, in the same order as supplied with the option `-d`. This is omitted if there is just one design.|
|`c{#}`|Contrast index, in the same order as in the files supplied with the `-t` and/or `-f` options. This is omitted if there is just one contrast.|
|`perm{#}`|Statistic map for the permutation `perm{#}` if the option `-saveperms` is used.|
|`{ext}`|File extension. This will typically be the same extension as the corresponding input files. |

In addition, other files are created: a small text file is saved containing the degrees of freedom associated with the statistic. For the v and G statistics, the degrees of freedom associated with the errors vary for each image point and are saved as a map. A configuration file with all the input options is also saved. This file can be used as the sole argument to palm to run the same analysis again.

### PALM or randomise?

Compared to randomise, PALM performs overall the same type of test, i.e., permutation, but offering various additional or novel features:

| Feature | Randomise | PALM |
| --- | --- | --- |
| Input formats | NIFTI. | NIFTI, CSV, GIFTI, FreeSurfer (surface formats, MGH files, curvature, including ASCII vertexwise and facewise). Partial support for CIFTI. |
| Univariate statistics | t and F. | t, F, Aspin–Welch's v, G, Pearson's r and R2. |
| Multivariate statistics | None | Hotelling's T2, Wilks' lambda, Pillai's trace, Lawley–Hotelling's trace, Roy's largest root (both variants), assessed parametrically and non-parametrically. Details [here](joint_inference.md). |
| Spatial statistics | Cluster extent, cluster mass, TFCE. All of these can be volume-based (voxelwise). | Cluster extent, cluster mass, TFCE. All of these can be volume-based (voxelwise) or surface-based (vertexwise or facewise). |
| Regression and permutation strategies | Freedman–Lane. | Draper–Stoneman, Still–White, Freedman–Lane, ter Braak, Manly, Huh–Jhun, Dekker and parametric (no permutation). |
| Shuffling possibilities | Permutations or sign-flippings. | Permutations, sign-flippings, and permutations with sign-flippings. |
| Block permutation | Either whole-block or within-block. | Whole-block and/or within-block, either in isolation or in arbitrary combinations, with multiple and complex levels of tree-like dependence between observations. Details [here](exchangeability_blocks.md). |
| Variance groups | One. Cannot be changed. | Multiple. Can be set automatically or defined by the user. Details [here](exchangeability_blocks.md). |
| Input datasets (modalities) | One. | Multiple. |
| Input designs matrices | One. | Multiple. |
| Joint permutation inference across modalities and/or contrasts (NPC) | No. | Yes. Various combining functions available: Tippett, Fisher, Pearson–David, Stouffer, Wilkinson, Winer, Edgington, Mudholkar–George, Friston, Darlington–Hayes, Zaykin, Dudbridge–Koeleman (2 variants), Taylor–Tibshirani and Jiang. Most of these can also be assessed parametrically. Details [here](joint_inference.md). |
| FWER correction | Within each contrast only. | Within contrast, across contrasts, within modality, across modalities, and across modalities and contrasts. For correction across modalities, these don't need to have the same resolution or geometry; in fact, they don't need to be imaging data (and images can be mixed with non-images). |
| FDR correction | Yes, separate command. Correction within each contrast only. | Yes, built-in. Correction can be within contrast, across contrasts, within modality, across modalities, and across modalities and contrasts. For correction across modalities, these don't need to have the same resolution or geometry; in fact, they don't need to be imaging data (and images can be mixed with non-images). |
| Non-stationarity correction | Yes. | No. |
| Voxelwise regressors | Yes. | Yes. Details [here](masks_and_voxelwise_regressors.md). |
| Fast approximations | No. | Yes, various methods described here. |
| Parametric p-values | No. | Yes, for univariate, multivariate, and most of the NPC combining functions (uncorrected and FDR-corrected). |
| Permutation metrics | None. | Average Hamming distance, Hubermann–Hogg complexity and entropy measurements. |
| Overall performance | Generally a bit slower, except for spatial statistics. | Generally a bit faster, except for spatial statistics, that can be extremely slow. Includes acceleration methods (details [here](faster_inference.md)). |
| Parallel processing | Yes, with a separate command (`randomise_parallel`). | No. |
| Language | Written in C++. | Written in Matlab/Octave. |
| License | [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html) (Free for academic use). | Free ([GPL](http://www.gnu.org/licenses/gpl.html)). |