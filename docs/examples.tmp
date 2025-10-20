# Examples

Some simple examples to get started are below. Note that by default, PALM saves the actual p-values, rather than 1-p. To have 1-p, use the option `-save1-p`. Alternatively, use the option `-savelogp`, which is good for visualisation.

## Example 1: Single modality

To run PALM for a single modality, use use something as:

```
palm -i dataset4D.nii -d design.mat -t design.con -m mask.nii -f design.fts -T -C 3.1 -n 5000 -corrcon -o myresults
```

Note that for these simple cases, the calls are nearly identical to what would be used in randomise. There are a few differences, though: the cluster forming threshold is supplied as the equivalent z-score for the t or F statistics (or whatever other statistic that will be computed internally, such as v, G, r or R2), and the value given with the option `-c` or `-C` applies to all these statistics (not just t). Also, while .nii.gz can be loaded directly, the current version prefers `.nii` (uncompressed), because it's faster to load but mainly, because for large datasets it doesn't risk occupying all the memory and making the system unstable.

The `design.mat`, `design.con` and `design.fts` can be replaced for `csv` files. The extension for these files needs to be `.csv` then. Likewise, the input doesn't need to be an image. It can be any `csv` file in which each row is a subject and each column is a measurement.

The option `-corrcon` does FWER-correction across all contrasts, taking into account any dependency that may exist between these contrasts.

If the design and contrast files are both omitted, a 1-sample t-test will be performed, with sign-flipping. If only the contrast files are omitted, and if there is a single EV in the design, it will be a 1-sample t-test; if there are more than one EV, it will run an F-test testing whether any EV has its regression coefficient different than zero.

If `-n 0` is used, the permutations are performed exhaustively. If the number of permutations is omitted, the default is 10000.

The option `-o` is to specify the output prefix (and it may include a path). If not given, the default output prefix is simply `palm`.

## Example 2: Permutations with sign-flipping

The example above, as shown, will be executed with permutations by default. To replace permutations for sign-flippings, the option `-ise` is included (this specifies that the errors are independent and symmetric, ISE)

```
palm -i dataset4D.nii -d design.mat -t design.con -m mask.nii -f design.fts -T -C 3.1 -n 5000 -corrcon -o myresults -ise
```

To run permutations together with sign-flippings, use the option `-ise` together with the option `-ee`, meaning that the errors are exchangeable (EE) and also independent and symmetric (ISE)

```
palm -i dataset4D.nii -d design.mat -t design.con -m mask.nii -f design.fts -T -C 3.1 -n 5000 -corrcon -o myresults -ee -ise
```

## Example 3: Multiple modalities, corrected across

To run PALM for multiple modalities, use something as:

```
palm -i modality1_4d.nii -i modality2_4d.nii -i modality3_4d.nii -d design.mat -t design.con -n 2000 -corrcon -corrmod
```

Each modality will be tested separately, and the option -corrmod will do FWER-correction across all of them, taking into account the dependence structure that may exist between them. These modalities don't need to be in the same resolution, don't need to be imaging, and can mix voxel-based with surface-based measurements.

## Example 4: Multiple designs

To run PALM with multiple designs use:

```
palm -i data_4d.nii -d design1.mat -t design1.con -d design2.mat -t design2.con -d design3.mat -t design3.con [...] -n 2000 -corrcon
```

Each contrast for each design will be tested separately, and the option -corrcon will do FWER-correction across all of them, taking into account the dependence structure that may exist between them.

## Example 5: Using multi-level block definitions

The file with the multi-level blocks is supplied with the option `-eb`. If the variances cannot be assumed to be the same for all observations, a file with variance groups can be supplied with the option `-vg`, or a conservative estimate for the groups can be defined, from the blocks, with `-vg auto`. If the option `-vg` is not used, all data is assumed as having the same variance.

If the file supplied with `-eb` has just one column, the options `-whole` and/or `-within` can be used to indicate how shuffling should happen. If he file has multiple columns, the information on how the blocks should be shuffled is extracted from the blocks themselves. More details [here](exchangeability_blocks.md).

```
palm -i modality1_4d.nii -d design.mat -t design.con -eb blocks.csv
```

## Example 6: Multiple modalities, joint inference with NPC

If the modalities have the same spatial resolution or geometry (for surfaces), it's possible to combine them and perform joint inference involving all of them. There are two methods available: one is NPC (Non-Parametric Combination), the other is MV (classic MANOVA/MANCOVA; next example).

For the NPC, the call is something as this:

```
palm -i modality1_4d.nii -i modality2_4d.nii -i modality3_4d.nii -d design.mat -t design.con -f design.fts -n 5000 -npc
```

This will combine the modalities using the Fisher method (default). There are various others available, and Fisher is probably the best, and it's also quite fast. The significance after the combination is assessed via permutations. More details [here](joint_inference.md).

## Example 7: MANOVA and MANCOVA (simple)

For the MV, the call is something as this:

```
palm -i modality1_4d.nii -i modality2_4d.nii -i modality3_4d.nii -d design.mat -t design.con -f design.fts -n 5000 -mv
```

For the F-contrasts, this runs the equivalent to a MANOVA/MANCOVA, computing the Wilks' Lambda (default for F-tests), and assessing its significance through permutations. Other multivariate statistics are also available (Roy, Lawley-Hotelling, Pillai). For the t-contrasts, the default Hotelling's \(T^2\). More details [here](joint_inference.md).

## Example 8: MANOVA and MANCOVA (with contrasts of response variables)

It is possible to further refine MANOVA/MANCOVA by providing contrasts of response variables. The null hypothesis is `H0 : C'*psi*D`, where `psi` are the estimated regression coefficients, `C` are contrasts as usual, and `D` are contrasts of response variables. These are supplied with the option `-con`:

```
palm -i modality1_4d.nii -i modality2_4d.nii -i modality3_4d.nii -d design.mat -con C.mset D.mset -n 5000 -mv
```

Each `*.mset` file contains multiple C and D contrasts, pairwise, that is, if `C.mset` contains 4 contrasts, the `D.mset` needs to contain also 4 contrasts, so that 4 different null hypotheses are tested, the first using the first contrast of each file, the second using the second of each file, and so on. More details [here](joint_inference.md) and [here](user_guide.md).

## Example 9: Faster inference

The most general acceleration method, that can be used in a wide variety of situations, is the tail approximation. A typical use is for FWER-corrected outputs, such that the uncorrected can be skipped. In the example below, 500 permutations are used. More details are [here](faster_inference.md).

```
palm [...] -n 500 -accel tail -nouncorrected
```

## Example 10: Using CIFTI files

**Case 1: Without spatial statistics (i.e., no TFCE or cluster-level inference):** If no spatial statistics will be used, i.e., no TFCE and no cluster-level inference, then CIFTI files are supported as any other input file, i.e., same as with NIFTI, CSV, FreeSurfer, or any other format.

Note that you can use `wb_shortcuts -cifti-concatenate` (from the Workbench tools) to conveniently merge/concatenate CIFTI files (e.g., `dtseries.nii`, `dscalar.nii`, `ptseries.nii`, and `pscalar.nii`) from individual subjects into a single combined file to serve as input to PALM. Similarly, `wb_shortcuts -metric-concatenate` can be used if you need to merge `func.gii` files.

**Case 2: With spatial statistics (i.e., TFCE or cluster-level inference):** Inference using TFCE and cluster-inference with dense CIFTI files that include both cortical surfaces and subcortical volumes are a bit more elaborate. These spatial statistics depend on the size, shape, and topology of the particular type of representation of the brain (surface or volume), such that a correction across these structures based on the distribution of the maximum statistic would lead to very conservative results for some of them. Instead, the CIFTI file needs to be broken down into its components, i.e., surface-based (GIFTI) and volume-based (NIFTI) structures, and PALM is then run separately for each, with the results only then corrected for the multiplicity of tests using either Bonferroni or Šidák methods. Moreover, for the analysis of the cortex, it is also necessary to supply a mesh indicating how the vertices are connected to each other, and optionally (but ideally) a file with the average area per vertex should also be provided.

To separate a CIFTI file into GIFTI (with the cortex of each cerebral hemisphere) and NIFTI (with subcortical structures) components, use the option `-cifti-separate` of the `wb_command` (this assumes that the `dscalar.nii` file has already been concatenated across subjects, using e.g., `wb_shortcuts -cifti-concatenate`):

```
wb_command -cifti-separate data.dscalar.nii COLUMN -volume-all data_sub.nii -metric CORTEX_LEFT data_L.func.gii -metric CORTEX_RIGHT data_R.func.gii
```

The above will produce three files that are used below as inputs in separate calls to PALM: `data_L.func.gii`, `data_R.func.gii`, and `data_sub.nii.` The two `.gii` files (GIFTI) are internally compressed (`GZIP_BASE64_BINARY`), but currently Octave (as of version 4.0.0) cannot do the uncompression. Thus we need to convert to uncompressed `BASE64_BINARY.` The `wb_command` has the option `-gifti-convert` that can be used for conversions as these (do it for both L and R):

```
wb_command -gifti-convert BASE64_BINARY data_L.func.gii data_L.func.gii
```

The surfaces (meshes) for the left and right are also necessary, so that the adjacency (neighbourhood) information between vertices can be computed. Use preferably a *midthickness*, *white*, or *pial* surface. Although the sphere can be used (the adjacency information is the same), it has severe areal distortions that make the use of average areas as input mandatory. The average areas are optional if the midthickness surfaces are used, but still recommended for better representing the sizes of clusters and of the support regions of TFCE. To obtain such areas, use the option `-surface-vertex-areas` of the `wb_command` for each subject, then merge all into a single file with the option `-metric-merge`, then finally compute the mean with the option `-metric-reduce.` That is, something as below (make the necessary changes to your file names and directory structure):

```
for subj in ${LIST_SUBJECTS} ; do
   wb_command -surface-vertex-areas ${subj}/L_midthick.surf.gii ${subj}/L_midthick_va.shape.gii
done

L_MERGELIST=""
for subj in ${LIST_SUBJECTS} ; do
   L_MERGELIST="${L_MERGELIST} -metric ${subj}/L_midthick_va.shape.gii"
end
wb_command -metric-merge L_midthick_va.func.gii ${L_MERGELIST}

wb_command -metric-reduce L_midthick_va.func.gii MEAN L_area.func.gii
```

Repeat all steps for the right hemisphere. For details of these options, syntax and examples, consult the [Workbench Command documentation](https://www.humanconnectome.org/software/workbench-command). Alternatively, if for any reason the average areas cannot be produced, consider using the sphere surface and the vertex-area average files available [here](https://github.com/Washington-University/Pipelines/tree/master/global/templates/standard_mesh_atlases/resample_fsaverage); these areas were produced from the S900 HCP release.

Once the input files, surfaces and their average areas have been produced, PALM can be run three separate times as:

```
palm -i data_sub.nii -d design.mat -t design.con -o results_sub -T -logp [...]
palm -i data_L.func.gii -d design.mat -t design.con -o results_L_cort -T -tfce2D -s L_midthickness.surf.gii L_area.func.gii -logp [...]
palm -i data_R.func.gii -d design.mat -t design.con -o results_R_cort -T -tfce2D -s R_midthickness.surf.gii R_area.func.gii -logp [...]
```

The results are merged with the option `-cifti-create-dense-from-template` of `wb_command`, creating a new CIFTI file for each matching set of p-value maps. For example:

```
wb_command -cifti-create-dense-from-template data.dscalar.nii results_merged_tfce_tstat_fwep_c1.dscalar.nii -volume-all results_sub_tfce_tstat_fwep_c1.nii -metric CORTEX_LEFT results_L_cort_tfce_tstat_fwep_c1.gii -metric CORTEX_RIGHT results_R_cort_tfce_tstat_fwep_c1.gii
```

These new CIFTI files can be opened with `wb_view`. The corrected significance threshold can be determined via Bonferroni, i.e., \(-log10(alpha/N) = -log10(0.05/3) = 1.7782\), or with Šidák, i.e., \(-log10(1-(1-alpha)1/N) = -log10(1-(1-0.05)1/3) = 1.7708\), the latter being therefore slightly more powerful.

As in the previous example, if the input data are subjects from the Human Connectome Project, it's also necessary to supply the exchangeability blocks file (`-eb EB.csv`). More details [here](exchangeability_blocks.md) (in the section about the Human Connectome Project).

## Example 11: Using a log and/or configuration file

The various options available can make it difficult to keep track of which were used for any call of the command. To help with this, PALM saves a small text file with all the options used, which works as a log-file. This file is named `{outputprefix}_palmconfig.txt`. For the Example 1 above, the file would be like this:

```
# Configuration file for PALM.
# 12-Mar-2014 11:29:47

-i dataset4D.nii
-d design.mat
-t design.con
-m mask.nii
-f design.fts
-T
-C 3.1
-n 5000
-corrcon
-o myresults
```

The same log-file can also be used as the input for a second run. In this case, PALM should be called with the file name as the sole input argument:

```
palm myresults_palmconfig.txt
```