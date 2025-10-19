# PALM: Permutation Analysis of Linear Models

* [[User guide]]
* [[Exchangeability blocks]]
* [[Joint inference]]
* [[Masks and voxelwise regressors]]
* [[Faster inference]]
* [[Viewing results]]
* [[Examples]]
* [[Frequently asked questions]]

**PALM — Permutation Analysis of Linear Models** — is a tool that allows inference using permutation methods, offering a number of features not available in other analysis software. These features currently include:

* Ability to work with volumetric and surface-based formats, including facewise data, as well as with non-imaging data;
* A range of various regression and permutation strategies;
* Statistics that are robust to heteroscedasticity;
* Shuffling of sets of observations, to allow, for instance, the analysis of certain designs with repeated measurements, with no missing data;
* Shuffling of observations with complex, tree-like covariance structure (such as for the Human Connectome Project);
* Permutation with sign-flipping (wild bootstrap);
* Modified Non-Parametric Combination (NPC) for joint inference over multiple modalities, or multiple contrasts, or both together, with various combining functions available;
* Classical multivariate statistics (MANOVA, MANCOVA) for joint inference over multiple modalities, assessed through robust permutation methods, and also parametrically when such approximations exist;
* Correction over multiple contrasts, multiple modalities, for images with or without the same size or geometry, including non-imaging data, controlling the FWER or FDR;
* Various acceleration methods based on the test statistics and their distributions.

PALM requires Matlab or Octave. It can be executed from inside either environment, or directly from the shell. It can also be called from scripts.

**PALM is experimental software.** As novel features are introduced, tested, verified, and validated, eventually they will be implemented and made available in randomise or in other tools. PALM is for users who are familiar with statistics and willing to use experimental analysis tools. Bugs (real or suspected) can be reported via [GitHub](https://github.com/andersonwinkler/PALM) or through the [FSL mailing list](https://fsl.fmrib.ox.ac.uk/fsl/docs/support.html).

To download PALM and for installation instructions, please visit the [[User guide]].

## References

The main reference for PALM is the same as for [randomise](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/randomise.html):

> Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. [Permutation inference for the general linear model](https://doi.org/10.1016/j.neuroimage.2014.01.060). NeuroImage, 2014;92:381-397 (Open Access)

For correction across contrasts (option `-corrcon`), the most detailed reference is:

> Alberton BAV, Nichols TE, Gamba HR, Winkler AM. [Multiple testing correction over contrasts for brain imaging.](https://dx.doi.org/10.1016/j.neuroimage.2020.116760) Neuroimage. 2020 Mar 19:116760. (Open Access)

For Non-Parametric Combination (NPC; options `-npc`, `-npcmod` and `-npccon`), classical multivariate tests (MANOVA, MANCOVA; option -mv) assessed with permutations, and for correction over contrasts and/or modalities (options -corrcon and -corrmod), the reference is:

> Winkler AM, Webster MA, Brooks JC, Tracey I, Smith SM, Nichols TE. [Non-Parametric Combination and related permutation tests for neuroimaging.](https://dx.doi.org/10.1002/hbm.23115) Hum Brain Mapp. 2016 Apr;37(4):1486-511. (Open Access)

For the multi-level exchangeability blocks (options `-eb`,
`-vg`, and for HCP data), the reference is:

> Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM. [Multi-level block permutation.](https://doi.org/10.1016/j.neuroimage.2015.05.092) Neuroimage. 2015;123:253-68. (Open Access)

For the accelerated permutation inference (options `-accel` or `-approx`), the reference is:

> Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM. [Faster permutation inference in brain imaging. Neuroimage.](https://doi.org/10.1016/j.neuroimage.2016.05.068) 2016 Jun 7;141:502-516. (Open Access)

For additional theory of permutation tests in neuroimaging, please see and cite:

> Nichols TE, Holmes AP. [Nonparametric permutation tests for functional neuroimaging: a primer with examples.](http://dx.doi.org/10.1002/hbm.1058) Hum Brain Mapp. 2002 Jan;15(1):1-25.

> Holmes AP, Blair RC, Watson JD, Ford I. [Nonparametric analysis of statistic images from functional mapping experiments.](http://dx.doi.org/10.1097/00004647-199601000-00002) J Cereb Blood Flow Metab. 1996 Jan;16(1):7-22.