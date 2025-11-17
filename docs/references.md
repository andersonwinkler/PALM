# References

The main reference for PALM is the same as for [randomise](https://fsl.fmrib.ox.ac.uk/fsl/docs/statistics/randomise.html):

> Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. [Permutation inference for the general linear model](https://doi.org/10.1016/j.neuroimage.2014.01.060). NeuroImage, 2014;92:381-397. (Open Access)

For correction across contrasts (option `-corrcon`), the most detailed reference is:

> Alberton BAV, Nichols TE, Gamba HR, Winkler AM. [Multiple testing correction over contrasts for brain imaging.](https://dx.doi.org/10.1016/j.neuroimage.2020.116760) Neuroimage. 2020 Mar 19:116760. (Open Access)

For Non-Parametric Combination (NPC; options `-npc`, `-npcmod` and `-npccon`), classical multivariate tests (MANOVA, MANCOVA; option `-mv`) assessed with permutations, and for correction over modalities and/or contrasts (options `-corrmod` and `-corrcon`), the reference is:

> Winkler AM, Webster MA, Brooks JC, Tracey I, Smith SM, Nichols TE. [Non-Parametric Combination and related permutation tests for neuroimaging.](https://dx.doi.org/10.1002/hbm.23115) Hum Brain Mapp. 2016 Apr;37(4):1486-511. (Open Access)

For the multi-level exchangeability blocks (options `-eb`,
`-vg`, and for HCP data), the reference is:

> Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM. [Multi-level block permutation.](https://doi.org/10.1016/j.neuroimage.2015.05.092) Neuroimage. 2015;123:253-68. (Open Access)

For the accelerated permutation inference (options `-accel` or `-approx`), the reference is:

> Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM. [Faster permutation inference in brain imaging. Neuroimage.](https://doi.org/10.1016/j.neuroimage.2016.05.068) 2016 Jun 7;141:502-516. (Open Access)

For the false discovery rate (FDR) approach used in PALM (option `-fdr`), the reference is:

> Winkler AM, Taylor PA, Nichols TE, Rorden C. [False Discovery Rate and Localizing Power](https://arxiv.org/abs/2401.03554). arXiv, 2024;2401.03554. (Preprint)

For additional theory of permutation tests in neuroimaging, please see and cite:

> Nichols TE, Holmes AP. [Nonparametric permutation tests for functional neuroimaging: a primer with examples.](http://dx.doi.org/10.1002/hbm.1058) Hum Brain Mapp. 2002 Jan;15(1):1-25.

> Holmes AP, Blair RC, Watson JD, Ford I. [Nonparametric analysis of statistic images from functional mapping experiments.](http://dx.doi.org/10.1097/00004647-199601000-00002) J Cereb Blood Flow Metab. 1996 Jan;16(1):7-22.