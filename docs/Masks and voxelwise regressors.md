It is possible to include one regressor (explanatory variable, EV) for each variable from the input files, that is, voxelwise EVs (or vertexwise, facewise, edgewise, columnwise, depending on the type of data). This is done with the option `-evperdat`. Below, the term "voxelwise EV" is used generically, even if it may actually be a vertexwise or any other type of columnwise EV.

With or without `-evperdat`, each input file (supplied with -i) can have their own mask. If `-mv` or `-npc` is used, an intersection mask is produced and used for all the tests. For the voxelwise EVs, further masks are generated internally to remove tests that would have constant EVs, as well as to remove the invalid values NaN and Inf if present. These masks created internally can be saved with the option `-savemask`, which will create one file for each effective mask used for each modality.

Whatever is the case, the option `-m` pairs with the option `-i`, in the same order. The option `-m` does not pair with the option `-evperdat`. Masks for the voxelwise EVs are defined based on the image modalities for which they will be used. Thus, the number of times `-m` can be given is zero, one, or the same number of times as `-i` is used. Any other number of masks entered with `-m` will cause an error.

## Without voxelwise EVs

If there are no voxelwise EVs, the possible mask cases are:

### Single input file

If there is just one input file (`-i`) and:

* **No masks are supplied:** A mask is created for the input file.
* **A mask is supplied:** This mask is used for the input file.

### Multiple input files

If there are multiple input files and:

* **No masks are supplied:** A mask is created for each input modality.
* **Only one mask is supplied:** This mask is used for all input modalities, which therefore need to be all with compatible sizes/geometries.
* **One mask is supplied for each input:** It is used for each respective input, and do not need to be all of the same geometry (except for MV and NPC, for which they need to be compatible).

## With voxelwise EVs

Voxelwise EVs are entered with the option `-evperdat`. This option can take up to three arguments:

```
palm [...] -evperdat <file> [evpos] [desnum] [...]
```

Where `<file>` is the file with one voxelwise EV, `[evpos]` is the column number (position) of this EV in the design matrix, and `[desnum]` is the design number, i.e., the design in which this EV should be added. Designs (`-d`) are numbered in the order as entered. If `[desnum]` is omitted, the default is 1. If `[evpos]` is also omitted, the default is 1. If `[evpos]` is 1 and a `[desnum]` is provided, but the respective design doesn't exist, a design containing a single EV is created.

If there are voxelwise EVs, the possible mask cases depend on the number of input files, and whether there is one design per input (option `-designperinput`) or not.

### One design per input

If there is one input (`-i`) and one design (`-d`), or if there are more than one input files and the option `-designperinput` is used, and:

* **No masks are supplied:** For each input modality, a mask is created. For each EV file for that input, a mask is also created. The intersection of the input mask and the EV masks is used to mask the input and the EVs.
* **Only one mask is supplied:** This mask is used for all input modalities and all EV files, which therefore need to be all of the same size.
* **One mask is supplied for each input:** It is used for each respective input and for the respective EV files. The size of the masks need to match the respective inputs and EV files, but don't otherwise have to match each other (except for MV and NPC).

### All designs vs. all inputs

If there is more than input or design, the option `-designperinput` is omitted, and:

* **No masks are supplied:** For each input modality, a mask is created. For each EV file, a mask is also created. The intersection of all these masks, for all the inputs and EV files, is used for all of them. This means all input files and voxelwise EVs need to be of same size.
* **Only one mask is supplied:** This mask is used for all input modalities and all EV files, which therefore need to be all of the same size.

Other combinations of masks, inputs, designs, and EVs will cause an error and PALM will stop.

Note that voxelwise EVs cause the runs to take longer to complete.

## References

The main reference for PALM is the same as for randomise:

> Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. [Permutation inference for the general linear model.](http://www.sciencedirect.com/science/article/pii/S1053811914000913) NeuroImage, 2014;92:381-397 (Open Access)