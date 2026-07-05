# File formats

Supported file formats are listed below. Formats are identified by the file extension. Different options may accept only certain types of files (for example, the option `-s` will only accept surface geometry formats).

| Extension | Read as |
| --- | --- |
| `.nii`, `.hdr`, `.img` | NIFTI file. Can be used to input data (option `-i`) and mask (`-m`). Note that the old ANALYZE format is not supported. |
| `.nii.gz` | Compressed NIFTI file. Should be uncompressed manually (with `gunzip`) and loaded as `nii` (it uses more disk space, but far less memory). Or can be read directly with the option `-noniiclass`, but beware that you need to have enough memory to load all the data. |
| `.mgh`, `.mgz` | FreeSurfer data, either volumetric or surface-based. Can be used to input data (option `-i`) or mask for the same kind of data (option `-m`). |
| `.dpv`, `.dpf`, `.dpx` | Data-per-vertex (vertexwise), data-per-face (facewise), or unspecified surface-based data. Can be used to specify a mask (option `-m`) for surface-based data. Multiple such files can be merged into a `csv` table, and then input with the option `-i`; see details below. |
| `.csv` | Table in `csv` format. Can be used to specify input data (option `-i`), mask for the same kind of data (option `-m`), design matrix (option `-d`), contrasts (options `-t` and `-f`), exchangeability blocks (`-eb`) and variance groups (`-vg`). |
|`.parquet`|Apache Parquet. Can be used to specify large tables. Make sure you have enough memory to load all data.|
|`.h5`, `.hdf5`|HDF5 (Hierarchical Data Format). Can be used to supply very large files. Make sure you have enough memory to load all data.|
| `.mat`, `.con`, `.fts`, `.grp` | Table in the FSL `vest` format. Can be used essentially in the same way as `csv` files, i.e., to specify input data (option `-i`), mask for the same kind of data (`-m`) design matrix (`-d`), contrasts (options `-t` and `-f`), exchangeability blocks (`-eb`) and variance groups (`-vg`). |
| `.mset` | Multiple tables (arrays) in a single ASCII file. This format can be used with the option `-con`. |
| `.srf` | Surface in ASCII format (option `-s`). |
| `.inflated`, `.nofix`, `.orig`, `.pial`, `.smoothwm`, `.sphere`, `.reg`, `.white` | Files with these extensions are read as FreeSurfer surface (option `-s`). |
| `.gii` | GIFTI file. |
| `.dtseries.nii`, `.ptseries.nii`, `.dscalar.nii`, `.pscalar.nii` | CIFTI file. Support is currently available for `.dscalar.nii`, `.dtseries.nii`, `.pscalar.nii` and `.ptseries.nii.` It is expected that in the future there will be complete support for CIFTI files. |

Files with extensions not listed above won't be read.

### Support for NIFTI files

Support for uncompressed NIFTI files (extension `.nii`) is provided, internally, by the publicly available NIFTI class. This allows reading and writing even huge files without using too much computer memory. However, the NIFTI class does not operate on compressed files, i.e., with extension `.nii.gz`. To read these files, it is recommended that they are uncompressed first (with `gunzip`).

Alternatively, if the datasets are small, the NIFTI class can be disabled with the option `-noniiclass`. This allows reading and writing `.nii.gz` files directly. However, if the files are too large, this can easily use all the computer memory and the system may become unstable/unusable. The option `-noniiclass` should be used with caution for large datasets. If the option `-noniiclass` is provided and PALM is running with MATLAB as the engine, then if the Image Processing Toolbox is installed, `.nii.gz` files will be read with the command `niftiread`; otherwise, i.e., if the PALM is running with Octave as the engine or if the Image Processing Toolbox is not available, then if the option `-noniiclass` is provided, `.nii.gz` files will be read using the command `load_nifti`, which is available internally within PALM (courtesy from the FreeSurfer developers).

The NIFTI class is used by default. It is provided with precompiled binaries for MATLAB for various platforms, and for Octave for most 64-bit Linux distributions. Nonetheless, if compilation is needed, use:

```
cd /full/path/to/palm/fileio/@file_array/private
./compile.sh
```

### Support for FreeSurfer files

FreeSurfer binary files (surfaces and curvatures) are read directly. Surface files in ASCII format are also read directly by PALM as long as their file extension is srf.

FreeSurfer "curvature" files converted to pseudo-volumes (with extension `.mgh` or `.mgz`) are read directly. For ASCII (`asc`/`dpv`/`dpf`/`dpx`), these need to be merged across subects and converted to `.csv` tables. To generate a valid `.csv` file, use the command `dpx2csv` (available [here](https://raw.githubusercontent.com/andersonwinkler/toolbox/master/bin/dpx2csv)); this will initially create a file with one column per subject and one row per vertex or face, then transpose the rows and columns of this file with the command `transpose` (available [here](https://raw.githubusercontent.com/andersonwinkler/toolbox/master/bin/transpose)), or use the option `-transposedata`.

### Support for MZ3 files

[MZ3](https://github.com/neurolabusc/surf-ice/tree/master/mz3) is a highly efficient, fast and compressed format to store triangular surfaces developed by Chris Rorden and colleagues, and is the default format used by SurfIce. Read/write support is provided natively.

### Support for HDF5 files

[HDF5 (Hierarchical Data Format)](https://www.hdfgroup.org/) is a high-performance, open-source file format designed for storing and managing complex, multi-dimensional scientific datasets. PALM can read and write HDF5 files natively with MATLAB; for Octave, the package `hdf5oct` (details [here](https://gnu-octave.github.io/packages/hdf5oct/)) must be installed. To install it, run from the Octave prompt:

```
pkg install -forge hdf5oct
```

The user needs to indicate which specific datablock is to be used from the file (a single HDF5 file can hold multiple multidimensional arrays). A valid specificication for a datablock is as:

```
/path/to/file.h5:/path/to/data:N
```

where `/path/to/file.h5` is the path to the HDF5 file (can use absolute or relative paths; if no path is provided, the file is assumed to exist in the current directory); `/path/to/data` is the full path to the multidimensional array that is intended to be used, and `N` is an integer that indicates dimension along which the data should be permuted.

If an HDF5 file is specified as output (e.g., `-o /my/directory/myresults.h5`), then all outputs will be stored into the same HDF5 file, which is a convenient way to store an entire analysis into a single file.

HDF5 files can be loaded with any extension (including `.mat` and `.nwb`), but PALM must be able to recognize these extensions. Make sure to list them in `palm_defaults.m`, under the variable `opts.hdf5`.

### Support for Parquet files

[Apache Parquet](https://parquet.apache.org/) files (`.parquet`) are supported natively with MATLAB. With Octave, support requires that [DuckDB](https://duckdb.org/) is installed in the system, and reading/writing uses a `.csv` file as intermediate.

### Support for CSV files

Although PALM was created with imaging in mind, it can be used for non-imaging data. Any kind of data that can be arranged in a table and saved as comma-separated values file (`.csv`) (i.e., any data) can be used as input. Each row constitutes a subject or observation, and each column represents a measurement. The `.csv` files must contain numeric fields only, i.e., title labels are not accepted and will cause errors.

If the data is arranged as space- or tab-separated values instead of comma-separated values, it can be converted quickly using `awk` (or `gawk`):

```
awk 'BEGIN { OFS="," } { $1=$1; print $0 }' oldfile.txt > newfile.csv
```

If your data has any other separator, e.g., a semicolon, the same applies with using a small modification:

```
awk 'BEGIN { FS="x"; OFS="," } { $1=$1; print $0 }' oldfile.txt > newfile.csv
```

where the `"x"` should be replaced by the separator in the original table (e.g., a semicolon).

### Support for MSET files

To allow multivariate contrasts, PALM uses a simple ASCII format that can contain multiple matrices. An example of such file is:

```
Matrix 1 3
1 -1 0

Matrix 2 3
1 -1 0
1 0 -1
```

In the example, the first matrix is defined as having 1 row and 3 columns (the numbers after the keyword `Matrix`). The second matrix is defined as having 2 rows and 3 columns. Files as these are meant to be used to test multivariate hypotheses as `H: C'*Beta*D` using the option `-con <file1> <file2>`, where `<file1>` is a file with multiple contrasts C and `<file2>` with multiple contrasts D. Each contrast in C pairs with a contrast in D.

### Support for CIFTI files

Both the surface and volume components in CIFTI files that are of the types `dscalar`, `pscalar`, `pconnscalar`, `dtseries`, `ptseries`, `pconnseries`, `dconn`, `pconn`, `pdconn`, `dpconn`, `dfan`, `dfibersamp`, `dfansamp`,  `dlabel`, and `merge` can be read directly. If no spatial statistics are requested (i.e., no TFCE or cluster-level inference), all “grayordinates” (whether surface or volume) are processed. However, when spatial statistics are requested, it is currently necessary to first manually split the CIFTI files into separate surface (GIFTI) and volume (NIFTI) components, which then can be loaded and processed (see examples [in this page](examples.md)). The [Connectome Workbench](http://www.humanconnectome.org/software/connectome-workbench.html) must be installed.

### Support for GIFTI files

GIFTI files are supported. However, note that the access to the actual data depends on parsing a potentially large XML tree, which can be slow in both MATLAB and Octave. Other equivalent formats are usually faster to load.

### Support for other file formats

If there is a file format that cannot be easily converted to `csv`, and which you'd like to be able to read directly in PALM, feel free to contact us. New formats may be added for future releases.