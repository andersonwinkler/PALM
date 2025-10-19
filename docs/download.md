
# Downloading PALM

To obtain the latest version, please visit [the repository on GitHub](https://github.com/andersonwinkler/PALM).

PALM can run as a standalone command (i.e., executed directly from the command line/terminal), or inside Octave or Matlab. In Linux and in Mac systems, it can be executed in any of these ways. In Windows, it can be executed inside Matlab or Octave.

### Running as a standalone command

It may be much simpler to run PALM as a command directly from the shell in Linux or in Mac, and it can easily be called from scripts. To do so:

1. Uncompress the downloaded file.
2. Open the file `palm` (not to be confused with `palm.m`). This is a script, inside which you can set whether the script should use Octave or Matlab. If the one that is chosen isn't in the `$PATH` variable, make sure to specify also the path to the directory that contains the executable for either of these (whichever you choose).
3. To invoke PALM, simply type `./palm` in the directory where it was installed. The path can also be added to the system's `$PATH` variable, so that it can be easily called from any directory just by typing `palm.`

If you are using Mac, and choose Octave, note that reading of NIFTI files need the option `-noniiclass` to work properly (more details below).

### Running inside Matlab

Uncompress the downloaded file, open Matlab, and add the newly created directory to the Matlab path (menu *File --> Set Path*). Typing `palm` at the prompt without arguments shows usage information. This works for Linux, Mac and Windows.

### Running inside Octave

Uncompress the downloaded file, start Octave, and add the newly created directory to the Octave path. This can be done with the command `addpath`:

```
addpath('/full/path/to/palm')
```

This line can be added to the `~/.octaverc` file, so that the change becomes permanent (if the `~/.octaverc` doesn't exist, an empty file can be created, then the line added). Typing `palm` at the prompt without arguments shows usage information. This works for Linux and Mac. It can work also with Windows if the paths are entered using the Windows filesystem convention (e.g., `C:\example`).

For Octave for Mac, note that reading of NIFTI files need the option `-noniiclass` to work properly (more details below).

For spatial statistics (cluster extent, cluster mass, and TFCE), the Octave package `image` is required.
