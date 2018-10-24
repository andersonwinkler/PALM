To learn about PALM, please visit **http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM**

On the link above there is no history about the package itself, so here it goes some. In one way or another, work on a permutation tool has been ongoing for a while, with some of the functions having been loosely written back in 2008, if not earlier. These functions began to be assembled and integrated by mid-2013. Between Oct/2013 and Mar/2014, Git was used locally for version control, but eventually PALM went on without it. Various early alpha versions circulated to collaborators, with the first public release in late Feb/2015 in the FSL website. These releases were, and continue to be, in the form of tarballs. Since all were kept, it was easy to retroactively apply the commits to the same local repository that had been neglected, and finally make it public on GitHub today, 04/Jul/2015.


----

### Docker usage

[Docker](https://www.docker.com/get-started) is an application that allows encapsulation of software environments. These environments are referred to as `containers` and `images`. [Click here](https://docs.docker.com/get-started/#containers-and-virtual-machines) for more information on Docker.

##### Build
PALM includes a `Dockerfile` with instructions to build the image.

1) Clone (`git clone`) or download this repository.
2) `cd` into this repo and run `docker build --rm -t mypalm .`

##### Usage
Once you have built (or pulled) the image, you can run it with local data. You can pass local directories into docker containers using the `-v` flag. Here's an example command for running PALM:

```bash
docker run --rm -it -v /path/to/my/data:/data mypalm -i /data/mydata.nii -d /data/design.mat -t /data/design.con
```
