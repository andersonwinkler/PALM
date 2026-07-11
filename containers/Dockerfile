# =============================================================================
# Dockerfile — Docker image for PALM
# =============================================================================
#
# Use this file if you want to build a container for PALM with Docker 
# (or Podman with Docker-compatible commands).
# For Apptainer / Singularity, use Apptainer.def in this directory instead
# (do not pass this Dockerfile to "apptainer build" or "singularity build").
#
# Contents of the image:
#   - GNU Octave (base: docker.io/gnuoctave/octave)
#   - PALM from the GitHub "dev" branch
#   - Octave packages: image, statistics, specfun
#   - Compiled lib/ extensions ("make octave", including hdf5oct)
#
# Build (run from this directory):
#   docker build -t palm .
#   podman build -t palm .
#
# Run (arguments after the image name are passed to palm):
#   docker run --rm palm -h
#   podman run --rm palm -h
#
# Mount a host data directory and run an analysis, e.g.:
#   docker run --rm -v /path/to/data:/data -w /data palm \
#     -i /data/mydata.nii -d /data/design.mat -t /data/design.con \
#     -o /data/results
#
# Override the entrypoint to open an Octave shell inside the image:
#   docker run --rm -it --entrypoint octave palm
#
# =============================================================================

FROM docker.io/gnuoctave/octave:7.3.0

RUN cd /opt \
 && git clone --branch dev "https://github.com/andersonwinkler/PALM" PALM \
 && sed -i \
      -e 's,^#\(OCTAVEBIN=/usr/bin/octave\),\1,' \
      -e 's,^\(OCTAVEBIN=/opt/local/bin/octave\),#\1,' \
      /opt/PALM/palm \
 && chmod +x /opt/PALM/palm \
 && echo -e "\n\naddpath('/opt/PALM')\n" \
      >> "$(octave-config --print LOCALSTARTUPFILEDIR)/octaverc" \
 && octave --quiet --eval "pkg install -forge image statistics specfun" \
 && make -C /opt/PALM/lib distclean \
 && make -C /opt/PALM/lib octave

WORKDIR /opt/PALM
ENTRYPOINT ["/opt/PALM/palm"]
