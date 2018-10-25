FROM neurodebian:stretch

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && apt-get install -y --no-install-recommends less \
                                                                     ca-certificates \
                                                                     wget \
                                                                     unzip \
                                                                     octave \
                                                                     octave-general \
                                                                     octave-signal \
                                                                     octave-image \
                                                                     liboctave-dev \
                       && apt-get clean \
                       && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# download HCP workbench (1.3.0) for CIFTI support
WORKDIR /opt
RUN wget -q https://ftp.humanconnectome.org/workbench/workbench-linux64-v1.3.0.zip -O wb.zip \
    && unzip wb.zip \
    && rm wb.zip

# copy repo into /opt/palm
WORKDIR palm
COPY . .

# recompile mex files to be compatible with our system
WORKDIR fileio
RUN cd \@gifti/private \
    && mkoctfile --mex zstream.c \ 
    && cd ../../\@xmltree/private \
    && mkoctfile --mex xml_findstr.c \
    && cd ../../\@file_array/private \
    && mkoctfile --mex file2mat.c \
    && mkoctfile --mex mat2file.c

# add palm and workbench to our PATH
ENV PATH=/opt/palm:/opt/workbench/bin_linux64:$PATH

# enter with call to palm executable
WORKDIR /working
ENTRYPOINT ["palm"]
