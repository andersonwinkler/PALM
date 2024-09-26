FROM bids/base_fsl

ENV DEBIAN_FRONTEND noninteractive

# install octave 3.6.2 (less needed as frontend)
RUN sudo apt-get update && sudo apt-get install -y less
RUN sudo apt-get install -y software-properties-common python-software-properties
RUN sudo add-apt-repository ppa:octave/stable
RUN sudo apt-get install -y octave octave-general octave-signal imagemagick
# cleanup package manager
RUN sudo apt-get autoclean && sudo apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


ENV FSLDIR /usr/share/fsl/5.0
ENV FSLOUTPUTTYPE NIFTI_GZ
ENV FSLMULTIFILEQUIT TRUE
ENV FSLTCLSH /usr/bin/tclsh
ENV FSLWISH /usr/bin/wish
ENV FSLBROWSER /etc/alternatives/x-www-browser
ENV LD_LIBRARY_PATH /usr/lib/fsl/5.0:${LD_LIBRARY_PATH}
ENV PATH ${FSLDIR}/bin:/home/ubuntu/bin:${PATH}

RUN mkdir /palm
COPY palm* /palm/
RUN mkdir /palm/fileio/
COPY fileio /palm/fileio/

RUN wget https://launchpad.net/ubuntu/+archive/primary/+files/octave-image_2.2.0-3_amd64.deb
RUN sudo dpkg -i octave-image_2.2.0-3_amd64.deb
WORKDIR /palm
CMD [".palm"]
