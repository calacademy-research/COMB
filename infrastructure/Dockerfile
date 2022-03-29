# Most contributors can ignore this file.
#
# It serves as documentation of the system-wide dependencies necessary to get
# from a bare Debian Linux distribution to a working development environment.
#
# Docker / podman / etc. are container management systems that, with this file,
# can create a somewhat sandboxed environment with this project's development
# setup, but most contributors will probably prefer the host OS for development
# work.

FROM debian:bullseye

LABEL org.opencontainers.image.source https://github.com/calacademy-research/COMB

RUN apt -y update
RUN apt -y install gnupg sudo git \
    libgdal-dev libnlopt-dev libudunits2-dev libssl-dev

RUN apt-key adv --keyserver keyserver.ubuntu.com \
    --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN echo -n 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/' \
    >> /etc/apt/sources.list
RUN apt -y update
RUN apt -y install r-base jags
RUN apt -y upgrade

# lubridate will call container-unfriendly timedatectl unless this is set
RUN echo -n 'export TZ=UTC' >> /root/.bashrc

WORKDIR /root
RUN git clone https://github.com/calacademy-research/COMB.git

WORKDIR /root/COMB
RUN Rscript -e 'renv::restore()'
