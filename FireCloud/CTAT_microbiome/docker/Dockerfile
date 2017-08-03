FROM ubuntu:16.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                                       unzip && \
    apt-get clean


WORKDIR /usr/local/bin

RUN wget "https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar"

WORKDIR /usr/local/src

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip && \
     unzip centrifuge-1.0.3-beta-Linux_x86_64.zip && \
     ln -s /usr/local/src/centrifuge-1.0.3-beta/centrifuge /usr/local/bin/.


RUN mkdir /workdir

WORKDIR /workdir

