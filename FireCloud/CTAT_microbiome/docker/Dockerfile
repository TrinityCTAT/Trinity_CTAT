FROM ubuntu:16.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
                                       texlive-latex-base \
                                       default-jre \
                                       python-pip python-dev \
                                       gfortran \
                                       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
                                       libblas-dev gfortran git unzip ftp libzmq3-dev nano ftp fort77 libreadline-dev \
                                       libcurl4-openssl-dev libx11-dev libxt-dev \
                                       x11-common libcairo2-dev libpng12-dev libreadline6-dev libjpeg8-dev pkg-config \
                   && apt-get clean


WORKDIR /usr/local/bin

RUN wget "https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar"


#############################

## samtools

WORKDIR /usr/local/src

RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
    tar xvf samtools-1.5.tar.bz2  && \
    cd samtools-1.5 && \
    ./configure && make && make install


############################

## Centrifuge


RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip && \
     unzip centrifuge-1.0.3-beta-Linux_x86_64.zip && \
     ln -s /usr/local/src/centrifuge-1.0.3-beta/centrifuge* /usr/local/bin/.


COPY centrifuge/centrifuge-kreport /usr/local/bin

###################

COPY util/* /usr/local/bin/


RUN mkdir /workdir

WORKDIR /workdir

