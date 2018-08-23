FROM trinityctat/ctatfusion:1.5.0preA

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

RUN mv /usr/local/src/STAR-Fusion*/* /usr/local/bin


###################################
## install picard for reverting bam to fastq

WORKDIR /usr/local/bin

RUN wget "https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar"

#############################

COPY util/* /usr/local/bin/


RUN mkdir /workdir

WORKDIR /workdir

