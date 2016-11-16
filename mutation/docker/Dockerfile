FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
ENV BCFTOOLS_VERSION 1.3.1
ENV SAMTOOLS_VERSION 1.3.1
ENV STAR_VERSION 2.5.2a
ENV SNPEFF_VERSION v4_3b_core
ENV PICARD_VERSION 2.0.1
### locations
ENV BIN /usr/local/bin
ENV SRC /usr/local/src
ENV MUTATION_HOME ${SRC}/Trinity_CTAT/mutation
ENV STAR_HOME ${BIN}/STAR
### URLS
ENV BCFTOOLS_URL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
ENV PICARD_URL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip
ENV SAMTOOLS_URL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
ENV SNPEFF_URL https://sourceforge.net/projects/snpeff/files/snpEff_${SNPEFF_VERSION}.zip/download
ENV STAR_URL https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz

## Set path
ENV PATH ${BIN}:${MUTATION_HOME}:${MUTATION_HOME}/src:${STAR_HOME}:${SRC}:${PATH}

######################
## Dependencies
######################
##############
## Helper tools
RUN apt-get update
RUN apt-get install -y git unzip wget software-properties-common

##############
## GenomeAnalysisTK-3.1-1-g07a4bf8
# This is pulled in during the build
# If you are using this container / image please make sure to register a license for GATK at
# https://software.broadinstitute.org/gatk/download/licensing.php
COPY GenomeAnalysisTK.jar /usr/local/bin

##############
## Java 1.7 for GATK
## /usr/lib/jvm/java-8-oracle/jre/bin/java
## JAVA 1.8 for PICARD tools
## /usr/lib/jvm/java-8-oracle/jre/bin/java
RUN echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee /etc/apt/sources.list.d/webupd8team-java.list && \
    echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list && \
    echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 && \
    apt-get update && \
    apt-get install -y oracle-java7-installer && \
    apt-get install -y oracle-java8-installer

##############
## PICARD tools 1.764
WORKDIR ${SRC}
RUN wget -P ${SRC} ${PICARD_URL} && \
    unzip ${SRC}/picard-tools-${PICARD_VERSION}.zip && \
    mv ${SRC}/picard-tools-${PICARD_VERSION}/* ${BIN} && \
    rm -r ${SRC}/picard-tools-${PICARD_VERSION} && \
    rm -r ${SRC}/picard-tools-${PICARD_VERSION}.zip

## Pulling in some old jars until the code is updated.
COPY AddOrReplaceReadGroups.jar /usr/local/bin
COPY MarkDuplicates.jar /usr/local/bin
COPY SortSam.jar /usr/local/bin

##############
## Python 2.7.9
RUN apt-get install -y python=2.7.11-1
## Requests
WORKDIR ${BIN}
RUN wget https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py && \
    pip install requests

##############
## R version 3.1.1
RUN apt-get install -y r-base=3.2.3-4

##############
## Samtools
WORKDIR ${BIN}
RUN wget -P ${BIN} ${SAMTOOLS_URL} && \
    tar -jxvf ${BIN}/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ${BIN}
WORKDIR ${BIN}/samtools-${SAMTOOLS_VERSION}
RUN make prefix=/usr/local install
WORKDIR ${BIN}/samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}
RUN make prefix=/usr/local install

##############
## STAR
WORKDIR ${SRC}
RUN wget -P ${SRC} ${STAR_URL} && \
    tar -xvf ${SRC}/${STAR_VERSION}.tar.gz -C ${SRC} && \
    mv ${SRC}/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR ${BIN} && \
    rm -r ${SRC}/STAR-${STAR_VERSION} && \
    rm ${SRC}/${STAR_VERSION}.tar.gz

##############
## BCFTOOLS
WORKDIR ${BIN}
RUN wget -P ${BIN} ${BCFTOOLS_URL} && \
    tar -xvf ${BIN}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C ${BIN}
WORKDIR ${BIN}/bcftools-${BCFTOOLS_VERSION}
RUN make && make install

#############
## SNPEFF
WORKDIR ${BIN}
RUN wget -P ${BIN} ${SNPEFF_URL} && \
    unzip ${BIN}/download && \
    rm download
WORKDIR ${BIN}/snpEff
RUN java -jar snpEff.jar download hg19


######################
## Tool installation
######################

###############
## CTAT Mutation
WORKDIR ${SRC}
RUN git clone --recursive https://github.com/NCIP/Trinity_CTAT.git

######################
# Specify default behavior
# -h ?
######################
CMD rnaseq_mutation_pipeline.py --help
