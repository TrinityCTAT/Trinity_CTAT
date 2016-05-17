FROM debian
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
				       texlive-latex-base \
                                       openjdk-7-jre \
				       python-pip python-dev \
    && apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install DB_File

RUN cpanm install Set::IntervalTree

RUN cpanm install URI::Escape

RUN pip install pysam


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

ENV DATA /usr/local/data
RUN mkdir $DATA


######################
## Tool installations:
######################

###############
## STAR-Fusion:

RUN cd $SRC && \
     git clone --recursive https://github.com/STAR-Fusion/STAR-Fusion.git

ENV STAR_FUSION_HOME $SRC/STAR-Fusion

##############
## STAR

RUN RELEASE="2.5.2a" && STAR_URL="https://github.com/alexdobin/STAR/archive/${RELEASE}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${RELEASE}.tar.gz -C $SRC && \
    mv $SRC/STAR-${RELEASE}/bin/Linux_x86_64_static/STAR /usr/local/bin


###################
## FusionInspector

RUN cd $SRC && \
    git clone --recursive https://github.com/FusionInspector/FusionInspector.git

ENV FUSION_INSPECTOR_HOME $SRC/FusionInspector


##########
## Trinity


RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz" && \
    wget -P $SRC $TRINITY_URL && \
    tar -xvf $SRC/v2.1.1.tar.gz -C $SRC && \
    cd $SRC/trinityrnaseq-2.1.1 && make

ENV TRINITY_HOME $SRC/trinityrnaseq-2.1.1


RUN cp $TRINITY_HOME/trinity-plugins/htslib/bgzip $BIN

RUN cp $TRINITY_HOME/trinity-plugins/BIN/samtools $BIN

RUN cp $TRINITY_HOME/trinity-plugins/htslib/tabix $BIN


#############
## Oases

RUN VELVET_URL="http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz" && \
    wget -P $SRC $VELVET_URL && \
    tar xvf $SRC/velvet_1.2.10.tgz -C $SRC && \
    ln -s $SRC/velvet_1.2.10 $SRC/velvet && \
    cd $SRC/velvet && \
    make && \
    cp velveth velvetg $BIN/


RUN OASES_URL="https://www.ebi.ac.uk/~zerbino/oases/oases_0.2.08.tgz" && \
    wget -P $SRC $OASES_URL && \
    tar -xvf $SRC/oases_0.2.08.tgz -C $SRC && \
    cd $SRC/oases_0.2.8 && \
    make && \
    cp oases $BIN/


##############
## DISCASM

RUN cd $SRC && \
    git clone --recursive https://github.com/DISCASM/DISCASM.git

ENV DISCASM_HOME $SRC/DISCASM


###############################
## Install

COPY PerlLib $SRC/

ENV PERL5LIB $SRC:${PERL5LIB}



COPY util/*.pl $BIN/

