FROM debian
MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make wget && \
    apt-get clean




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



RUN STAR_FUSION_URL="https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v0.7.0/STAR-Fusion_v0.7.0_FULL.tar.gz" && \
       cd $SRC && \
       wget $STAR_FUSION_URL && \
       tar xvf STAR-Fusion_v0.7.0_FULL.tar.gz

ENV STAR_FUSION_HOME $SRC/STAR-Fusion_v0.7.0_FULL


##############
## STAR


RUN STAR_URL="https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/2.5.1b.tar.gz -C $SRC && \
    mv $SRC/STAR-2.5.1b/bin/Linux_x86_64_static/STAR /usr/local/bin


###################
## FusionInspector

RUN FI_URL="https://github.com/FusionInspector/FusionInspector/releases/download/v0.5.0/FusionInspector_v0.5.0_FULL.tar.gz" && \
	wget -P $SRC $FI_URL && \
    tar -xvf $SRC/FusionInspector_v0.5.0_FULL.tar.gz -C $SRC

ENV FUSION_INSPECTOR_HOME $SRC/FusionInspector_v0.5.0_FULL

RUN apt-get install -y curl  ## move this to top

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN apt-get install -y libdb-dev

RUN cpanm install DB_File

RUN cpanm install Set::IntervalTree


##########
## Trinity

RUN apt-get install -y zlib1g-dev bzip2 libncurses5-dev

RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz" && \
    wget -P $SRC $TRINITY_URL && \
    tar -xvf $SRC/v2.1.1.tar.gz -C $SRC && \
    cd $SRC/trinityrnaseq-2.1.1 && make

ENV TRINITY_HOME $SRC/trinityrnaseq-2.1.1




RUN cp $TRINITY_HOME/trinity-plugins/htslib/bgzip $BIN

RUN cp $TRINITY_HOME/trinity-plugins/BIN/samtools $BIN

RUN cpanm install URI::Escape

###############################
## Install
RUN apt-get install -y openjdk-7-jre

COPY util/*.pl $BIN/

