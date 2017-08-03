FROM trinityctat/ctatfusion:1.1.0

RUN mv /usr/local/src/STAR-Fusion*/* /usr/local/bin

COPY util/fusion_pipe_runner.py /usr/local/bin

## install picard for reverting bam to fastq

WORKDIR /usr/local/bin

RUN wget "https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar"


RUN mkdir /workdir

WORKDIR /workdir

