FROM trinityctat/ctatfusion:1.1.0

RUN mv /usr/local/src/STAR-Fusion*/* /usr/local/bin

COPY util/fusion_pipe_runner.py /usr/local/bin

RUN mkdir /workdir

WORKDIR /workdir

