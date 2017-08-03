FROM ubuntu:16.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                                       unzip && \
    apt-get clean


