#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build -t trinityctat/firecloud_ctatmicrobiome:${VERSION} .
docker build -t trinityctat/firecloud_ctatmicrobiome:latest .

