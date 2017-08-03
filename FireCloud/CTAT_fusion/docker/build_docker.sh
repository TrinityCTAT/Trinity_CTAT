#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build -t trinityctat/firecloud_ctatfusion:${VERSION} .
docker build -t trinityctat/firecloud_ctatfusion:latest .

