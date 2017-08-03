#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/firecloud_ctatfusion:${VERSION}
docker push trinityctat/firecloud_ctatfusion:latest

