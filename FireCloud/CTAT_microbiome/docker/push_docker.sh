#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/firecloud_ctatmicrobiome:${VERSION}
docker push trinityctat/firecloud_ctatmicrobiome:latest

