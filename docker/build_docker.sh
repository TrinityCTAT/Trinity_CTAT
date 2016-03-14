#!/bin/bash

VERSION=$(cat VERSION)

docker build -t trinityctat/ctatfusion:$VERSION --rm=true .

