#!/bin/bash

VERSION=$(cat VERSION)

sudo docker build -t trinityctat/ctatfusion:$VERSION --rm=true .

