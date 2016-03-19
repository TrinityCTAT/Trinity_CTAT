#!/bin/bash

VERSION=$(cat VERSION)

sudo docker push trinityctat/ctatfusion:$VERSION
