#!/bin/bash

VERSION=$(cat VERSION)

docker push trinityctat/ctatfusion:$VERSION
