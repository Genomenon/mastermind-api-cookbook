#!/bin/bash

# (sets $HERE to this dir)
pushd `dirname ${BASH_SOURCE[0]}` > /dev/null; HERE=`pwd`; popd > /dev/null

# mastermind-api-cookbook dir to mount in the container:
cookbook="$HERE/.."

# get the script name
pyscript=$1
shift

# Use python2 image from docker hub
docker_img="qzkc/python2.7:v2"

# Run python in the container, with whatever args.
# !! NB !! ANY FILENAMES PASSED AS ARGS NEED TO BE PREFIXED WITH /data/ SO THE SCRIPT CAN FIND THEM IN THE CONTAINER
docker run -v "${cookbook}":/data --rm ${docker_img} python data/${pyscript} $@

