#!/bin/bash

set -euxo pipefail

if [ -f local.mk ]
then
	echo "ERROR: local.mk already exists."
	exit 1
fi

echo "Generating local.mk"
m4 -D CONDA_BUILD_PREFIX="${BUILD_PREFIX}" conda/local.mk.m4 > local.mk

echo "Running make"
make -j "${CPU_COUNT}"

echo "Copying build products"
dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp seed_finder "${dst_bin}/"
