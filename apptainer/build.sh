#!/bin/bash

set -euo pipefail

dst_dir=/opt/src/seed-finding
apptainer exec --bind "$(pwd)/..:${dst_dir}" seed-finding.sif bash -c "mkdir -p ${dst_dir} && cd ${dst_dir} && make -j $(nproc --all)"
