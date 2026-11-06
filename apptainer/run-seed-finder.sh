#!/bin/bash

set -euo pipefail

dst_dir=/opt/src/seed-finding
apptainer exec --bind "$(pwd)/..:${dst_dir}" seed_finding "${dst_dir}/seed_finder" "${@}"
