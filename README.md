# seed_finding

Build instructions on linux system with conda, installed:
```shell
git clone git@github.com:saskeli/seed_finding.git
conda create -n cxx14 -c conda-forge cxx-compiler gsl zlib gxx=14 make
conda activate cxx14
cd seed_finding
make
```
