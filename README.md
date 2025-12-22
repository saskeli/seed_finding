# seed_finding

## Building

To clone the repository with submodules, please use `git clone --recursive https://github.com/saskeli/seed_finding.git`.

### With [Conda-build](https://docs.conda.io/projects/conda-build/)

A conda package can be built with Conda-build as follows. The build script has been tested with Conda-build 25.11.1.

1. `mamba create -n conda-build -c conda-forge conda-build conda-verify` (Conda-verify is not strictly required.)
2. `mamba activate conda-build`
3. `cd conda`
4. `./conda-build.sh`

Conda-build will then report the location of the package. An environment can be created with e.g. the following command:
```
mamba create -n seed-finding -c local -c conda-forge seed_finder
```

### By Creating a Build Environment with [Apptainer](https://apptainer.org/) (Singularity)

Apptainer can be used to create a build environment as follows.

1. `cd apptainer`
2. `apptainer build seed-finding.sif seed-finding.def`
3. `./build.sh`

The resulting binaries may not necessarily execute on your system. To run `seed_finder` from the build environment, use `./run-seed-finder.sh`. Any command line parameters given to the script will be passed to `seed_finder`.

### By Running Make Directly

The following software and libraries are required to build PanVC 3. The tested versions are also listed.

- C and C++ compilers. The C++ compiler needs to support C++23. ([GCC 14.2.0](https://gcc.gnu.org/) on Linux, [LLVM 21](https://llvm.org/) on macOS.)
- [autoconf 2.71](https://www.gnu.org/software/autoconf/)
- [aclocal (part of automake 1.16.5)](https://www.gnu.org/software/automake/)
- [Boost 1.83.0](https://www.boost.org)
- [GNU Make 4.3](https://www.gnu.org/software/make/)
- [GNU Scientific Library 2.8](https://www.gnu.org/software/gsl/)
- [Ragel 6.10](http://www.colm.net/open-source/ragel/)
- [zlib 1.2.13](https://zlib.net)

On Ubuntu, the following command should install all the required packages.

```
apt-get install -y --no-install-recommends autoconf automake build-essential gcc-14 g++-14 git libboost-all-dev libgsl-dev make pkg-config ragel zlib1g-dev
```

After installing the prerequisites, the software can be built with e.g. `make -j16`.

## Running

Signal and background reads can be given as FASTA or FASTQ. Both (b)gzip-compressed and uncompressed inputs are accepted. Given *signal.fastq.gz* and *background.fastq.gz*, seed *k*-mers may be searched with a command similar to the following.

```
seed_finder -p 0.0001 --lf 2 -s --max-s 0 --pruning -b background.fastq.gz signal.fastq.gz
```

Please see `seed_finder --help` for all available command line parameters.
