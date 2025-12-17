# Clang 20.0.0 reports syntax errors in libc++ that have
# to do with availability macros. (This could be a Conda-specific
# issue.) We try to mitigate by disabling the checks.
CPPFLAGS += -D_LIBCPP_DISABLE_AVAILABILITY -isystem CONDA_BUILD_PREFIX/include
LDFLAGS += -L CONDA_BUILD_PREFIX/lib
