SUNDIALS SWIG-Fortran
====

The SUNDIALS SWIG-Fortran code is used to generate Fortran 2003 bindings
to the SUNDIALS C API in order to provide a scalable and sustainable
Fortran 2003 interface to SUNDIALS. The intent is to closely mimic the
C API while providing an idiomatic Fortran interface.

## Getting SWIG-Fortran

We use the SWIG-Fortran fork of SWIG created by Seth R. Johnson @ ORNL.
The upstream repository is maintained on [GitHub](https://github.com/swig-fortran/swig).
We maintain [a fork of SWIG-Fortran](https://github.com/sundials-codes/swig)
that is held at the last working commit and includes any of our own bug fixes.
So if the the latest swig obtained from the actual SWIG-Fortran repository
doesn't work and the fixes required to make it work are non-trivial, you can
clone our fork.

To build SWIG-Fortran (and optionally install it on your system), first complete
the following commands:

```bash
$ git clone https://github.com/sundials-codes/swig
$ cd swig
$ ./autogen.sh
$ ./configure --prefix=/my/install/location
```

*Note, for the complex support prototype, you need the backport-complex-support branch
of sundials-codes/swig.*

At this point you should check and make sure that autoconf will in fact build
the Fortran generator. The final line of the configure output should say
something like:

```bash
The SWIG test-suite and examples are configured for the following languages:
fortran
```

If it does not report back fortran. Try rerunning configure like so:

```bash
$ ./configure --with-fortran=/path/to/fortran/compiler --prefix=/my/install/location
```

Finally, proceed to make and optionally install SWIG.

```bash
$ make
$ make install # optional
```

## How to regenerate the interfaces

To regenerate the interfaces that have already been created, configure the
standalone CMake project in `sundials/swig` with the desired
`SUNDIALS_INDEX_SIZE` and `SUNDIALS_SCALAR_TYPE`, then build the default target:

```bash
$ cmake -S swig -B build-swig
$ cmake --build build-swig
```

The standalone project defaults to `SUNDIALS_INDEX_SIZE=64` and
`SUNDIALS_SCALAR_TYPE=REAL`. To regenerate a different variant, configure those
cache variables explicitly. For example, to generate the complex 32-bit
wrappers:

```bash
$ cmake -S swig -B build-swig-c32 \
    -DSUNDIALS_INDEX_SIZE=32 \
    -DSUNDIALS_SCALAR_TYPE=COMPLEX
$ cmake --build build-swig-c32
```

If SWIG-Fortran is not on `PATH`, set `SWIG_EXECUTABLE` when configuring:

```bash
$ cmake -S swig -B build-swig -DSWIG_EXECUTABLE=/path/to/swig
```

It is also possible to force the SWIG generation to run as part of the main
SUNDIALS build by enabling the optional top-level switch `SUNDIALS_ENABLE_SWIG`
Once enabled, the SWIG generation step will run during the CMake build phase
(before compilation steps). You can still invoke the SWIG generation step explicitly
with the `sundials_swig` target if needed:

```bash
$ cmake -S . -B build -DSUNDIALS_ENABLE_FORTRAN=ON -DSUNDIALS_ENABLE_SWIG=ON
$ cmake --build build
# or explicitly:
$ cmake --build build --target sundials_swig
```

In that mode `sundials_swig` follows the configured `SUNDIALS_INDEX_SIZE` and
`SUNDIALS_SCALAR_TYPE`. Both the standalone `swig/` project and the top-level
integration now follow the configured cache variables rather than exposing
separate per-variant build targets.

**This will replace all the generated files in `sundials/src`.**


## Creating a new interface

To create an interface to a new SUNDIALS module or package, the easiest thing
to do is copy one of the existing `.i` files for a module that is similar.
Then add the file to the `SUNDIALS_SWIG_INTERFACE_SPECS` list in
`swig/CMakeLists.txt`.

It may be useful to first read the "SUNDIALS Fortran 2003 interface" section
of the  user guide before trying to develop new interfaces.


## SWIG-Fortran documentation

The SWIG-Fortran documentation is in the SWIG repository: `swig/Doc/Manual`.
The Fortran specific section is in the file `swig/Doc/Manual/Fortran.html`.

## Other notes

The `_SUNDIALS_STRUCT_` macro (defined in `sundials_types.h`) must be used when
declaring a `struct` which will be interfaced to in Swig
(e.g. the `_generic_N_Vector` structure). The macro is defined as a `struct`
unless generating the SWIG interfaces - in that case it is defined as nothing.
This is needed to work around a bug in SWIG which prevents it from properly parsing
our generic module structures.
