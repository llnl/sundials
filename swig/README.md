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

To regenerate the interfaces that have already been created. Simply run
`make all32 all64` from the `sundials/swig` directory.
**This will replace all the generated files in `sundials/src`.**


## Creating a new interface

To create an interface to a new SUNDIALS module or package, the easiest thing
to do is copy one of the existing `.i` files for a module that is similar.
Then add the file to the appropriate section of the Makefile.

It may be useful to first read the "SUNDIALS Fortran 2003 interface" section
of the  user guide before trying to develop new interfaces.


## SWIG-Fortran documentation

The SWIG-Fortran documentation is in the SWIG repository: `swig/Doc/Manual`.
The Fortran specific section is in the file `swig/Doc/Manual/Fortran.html`.

## Other notes

The `SUNDIALS_STRUCT` macro (defined in `sundials_types.h`) must be used when
declaring a `struct` which will be interfaced to in Swig
(e.g. the `N_Vector_` structure). The macro is defined as a `struct`
unless generating the SWIG interfaces - in that case it is defined as nothing.
This is needed to work around a bug in SWIG which prevents it from properly parsing
our generic module structures.
