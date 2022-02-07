# CapacitatedFacilityLocationBlock

This project centers on CapacitatedFacilityLocationBlock, an
implementation of the Block concept for a "pretty basic version" of
the Capacitated Facility Location (CFL) problem, a.k.a. the
Capacitated Warehouse Location (CWL) problem.

This class only represent the "basic version" of CFL and it is
primarily intended as a "didactic" implementation for showing some
of the features of SMS++, among which:

- CapacitatedFacilityLocationBlock supports a number of different
  formulations of the problem, among which ones suitable for the
  use of decomposition approaches;

- CapacitatedFacilityLocationBlock supports
  reformulations/relaxations of the problem via the "R3Block"
  mechanism;

- .. possibly others.
 

## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(CapacitatedFacilityLocationBlock)
target_link_libraries(<my_target> SMS++::CapacitatedFacilityLocationBlock)
```


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those
unwilling to use CMake. General instructions are:

- The arrangements of folders must be that envisioned by the
  [Umbrella SMS++ Project](https://gitlab.com/smspp/smspp-project)

- The main step is to edit the makefiles into ../extlib/. There is
  one for each of the external libraries that any module requires,
  starting with

  = [Boost](https://www.boost.org)

  = [Eigen](http://eigen.tuxfamily.org)

  = [netCDF-C++](https://www.unidata.ucar.edu/software/netcdf)

  that are required by the "core" SMS++ library and therefore by
  everyone. Setting the

```make
lib*INC = -I<paths to include files directories>
lib*LIB = -L<paths to lib files directories> -l<libs>
```

  in each allows one to set any non-standard path if the library is
  not installed in the system (or leave them empty if they are).

- The "core" SMS++ classes have a makefile for building the
  corresponding library in

```sh
SMS++/lib/makefile-lib
```

  The makefile allow to choose the compiler name and the
  optimization/debug. This builds the lib/libSMS++.a that can be
  linked upon. Also, the

```sh
SMS++/lib/makefile-inc
```

  file is provided for allowing external makefiles to ensure that
  the library is up-to-date (useful in case one is actually
  developing it). The simplest way to learn how to use it is to
  check the makefiles of the tools

```sh
CapacitatedFacilityLocationBlock/tools/makefile
```

  Note that the "basic" makefile macros

```make
CC =
SW =
```

  for setting the c++ compiler and its options are "automatically
  forwarded" from the makefile to these of the other SMS++
  components, and therefore (possibly at the cost of a make clean)
  ensure consistency during the building process.


## Tools

We provide a simple tool that converts CFL instances written in a
"non-standard" format, i.e., different from that of the
[ORLib](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/capinfo.html),
into either the "standard" format or into netCDF files.

You can run the tool from the `<build-dir>/tools` directory or
install it with the library (see above).
Run the tool without arguments for info on its usage:

```sh
cfl2nc4
```

## Getting help

If you need support, you want to submit bugs or propose a new
feature, you can
[open a new issue](https://gitlab.com/smspp/mcfblock/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our
code of conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Universita' di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General
Public License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source
license. As such, it is provided "*as is*", without any explicit or
implicit warranty that it will properly behave or it will suit your
needs. The Authors of the code cannot be considered liable, either
directly or indirectly, for any damage or loss that anybody could
suffer for having used it. More details about the non-warranty
attached to this code are available in the license description file.
