# MCFBlock / MCFSolver

This project covers two conceptually different things (which may one day be
split to two different projects):

- `MCFBlock`, a SMS++ :Block for Linear Min-Cost Flow Problems

- `MCFSolver`, a SMS++ :Solver for MCFBlock based on forwarding the interface
  of (objects derived from the general abstract) [MCFClass of the
  MCFClass project](http://www.di.unipi.it/optimize/Software/MCF.html)


## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

- [MCFClass](https://github.com/frangio68/Min-Cost-Flow-Class) and its
  requirements (depending on the actual :MCFClass solvers built).


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
Moreover, you can choose the solver with the variable `MCFBlock_SOLVER`.
Available values are:

| Value     | Solver             |
| --------- | ------------------ |
| `relax`   | RelaxIV            |
| `cplex`   | CPLEX (default)    |
| `mfsmx`   | Simplex            |
| `sptree`  | Shortest Path Tree |

Optionally, install the library in the system with:

```sh
sudo make install
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(MCFBlock)
target_link_libraries(<my_target> SMS++::MCFBlock)
```


### Running the tests with CMake

A unit test will be built with the library.
To disable it, set the option `BUILD_TESTING` to `OFF`.

The test takes an instance of a MCF in DIMACS or NC4 format. The MCF problem
is then repeatedly solved with several changes in costs/capacities/deficits,
arcs openings/closures and arcs additions/deletions. The same operations are
performed on the two solvers, and the results are compared.

The test can be run manually, using the provided batch file,
or using `ctest` from the build directory.


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. General instructions are:

- The arrangements of folders must be that envisioned by the
  [Umbrella SMS++ Project](https://gitlab.com/smspp/smspp-project)

- The main step is to edit the makefiles into ../extlib/. There is one for
  each of the external libraries that any module requires, starting with

  = [Boost](https://www.boost.org)

  = [Eigen](http://eigen.tuxfamily.org)

  = [netCDF-C++](https://www.unidata.ucar.edu/software/netcdf)

  that are required by the "core" SMS++ library and therefore by everyone.
  Setting the

```make
lib*INC = -I<paths to include files directories>
lib*LIB = -L<paths to lib files directories> -l<libs>
```

  in each allows one to set any non-standard path if the library is not
  installed in the system (or leave them empty if they are).

- The "core" SMS++ classes have a makefile for building the corresponding
  library in

```sh
SMS++/lib/makefile-lib
```

  The makefile allow to choose the compiler name and the optimization/debug.
  This builds the lib/libSMS++.a that can be linked upon. Also, the

```sh
SMS++/lib/makefile-inc
```

  file is provided for allowing external makefiles to ensure that the library
  is up-to-date (useful in case one is actually developing it). The simplest
  way to learn how to use it is to check the makefiles of the tester

```sh
MCFBlock/test/makefile
```

  Note that the "basic" makefile macros

```make
CC =
SW =
```

  for setting the c++ compiler and its options are "automatically forwarded"
  from the makefile to these of the other SMS++ components, and therefore
  (possibly at the cost of a make clean) ensure consistency during the
  building process.

- The [MCFClass project](https://github.com/frangio68/Min-Cost-Flow-Class)
  has a similar arrangement with its own extlib/ folder that must be
  independently edited in an analogous way.


## Tools

We provide a simple tool that converts MCF instances written in the DIMACS
standard into netCDF files. Optionally it hacks into the netCDF file to
change the number of static and dynamic nodes and arcs, as well as the
maximum number of nodes and arcs.

You can run the tool from the `<build-dir>/tools` directory or install it
with the library (see above).
Run the tool without arguments for info on its usage:

```sh
dmx2nc4
```

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/mcfblock/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
