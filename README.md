# CapacitatedFacilityLocationBlock

This project provides `CapacitatedFacilityLocationBlock`, an implementation
of the Block concept for a "pretty basic version" of the Capacitated Facility
Location (CFL) problem, a.k.a. the Capacitated Warehouse Location (CWL)
problem.

This class only represent the "basic version" of CFL and it is primarily
intended as a "didactic" implementation for showing some of the features of
SMS++, among which:

* `CapacitatedFacilityLocationBlock` supports both the unsplittable version
  of the problem, where each customer need be served by exactly one
  facility, and the splittable version where each customer can be served
  by any number of facilities.

* `CapacitatedFacilityLocationBlock` supports three different formulations of
  the problem:

  - The "natural formulation" (NF) in which the standard X[ j , i ] and
    Y[ i ] variables and the corresponding constraints are added to the
    `CapacitatedFacilityLocationBlock`. The formulation only has the two
	natural (static) groups of constraints corresponding to customer demand
	satisfaction and facility capacity + construction (aggregated), but an
	option is provided for the strong linking constraints (X[ j , i ] <=
	Y[ i ]) to be dynamically separated (the option actually applies to
	all formulations, although it makes no sense for the KF one below).

  - The Lagrange-friendly "knapsack formulation" (KF), where the
    `CapacitatedFacilityLocationBlock` "grows" m (number of facilities)
    sub-Block, each of type `BinaryKnapsackBlock` and n (number of customers)
	+ 1 variables. Sub-Block i corresponds to facility i: the first n
	variables correspond to the transportation variables X[ j , i ] between
	i and all the customers (in the natural order), while the last variable
	correspond to the design variable Y[ i ]. The first are set to be binary
	if and only if the problem is splittable, the last is always binary. The
	linking (customer satisfaction) constraints are the only static group of
	Constraint in the `CapacitatedFacilityLocationBlock`. This formulation
	has the right form so that a `LagrangianDualSolver` can be attached to
	it in order to compute tight lower bounds and the corresponding
	convexified primal solutions.

  - The Benders-friendly "flow formulation" (FF), where
    `CapacitatedFacilityLocationBlock` "grows" two sub-Block. The first one
    is an `AbstractBlock` that only has the m kBinary design variables
	Y[ i ]. The second is instead a `MCFBlock` representing the continuous
	relaxation of the problem as produced by get\_R3\_Block() (non-copy
	option), except the costs of the "facility arcs" are set to 0. The
	`CapacitatedFacilityLocationBlock` contains the linking constraints
	"arc\_flow[ i ] \leq Q[ i ] Y[ i ]", where arc\_flow[ i ] is the flow on
	the "facility arc" corresponding to facility i in the `MCFBlock`. Two
	different versions of the `MCFBlock` are supported, one in which the
	reformulation is "exact" (the problem is completely equivalent to that
	of all the other formulations), and	another that is "approximate" in
	that extra high-cost "slack arcs" are added that ensure that the
	instance is always feasible even if the aggregate demand is larger than
	the aggregate facility capacity. This is useful to apply Benders'
	decomposition to the problem without having to worry about vertical
	linearizations (feasibility cuts), which may be needed because the
	ability of current MCF solvers to produce them is severely limited. Note
	that the unsplittable version of the problem cannot be represented in
	the FF.

  The methods for loading, reading and changing the data of the instance in
  `CapacitatedFacilityLocationBlock` can all be used independently from which
  of the formulations is employed.

* `CapacitatedFacilityLocationBlock` supports reformulations/relaxations of
  the problem via the "R3Block" mechanism; in particular, besides the
  "copy" R3Block, also the "flow relaxation" of the
  CapacitatedFacilityLocationBlock is supported where the R3Block is a
  `MCFBlock` that contains the "second sub-Block of the FF" (except with the
  right costs on the "facility arcs"). Both the "exact" and "approximate"
  relaxations are supported, as described in the flow formulation.

* `CapacitatedFacilityLocationBlock` supports reading the data from three
  different text input formats (as well as from its onw netCDF one).

`CapacitatedFacilityLocationBlock` currently lacks some capabilities:

* The transportation graph is fixed and complete, there is no way to
  specify that a specific user cannot be served by a specific facility
  (save by placing a huge cost on the corresponding arc).

* Changing customers' demands via the abstract representation is not
  allowed in the SF and the KF, since the same demand is replicated in
  multiple constraints; one could ask that all the changes happen at the
  same time and that the corresponding Modification are bunched together in
  a GroupModification, which is possible but complex and not implemented
  yet. The change is instead possible in the Flow Formulation where demands
  are node deficits.

* Changing the splittable/unsplittable form of the problem, i.e., the
  integrality of all variables x[ i ][ j ], via the abstract
  representation is never allowed.

* map\_[forward/back]\_[Modification/Solution]() are fully implemented
  for both types of R3Block, except "back Modification" that is not
  implemented for the MCF R3Block.


## Getting started

These instructions will let you build `CapacitatedFacilityLocationBlock` on
your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its requirements.

- [BinaryKnapsackBlock](https://gitlab.com/smspp/binaryknapsackblock)


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

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `CapacitatedFacilityLocationBlock`, such as the
[tester for multiple features of `CapacitatedFacilityLocationBlock`](https://gitlab.com/smspp/tests/-/blob/develop/CapacitatedFacilityLocation/test.cpp?ref_type=heads),
has to include a "main makefile" of the module, which typically is either
[makefile-c](makefile-c) including all necessary libraries comprised the
"core SMS++" one, or [makefile-s](makefile-s) including all necessary
libraries but not the "core SMS++" one (for the common case in which this is
used together with other modules that already include them). These in turn
recursively include all the required other makefiles, hence one should only
need to edit the "main makefile" for compilation type (C++ compiler and its
options) and it all should be good to go. In case some of the external
libraries are not at their default location, it should only be necessary to
create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tools

We provide a simple tool that reads CFL instances written in three different
formats and convert them to the SMS++ native netCDF format supported by
CapacitatedFacilityLocationBlock, or read a netCDF file and produce the
corresponding text one.

You can run the tool from the `<build-dir>/tools` directory or install it
with the library (see above). Run the tool without arguments for info on
its usage:

```sh
txt2nc4
```

A batch file is provided to build netCDF files for a test bed composed
by three different sets of instances. First decompress `data/txt.tgz`
in place and then run `data/batch` to have the instances produced in
`data/nc4'.


## Getting help

If you need support, you want to submit bugs or propose a new feature,
you can [open a new
issue](https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our
code of conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

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
