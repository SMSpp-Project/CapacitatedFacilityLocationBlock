# CapacitatedFacilityLocationBlock

This project provides `CapacitatedFacilityLocationBlock`, an implementation
of the Block concept for a "pretty basic version" of the Capacitated Facility
Location (CFL) problem, a.k.a., the Capacitated Warehouse Location (CWL)
problem.

This class only represents the "basic version" of CFL and it is primarily
intended as a "didactic" implementation for showing some of the features of
SMS++, among which:

- `CapacitatedFacilityLocationBlock` supports both the unsplittable version
  of the problem, where each customer need be served by exactly one
  facility, and the splittable version where each customer can be served
  by any number of facilities.

- `CapacitatedFacilityLocationBlock` supports three different formulations of
  the problem:

  * The "natural formulation" (NF) in which the standard X[ j , i ] and
    Y[ i ] variables and the corresponding constraints are added to the
    `CapacitatedFacilityLocationBlock`. The formulation only has the two
    natural (static) groups of constraints corresponding to customer demand
    satisfaction and facility capacity + construction (aggregated), but an
    option is provided for the strong linking constraints (X[ j , i ] <=
    Y[ i ]) to be dynamically separated (the option actually applies to
    all formulations, although it makes no sense for the KF one below).

  * The Lagrange-friendly "knapsack formulation" (KF), where the
    `CapacitatedFacilityLocationBlock` "grows" m (number of facilities)
    sub-Block, each of type `BinaryKnapsackBlock` and n (number of customers)
    + 1 variables. Sub-Block i corresponds to facility i: the first n
    variables correspond to the transportation variables X[ j , i ] between
    i and all the customers (in the natural order), while the last variable
    corresponds to the design variable Y[ i ]. The first are set to be binary
    if and only if the problem is splittable, the last is always binary. The
    linking (customer satisfaction) constraints are the only static group of
    Constraint in the `CapacitatedFacilityLocationBlock`. This formulation
    has the right form so that a `LagrangianDualSolver` can be attached to
    it in order to compute tight lower bounds and the corresponding
    convexified primal solutions.

  * The Benders-friendly "flow formulation" (FF), where
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
    of all the other formulations), and another that is "approximate" in
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

- `CapacitatedFacilityLocationBlock` supports reformulations/relaxations of
  the problem via the "R3Block" mechanism; in particular, besides the
  "copy" R3Block, also the "flow relaxation" of the
  CapacitatedFacilityLocationBlock is supported where the R3Block is a
  `MCFBlock` that contains the "second sub-Block of the FF" (except with the
  right costs on the "facility arcs"). Both the "exact" and "approximate"
  relaxations are supported, as described in the flow formulation.

- `CapacitatedFacilityLocationBlock` supports reading the data from three
  different text input formats (as well as from its own netCDF one).

However, `CapacitatedFacilityLocationBlock` is also useful to represent
*scenario reduction* problems, whereby one wants to choose a (small)
subset (facilities) of the (many) available scenarios (customers) so as to
minimize the expected loss of accuracy in the corresponding stochastic
optimization problem. Basically, a set of scenarios is replaced by a single
one, assigning it the total probability of the original ones. For this
application, all three formulations above can optionally include a
constraint on the maximum number of facilities that can be opened (when
wc & 4 in generate\_abstract\_constraints()). Furthermore, the
`ScenarioReductionSolver` is provided that implements a bunch of fast
heuristics for the specific version of CFL used in scenario reduction
applications.

Still, `CapacitatedFacilityLocationBlock` currently lacks some capabilities:

- The transportation graph is fixed and complete, there is no way to
  specify that a specific user cannot be served by a specific facility
  (save by placing a huge cost on the corresponding arc).

- Changing customers' demands via the abstract representation is not
  allowed in the NF and the KF, since the same demand is replicated in
  multiple constraints; one could ask that all the changes happen at the
  same time and that the corresponding Modification are bunched together in
  a GroupModification, which is possible but complex and not implemented
  yet. The change is instead possible in the Flow Formulation where demands
  are node deficits.

- Changing the splittable/unsplittable form of the problem, i.e., the
  integrality of all variables x[ i ][ j ], via the abstract
  representation is never allowed.

- map\_[forward/back]\_[Modification/Solution]() are fully implemented
  for both types of R3Block, except "back Modification" that is not
  implemented for the MCF R3Block.


### ScenarioReductionSolver

`ScenarioReductionSolver` is a specialized solver for scenario reduction
problems formulated as CFL instances, i.e., such that

- all facility capacities must equal 1.0 (allowing full probability mass
  assignment to any selected scenario)

- customer demands represent scenario probabilities (automatically
  normalized if they don't sum to 1.0)

- number of facilities must equal number of customers (square distance
  matrix)

- transportation costs represent pairwise scenario distances

It implements algorithms to select a representative subset of scenarios
that aims to minimize a Wasserstein distance between the full scenario set
and the chosen representative subset.

**Available algorithms:**

- **Baseline**: Simple greedy selection based on scenario probabilities

- **Dupacova**: Forward selection algorithm that iteratively adds scenarios
  to minimize Wasserstein distance (default).

- **BestFit**: Local search with best improvement selection among all
  possible pairs of choices.

- **FirstFit**: Local search with first improvement selection.

The solver can be configured through parameters:

- `intAlgorithm`: Algorithm selection (0=Baseline, 1=Dupacova, 2=BestFit,
  3=FirstFit)

- `dblEll`: Power parameter for Wasserstein distance (default: 2.0)

- `intShuffle`: Enable shuffling for FirstFit algorithm (currently only
  implemented for FirstFit, but could be extended to BestFit)

- `intUseWarmstart`: Enable warm start for local search algorithms


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
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(CapacitatedFacilityLocationBlock)
target_link_libraries(<my_target> SMS++::CapacitatedFacilityLocationBlock)
```


### Running the tests with CMake

A unit test will be built with the library; set `BUILD_TESTING` to `OFF` to
disable it. It exercises the `ScenarioReductionSolver` on
`CapacitatedFacilityLocationBlock` instances.


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

## Data

We provide some data sets that are used, among other things, by some of the
testers of the [test repo](https://gitlab.com/smspp/tests). Since they are
large they are not included in the repo. They are automatically downloaded
by Cmake if the test repo is included, but if you are not using Cmake to
build the system you need to do it by hand, via

```sh
cd data
wget https://gitlab.com/api/v4/projects/24780184/packages/generic/txt/latest/txt.tgz
tar xzvf txt.tgz
```

This builds the following folders:

- [data/txt/ORLib](data/txt/ORLib) that contains the [original
  ORLib instances](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/capinfo.html)
  (with an addition and minor tweaks) in the ORLib format, see
  [data/txt/ORLib/doc](data/txt/ORLib/doc) for details

- [data/txt/TBED](data/txt/TBED) that contains instances that have
  been downloaded from
  [here](https://or-brescia.unibs.it/instances/instances_sscflp) and are in
  facility-oriented, demands-first format: see
  [data/txt/TBED/Reame.txt](data/txt/TBED/Reame.txt) for details

- [data/txt/Yang](data/txt/Yang) that contains instances that have
  been downloaded from
  [here](https://or-brescia.unibs.it/instances/instances_sscflp) and are in
  facility-oriented, demands-last format: see
  [data/txt/Yang/format.pdf](data/txt/Yang/format.pdf) for details


## Tools

We provide a simple tool that reads CFL instances written in the three
different formats (see `Data` above) and converts them to the SMS++ native
netCDF format supported by `CapacitatedFacilityLocationBlock`, or reads a
netCDF file and produces the corresponding text one.

You can run the tool from the `<build-dir>/tools` directory or install it
with the library (see above). Run the tool without arguments for info on
its usage:

```sh
cfl2nc4
```

A batch file is provided to build netCDF files for a test bed composed
by three different sets of instances. First obtain (see `Data` above) and
decompress `data/txt.tgz` in place, then run `data/batch` to have the
instances produced in `data/nc4`.


## Tests

The [test](test) folder contains a tester for the `ScenarioReductionSolver`,
exercising its parameter management and solution handling on
`CapacitatedFacilityLocationBlock` instances. The integration suites that
compare the block against a `:MILPSolver` live in the
[tests repo](https://gitlab.com/smspp/tests).


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

- **Benoit Tran**  
  Dipartimento di Informatica  
  Università di Pisa

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
