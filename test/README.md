# test

A tester for `ScenarioReductionSolver`, the solver that selects a
representative subset of scenarios for a `CapacitatedFacilityLocationBlock`.

`SR_test` exercises the solver end to end: parameter management (setting and
getting the algorithm, shuffle and random-seed parameters, boundary values,
validation and defaults), solution handling (facility-selection and
customer-assignment variables, warm starting, solution independence and
persistence), and a comparison of the reduction algorithms (for instance the
baseline against Dupačová's). The scenarios and the solver configuration are
read from `scenarios.txt` and `BSConfig_SR.txt`.

The `makefile` builds the executable including the
`CapacitatedFacilityLocationBlock` module and the core SMS++ library.


## Authors

- **Benoît Tran**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
