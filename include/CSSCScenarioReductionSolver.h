/*--------------------------------------------------------------------------*/
/*-------------- File CSSCScenarioReductionSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header of CSSCScenarioReductionSolver.
 *
 * ### Design
 *
 * This solver uses TWO blocks with distinct roles:
 *
 *   f_Block     (inherited from Solver, set via set_Block())
 *               Must be an N×N CapacitatedFacilityLocationBlock where:
 *                 NCustomers     = N  -> nb_atoms
 *                 NMaxFacilities = K  -> nb_reduced
 *                 demands[i]     = scenario weights
 *                 tcosts[i][j]   = pairwise distances (for Wasserstein)
 *               Used only by the base class ScenarioReductionSolver to
 *               initialise nb_atoms, nb_reduced, weights,
 *               f_transportation_costs.
 *
 *   f_sub_block (set via set_sub_problem_block())
 *               The real two-stage CFL block (nf facilities × nc customers).
 *               Must have a StochasticBlock as its parent (get_f_Block())
 *               so that scenarios can be injected via set_data().
 *               compute_V_matrix() attaches the MILP solver here, solves
 *               N^2 sub-problems, and reads/fixes/restores the y variables.
 *
 * \author Minh Duc Pham \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/

#ifndef __CSSCScenarioReductionSolver
#define __CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioReductionSolver.h"
#include "BlockSolverConfig.h"
#include "DiscreteScenarioSet.h"
#include "ThinComputeInterface.h"    // ComputeConfig

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS CSSCComputeConfig -----------------------------*/
/*--------------------------------------------------------------------------*/
/** ComputeConfig for CSSCScenarioReductionSolver.
 *
 * Carries:
 *   f_extra_Configuration: BlockSolverConfig for Step 1 sub-problems
 *                          (LP relaxation recommended for speed).
 *                          Also used for Step 2 if f_milp_config is null.
 *
 *   f_milp_config: Optional BlockSolverConfig for Step 2 MILP. If null, f_extra_Configuration is used instead.
 *
 *   f_scenario_set: pointer to the DiscreteScenarioSet. Not owned.
 */

class CSSCComputeConfig : public ComputeConfig {

public:

 CSSCComputeConfig() = default;

 ~CSSCComputeConfig() override {
  delete f_owned_dss;
  delete f_milp_config;
 }

 void serialize( netCDF::NcGroup & group ) const override;
 void deserialize( const netCDF::NcGroup & group ) override;
 void load( std::istream & input ) override;

 CSSCComputeConfig * clone() const override;

 /// BlockSolverConfig for Step 1 sub-problems (LP relaxation recommended).
 /// Stored in f_extra_Configuration (inherited from ComputeConfig).

 /// Optional separate BlockSolverConfig for Step 2 partitioning MILP.
 /// If null, f_extra_Configuration is used for Step 2 as well.
 BlockSolverConfig * f_milp_config = nullptr;

 /// Pointer to the DiscreteScenarioSet (not owned unless f_owned_dss != null)
 const DiscreteScenarioSet * f_scenario_set = nullptr;

 /// If the DSS was constructed internally (load/deserialize), we own it
 DiscreteScenarioSet * f_owned_dss = nullptr;

};  // end class CSSCComputeConfig

/*--------------------------------------------------------------------------*/
/*---------------- CLASS CSSCScenarioReductionSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/** Scenario reduction using the Cost-Space Scenario Clustering (CSSC)
 *  algorithm of Keutchayan, Ortmann & Rei (2023).
 *
 * Usage:
 *   1. cssc.set_Block( sr_cfl ): N×N block for base class
 *   2. cssc.set_sub_problem_block( b ): real CFL block for Step 1
 *   3. cssc.set_ComputeConfig( cfg ): CSSCComputeConfig with MILP solver
 *                                         + DiscreteScenarioSet
 *   4. cssc.compute()
 *   5. cssc.get_ind_red() / get_reduced_atoms()
 */

class CSSCScenarioReductionSolver : public ScenarioReductionSolver {

public:

 CSSCScenarioReductionSolver() = default;

 ~CSSCScenarioReductionSolver() override {
  delete f_milp_config;
  delete f_milp_config_step2;
 }

 /*-----------------------------------------------------------------------*/
 /** Set the real two-stage CFL block used for sub-problem solves in
  *  Step 1 (compute_V_matrix).
  *
  *  Requirements:
  *    - block must be a CapacitatedFacilityLocationBlock
  *    - block->get_f_Block() must return a StochasticBlock that wraps it,
  *      with a DataMapping for scenario injection via set_data()
  *
  *  This is separate from set_Block() which receives the N×N synthetic
  *  block used only to communicate N and K to the base class.
  */
 void set_sub_problem_block( Block * block ) {
  auto * cfl = dynamic_cast< CapacitatedFacilityLocationBlock * >( block );
  if( ! cfl )
   throw std::invalid_argument(
     "CSSCScenarioReductionSolver::set_sub_problem_block: "
     "block must be a CapacitatedFacilityLocationBlock." );
  f_sub_block = cfl;
 }

 /*-----------------------------------------------------------------------*/

 int compute( bool changedvars = false ) override;

 void set_ComputeConfig( const ComputeConfig * cfg ) override;

/*--------------------------------------------------------------------------*/

protected:

 /// The real CFL block used for sub-problem solves in Step 1.
 /// Not owned. Set via set_sub_problem_block().
 CapacitatedFacilityLocationBlock * f_sub_block = nullptr;

 /// Non-owning pointer to the scenario set. Set via CSSCComputeConfig.
 const DiscreteScenarioSet * f_scenario_set = nullptr;

private:

 /// BlockSolverConfig for Step 1 sub-problems. Owned.
 BlockSolverConfig * f_milp_config = nullptr;

 /// BlockSolverConfig for Step 2 MILP. Owned.
 /// If null, f_milp_config is used for Step 2 as well.
 BlockSolverConfig * f_milp_config_step2 = nullptr;

 /*-----------------------------------------------------------------------*/

 /// Step 1: build the N×N opportunity-cost matrix V.
 std::vector< std::vector< double > > compute_V_matrix();

 /// Step 2: solve the CSSC MILP partitioning problem.
 void solve_cssc_milp( const std::vector< std::vector< double > > & V );

public:

 /// Cluster assignment from Step 2 MILP: f_scenario_assignment[i] = index of
 /// the representative that scenario i is assigned to.
 /// Populated by solve_cssc_milp() from the x_ij variables before the
 /// internal AbstractBlock is destroyed.
 std::vector< Index > f_scenario_assignment;

 /*-----------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

};  // end class CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/

}  // namespace SMSpp_di_unipi_it

#endif  //CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/
/*----------- End File CSSCScenarioReductionSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/