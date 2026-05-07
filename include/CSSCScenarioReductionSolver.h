/*--------------------------------------------------------------------------*/
/*-------------- File CSSCScenarioReductionSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header for CSSCScenarioReductionSolver: a separate Solver class that
 * implements the Cost-Space Scenario Clustering (CSSC) algorithm for
 * scenario reduction.
 *
 * ### Why a separate class?
 *
 * CSSC is fundamentally different from the heuristics in
 * ScenarioReductionSolver (Dupacova, BestFit, FirstFit):
 *
 *  - It requires a MILP solver (attached via BlockSolverConfig), whereas the
 *    heuristics need no external solver at all.
 *
 *  - It requires a DiscreteScenarioSet to iterate over individual scenario
 *    vectors for the N×N opportunity-cost matrix (Step 1).
 *
 *  - It requires the CapacitatedFacilityLocationBlock to be the inner block
 *    of a StochasticBlock, so that individual scenarios can be injected via
 *    StochasticBlock::set_data().
 *
 * These requirements would pollute the interface of ScenarioReductionSolver
 * if CSSC were kept there. Per the SMS++ design principle, all Solvers should
 * share the same interface; solver-specific configuration is passed through
 * the standard set_ComputeConfig() / set_*_par() methods.
 *
 * ### Configuration
 *
 * The MILP solver config and the scenario set pointer are passed via
 * set_ComputeConfig() using the "extra Configuration" field of ComputeConfig:
 *
 *   ComputeConfig cc;
 *   cc.f_extra = new SimpleConfiguration< std::pair<
 *       BlockSolverConfig * , const DiscreteScenarioSet * > >( { bsc , dss } );
 *   solver.set_ComputeConfig( &cc );
 *
 * ### Shared infrastructure
 *
 * CSSCScenarioReductionSolver inherits all output/query methods from
 * ScenarioReductionSolver (get_reduced_atoms, get_var_value, get_var_solution,
 * has_var_solution, etc.) and reuses the shared fields (nb_atoms, nb_reduced,
 * weights, ind_red, reduced_atoms, f_solution_value, f_transportation_costs).
 * Only compute() and configuration are overridden.
 *
 * \author Minh Duc Pham \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CSSCScenarioReductionSolver
#define __CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioReductionSolver.h"
#include "BlockSolverConfig.h"
#include "ThinComputeInterface.h"   // ComputeConfig
#include "DiscreteScenarioSet.h"    // full type for delete in destructor

// Forward declarations
namespace SMSpp_di_unipi_it {
 class DiscreteScenarioSet;
 class StochasticBlock;
}

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS CSSCComputeConfig ----------------------------*/
/*--------------------------------------------------------------------------*/
/** Derived ComputeConfig for CSSCScenarioReductionSolver.
 *
 * Extends the standard ComputeConfig with two CSSC-specific fields:
 *
 *  - f_extra_Configuration: (inherited) holds a BlockSolverConfig* describing
 *    which MILP solver to use for sub-problems (Step 1) and the partitioning
 *    MILP (Step 2). CSSCScenarioReductionSolver::set_ComputeConfig() reads
 *    this field and clones it.
 *
 *  - f_scenario_set: pointer to the DiscreteScenarioSet providing the N
 *    scenario vectors for compute_V_matrix(). NOT owned, caller retains
 *    ownership.
 *
 * Usage:
 * @code
 *   CSSCComputeConfig cfg;
 *   cfg.f_extra_Configuration = bsc;       // BlockSolverConfig*
 *   cfg.f_scenario_set        = dss;       // const DiscreteScenarioSet*
 *   solver.set_ComputeConfig( &cfg );
 * @endcode
 */
class CSSCComputeConfig : public ComputeConfig {
public:

 /// Pointer to the scenario set. May be owned (if loaded via deserialize/load)
 /// or non-owned (if set directly by caller). Ownership tracked by f_owned_dss.
 const DiscreteScenarioSet * f_scenario_set = nullptr;

 /// Owned DiscreteScenarioSet (created during deserialize/load, nullptr otherwise)
 DiscreteScenarioSet * f_owned_dss = nullptr;

 CSSCComputeConfig() = default;

 /// Destructor: deletes f_owned_dss if owned
 ~CSSCComputeConfig() override { delete f_owned_dss; }

 /// Clone this config
 [[nodiscard]] CSSCComputeConfig * clone( void ) const override;

 /// Serialize to netCDF: calls base class + serializes DiscreteScenarioSet
 void serialize( netCDF::NcGroup & group ) const override;

 /// Deserialize from netCDF: calls base class + reconstructs DiscreteScenarioSet
 void deserialize( const netCDF::NcGroup & group ) override;

 /// Load from txt stream: calls base class + loads DiscreteScenarioSet
 void load( std::istream & input ) override;
};

/*--------------------------------------------------------------------------*/
/*---------------- CLASS CSSCScenarioReductionSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/** Solver implementing the CSSC scenario reduction algorithm.
 *
 * Inherits all output/query methods from ScenarioReductionSolver and
 * overrides compute() with the two-step CSSC procedure:
 *
 *   Step 1 for compute_V_matrix(): solve N deterministic CFL sub-problems
 *             to build the N×N opportunity-cost matrix V.
 *
 *   Step 2 for solve_cssc_milp(V): solve a MILP partitioning problem
 *             (equations 24-29 of Keutchayan et al. 2023) to select the
 *             K best representative scenarios.
 *
 * Reference: Keutchayan, Ortmann & Rei, Computational Management Science,
 * 2023, Section 4.3.
 */
class CSSCScenarioReductionSolver : public ScenarioReductionSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// Default constructor
 CSSCScenarioReductionSolver() = default;

 /// Destructor, cleans up owned BlockSolverConfig
 ~CSSCScenarioReductionSolver() override { delete f_milp_config; }

/*--------------------------------------------------------------------------*/
/*------------------------ MAIN COMPUTATION METHOD -------------------------*/
/*--------------------------------------------------------------------------*/

 /** @brief Runs the full CSSC pipeline (Step 1 + Step 2).
  *
  * Requires that set_ComputeConfig() has been called with a ComputeConfig
  * whose f_extra field carries a BlockSolverConfig* (MILP solver) and a
  * const DiscreteScenarioSet* (scenario data). Throws std::logic_error if
  * either is missing.
  *
  * @return Solver::kOK on success.
  */
 int compute( bool changedvars = false ) override;

/*--------------------------------------------------------------------------*/
/*----------------------- CONFIGURATION METHOD -----------------------------*/
/*--------------------------------------------------------------------------*/

 /** @brief Accepts configuration via the standard SMS++ interface.
  *
  * Reads the "extra Configuration" field of the ComputeConfig to extract:
  *  - a BlockSolverConfig* for the MILP solver used in both Step 1 and Step 2
  *  - a const DiscreteScenarioSet* providing the N scenario vectors
  *
  * Expected extra type:
  *   SimpleConfiguration< std::pair< BlockSolverConfig *,
  *                                   const DiscreteScenarioSet * > >
  *
  * The BlockSolverConfig is cloned internally; the DiscreteScenarioSet is
  * not owned (caller retains ownership).
  *
  * @param cfg pointer to a ComputeConfig (may be nullptr to clear).
  */
 void set_ComputeConfig( const ComputeConfig * cfg ) override;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// MILP solver configuration (owned, cloned from what user provides).
 /** Used to attach a MILPSolver to:
  *  (a) f_Block for the N*(N-1) sub-problem solves in compute_V_matrix(), and
  *  (b) the AbstractBlock holding the CSSC MILP in solve_cssc_milp(). */
 BlockSolverConfig * f_milp_config = nullptr;

 /// Pointer to the DiscreteScenarioSet (NOT owned).
 /** Provides the N scenario vectors needed by compute_V_matrix() to inject
  *  each scenario into the CFL block via StochasticBlock::set_data(). */
 const DiscreteScenarioSet * f_scenario_set = nullptr;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /** @brief Builds the N×N opportunity-cost matrix V (CSSC Step 1).
  *
  * For each scenario i:
  *   (a) Load ξ_i via StochasticBlock::set_data() and solve the CFL → x*_i
  *   (b) For each j≠i: fix y=x*_i, load ξ_j, solve → V[i][j] = F(x*_i, ξ_j)
  *
  * @return V[i][j] = cost of first-stage solution x*_i under scenario j.
  * @throws std::logic_error  if f_Block has no StochasticBlock parent.
  * @throws std::runtime_error if any sub-problem solve fails.
  */
 std::vector< std::vector< double > > compute_V_matrix();

 /** @brief Builds and solves the CSSC MILP partitioning problem (Step 2).
  *
  * Constructs an AbstractBlock with variables x_ij, u_j, t_j and the
  * constraints from equations (24)-(28b) of the paper, then solves it
  * using f_milp_config. Populates ind_red, indices_to_choose, reduced_atoms,
  * and f_solution_value.
  *
  * @param V  the N×N opportunity-cost matrix from compute_V_matrix().
  * @throws std::runtime_error if the MILP solve fails.
  */
 void solve_cssc_milp( const std::vector< std::vector< double > > & V );

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY REGISTRATION --------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end class CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/

}  // namespace SMSpp_di_unipi_it

/*--------------------------------------------------------------------------*/

#endif  // __CSSCScenarioReductionSolver

/*--------------------------------------------------------------------------*/
/*------------ End File CSSCScenarioReductionSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/