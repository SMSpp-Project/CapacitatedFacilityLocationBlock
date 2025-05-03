/*--------------------------------------------------------------------------*/
/*------------------ File ScenarioReductionSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ScenarioReductionSolver class, which solves
 * the Discrete Scenario Reduction problem formulated as a Capacitated 
 * Facility Location Problem.
 * 
 * It assumes that the user loaded a CapacitatedFacilityLocationBlock 
 * instance and uses its internal variables and methods. 
 * ScenarioReductionSolver checks that the instance can be correctly 
 * interpreted as a MILP reformulation of the (Discrete) Scenario Reduction
 * Problem. 
 * 
 * For now, this solver works directly with the physical representation of 
 * the problem.
 * 
 * This implementation includes multiple scenario reduction algorithms:
 * - Dupacova's forward algorithm (default)
 * - Local search with BestFit strategy
 * - Local search with FirstFit strategy
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ScenarioReductionSolver
 #define __ScenarioReductionSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <random>
#include <unordered_map>

#include "Solver.h"
#include "CapacitatedFacilityLocationBlock.h"
#include "DiscreteScenarioSet.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioReductionSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * A solver for Capacitated Facility Location (CFL) problems that works directly
 * with the physical representation of the problem, forgetting the abstract
 * representation for now. This solver solves the Scenario Reduction Problem
 * as a specific instance of CFL.
 * 
 * Multiple scenario reduction algorithms are implemented:
 * - Dupacova's forward algorithm (default)
 * - Local search with BestFit strategy 
 * - Local search with FirstFit strategy
 */
class ScenarioReductionSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

  // Using types from CapacitatedFacilityLocationBlock for clarity
  using Index = CapacitatedFacilityLocationBlock::Index;
  using Demand = CapacitatedFacilityLocationBlock::Demand;
  using Cost = CapacitatedFacilityLocationBlock::Cost;
  using DVector = CapacitatedFacilityLocationBlock::DVector;
  using CVector = CapacitatedFacilityLocationBlock::CVector;
  using CMatrix = CapacitatedFacilityLocationBlock::CMatrix;
  using IntSolution = CapacitatedFacilityLocationBlock::IntSolution;
  using CntSolution = CapacitatedFacilityLocationBlock::CntSolution;

  using OFValue = RealObjective::OFValue;
  using ScenarioIndex = ScenarioGenerator::ScenarioIndex;

  /**
   * Enumeration of available scenario reduction algorithms
   */
  enum class Algorithm {
    Dupacova,  // Dupacova's forward algorithm (default)
    BestFit,   // Local search with BestFit strategy
    FirstFit   // Local search with FirstFit strategy
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

  /// Default constructor
  ScenarioReductionSolver();

  /// Destructor
  ~ScenarioReductionSolver() override;

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  /// Sets the Block that the Solver has to solve and caches physical data
  void set_Block(Block* block) override;

  /// Helper method to refresh cached data from the block
  void refresh_cached_data(int k);

  /// Solve the scenario reduction problem from a choice of methods 
  virtual int compute(bool changedvars = false) override;

  /// Get the power ell in the objective ell-Wasserstein distance
  [[nodiscard]] float get_ell() const { return ell; }

  /// Returns the value of the current solution, if any
  [[nodiscard]] OFValue get_var_value() override;

  /// Writes the current solution in the Block
  void get_var_solution(Configuration* solc = nullptr) override;

  const IntSolution& get_reduced_atoms() const { return reduced_atoms; }

  void set_ell(float exponent) { ell = exponent; }

  /**
   * Set the algorithm to use for scenario reduction
   * 
   * @param alg The algorithm to use
   */
  void set_algorithm(Algorithm alg) { algorithm = alg; }

  /**
   * Get the current algorithm being used
   * 
   * @return The current algorithm
   */
  Algorithm get_algorithm() const { return algorithm; }

  /**
   * Set the rho parameter for local search algorithms (minimum improvement threshold)
   * 
   * @param value The rho value (should be non-negative)
   */
  void set_rho(double value) { 
    if (value < 0.0) {
      throw std::invalid_argument("rho should be non-negative");
    }
    rho = value; 
  }

  /**
   * Get the rho parameter
   * 
   * @return The current rho value
   */
  double get_rho() const { return rho; }

  /**
   * Set whether to use shuffling in FirstFit local search
   * 
   * @param enable Whether to enable shuffling
   */
  void set_shuffle(bool enable) { shuffle = enable; }

  /**
   * Get whether shuffling is enabled for FirstFit
   * 
   * @return True if shuffling is enabled
   */
  bool get_shuffle() const { return shuffle; }

  /**
   * Set random seed for shuffling
   * 
   * @param seed The random seed
   */
  void set_random_seed(unsigned int seed) { rng.seed(seed); }

  /// Selecting a pair of atoms to swap
  /**
   * Virtual method which implements Dupacova greedy heuristic by default.
   * Abstract method to select a pair of indices (i, j) for swapping atoms 
   * in a candidate reduced distribution.
   * This method should be implemented in concrete subclasses to define 
   * the strategy for selecting a pair of atoms to swap.
   */
  virtual std::tuple<int, int> pick_candidate(const std::vector<bool>& ind, 
    const std::vector<double>& min_cost);

  const ScenarioIndex& get_nb_atoms() const { return nb_atoms; }
  const ScenarioIndex& get_nb_reduced() const { return nb_reduced; }
  const DVector* get_weights() const { return weights; }

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

  /**
   * Pair of Index and value to store distances
   */
  using IndexDistancePair = std::pair<Index, double>;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

  /// Algorithm selection
  Algorithm algorithm = Algorithm::Dupacova;

  /// Parameters for local search
  double rho = 0.0;         // Minimum improvement threshold
  bool shuffle = false;     // Whether to shuffle indices (for FirstFit)
  std::mt19937 rng{std::random_device{}()};  // Random number generator
  double dist_dupa = std::numeric_limits<double>::infinity(); // Distance from Dupacova algorithm

  /// Problem dimensions
  ScenarioIndex nb_atoms;
  ScenarioIndex nb_reduced;

  std::vector<ScenarioIndex> indices_to_choose;  // Indices not yet chosen
  std::vector<Index> ind_red;                   // Indices in reduced set (for local search)

  // Power in the ell-Wasserstein distance
  float ell = 2.0;

  /// Problem data - these are references to the data in the block
  const DVector* weights;
  const CMatrix* f_transportation_costs;

  /// Solution data
  IntSolution reduced_atoms;   // Binary values indicating which atoms are selected
  double f_solution_value;     // Objective function value of current solution

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

  /**
   * Perform scenario reduction using Dupacova's forward algorithm
   * 
   * @return Status code (kOK if successful)
   */
  int compute_dupacova();

  /**
   * Perform scenario reduction using local search (BestFit or FirstFit)
   * 
   * @return Status code (kOK if successful)
   */
  int compute_local_search();

  /**
   * Initialize the reduced set of atoms for local search
   * 
   * @return The distance from Dupacova if rho > 0
   */
  double init_local_search();

  /**
   * Determine whether a candidate update improves the distance sufficiently
   * 
   * @param trial_d The candidate distance
   * @param curr_d The current distance
   * @param dist_dupa The distance from Dupacova's algorithm
   * @return True if the improvement is sufficient
   */
  bool improvement_condition(double trial_d, double curr_d, double dist_dupa) const;

  /**
   * Select a pair of atoms to swap using BestFit strategy
   * 
   * @param curr_d Current distance
   * @return A tuple containing indices to swap and resulting distance
   */
  std::tuple<Index, Index, double> pick_ij_bestfit(double& curr_d);

  /**
   * Helper method for BestFit to select the best atom to add
   * 
   * @param curr_indices Current indices in the reduced set
   * @return A pair with the selected index and distance
   */
  IndexDistancePair bestfit_selection(const std::vector<Index>& curr_indices);

  /**
   * Select a pair of atoms to swap using FirstFit strategy
   * 
   * @param curr_d Current distance
   * @return A tuple containing indices to swap and resulting distance
   */
  std::tuple<Index, Index, double> pick_ij_firstfit(double& curr_d);

  /**
   * Helper method for FirstFit to select the first satisfactory atom
   * 
   * @param curr_indices Current indices in the reduced set
   * @param curr_d Current distance
   * @return A pair with the selected index and distance
   */
  IndexDistancePair firstfit_selection(
    const std::vector<Index>& curr_indices, 
    double curr_d);

  /**
   * Swap atoms in the reduced distribution and update the indices
   * 
   * @param i Index of the atom to remove
   * @param j Index of the atom to add
   */
  void swap_indices(Index i, Index j);

  /**
   * Update the reduced_atoms binary vector from ind_red
   */
  void update_reduced_atoms();

  // /// Helper method to process modifications
  // void process_pending_modifications() 
  // {
  //   // Process all pending modifications
  //   sp_Mod mod;
  //   while ((mod = pop()) != nullptr) {
  //     // Handle different types of modifications
      
  //     // Check for changes in facility costs
  //     if (auto cflMod = std::dynamic_pointer_cast<CapacitatedFacilityLocationBlockRngdMod>(mod)) {
  //       if (cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgFCost ||
  //           cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgTCost ||
  //           cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgCap ||
  //           cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgDem) {
  //         // Data has changed, refresh our cached references
  //         refresh_cached_data();
  //       }
  //     }
  //     // Check for subset-based modifications
  //     else if (auto cflSubsetMod = std::dynamic_pointer_cast<CapacitatedFacilityLocationBlockSbstMod>(mod)) {
  //       // Similar handling as above
  //       refresh_cached_data();
  //     }
  //     // Check for nuclear option (reload)
  //     else if (std::dynamic_pointer_cast<NBModification>(mod)) {
  //       // Complete reload of the problem
  //       refresh_cached_data();
        
  //       // Reinitialize solution structures
  //       reduced_atoms.resize(nb_atoms, false);
  //       f_solution_value = 0.0;
  //     }
  //   }
  // }

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY REGISTRATION --------------------------*/
/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class ScenarioReductionSolver )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ScenarioReductionSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------- End File ScenarioReductionSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/