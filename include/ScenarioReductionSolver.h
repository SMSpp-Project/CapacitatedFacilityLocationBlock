/*--------------------------------------------------------------------------*/
/*------------------ File ScenarioReductionSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ScenarioReductionSolver class, a specialized Solver for
 * the Discrete Scenario Reduction problem formulated as a Capacitated Facility
 * Location problem.
 *
 * The ScenarioReductionSolver class provides efficient heuristic algorithms
 * for scenario reduction in stochastic optimization. It interprets a
 * CapacitatedFacilityLocationBlock instance as a scenario reduction problem,
 * where facilities represent selected scenarios and customers represent all
 * scenarios, with transportation costs encoding the distance between scenarios.
 *
 * The solver works directly with the physical representation of the Block and
 * implements several scenario reduction algorithms including Dupačová's forward
 * selection, baseline selection by probability weight, and local search methods
 * with BestFit and FirstFit strategies.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Benoît Tran
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

#include "Solver.h"
#include "CapacitatedFacilityLocationBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioReductionSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// specialized Solver for scenario reduction via heuristic algorithms
/** The ScenarioReductionSolver class is a specialized Solver that interprets
 * a CapacitatedFacilityLocationBlock as a Discrete Scenario Reduction problem
 * and solves it using efficient heuristic algorithms.
 *
 * ### Problem Interpretation
 *
 * The solver maps the CFL problem to scenario reduction as follows:
 * - Facilities represent selected scenarios (reduced set)
 * - Customers represent all scenarios (original set)
 * - Transportation costs encode distances (already raised to power ell) between scenarios
 * - The objective minimizes the \f$\ell\f$-th power of the Wasserstein distance between
 *   the original and reduced probability distributions
 *
 * ### Requirements
 *
 * The CapacitatedFacilityLocationBlock must satisfy:
 * - All facility capacities must equal 1.0 (allowing full probability mass
 *   assignment to any selected scenario)
 * - Customer demands represent scenario probabilities (automatically normalized
 *   if they don't sum to 1.0)
 * - Number of facilities must equal number of customers (square distance matrix)
 * - Transportation costs represent pairwise scenario distances
 *
 * ### Implemented Algorithms
 *
 * The solver provides four heuristic algorithms:
 * - **Baseline**: Selects the k scenarios with highest probability weights
 * - **Dupačová** (default): Forward selection algorithm minimizing Wasserstein
 *   distance at each step
 * - **BestFit**: Local search exploring all possible swaps and selecting the
 *   best improvement
 * - **FirstFit**: Local search accepting the first improvement found (with
 *   optional shuffling)
 *
 * ### Thread Safety
 *
 * The solver is thread-safe and uses RAII locking for all public methods.
 * Multiple threads can safely call compute() and parameter methods concurrently.
 *
 * @note For exact MILP-based scenario reduction, attach a MILPSolver to the
 *       CapacitatedFacilityLocationBlock instead of using this heuristic solver.
 */
class ScenarioReductionSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

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

  /**
   * Enumeration of available scenario reduction algorithms
   */
  enum class Algorithm {
    Baseline,  // Select scenarios with the most pb. weights
    Dupacova,  // Dupacova's forward algorithm (default)
    BestFit,   // Local search with BestFit strategy
    FirstFit   // Local search with FirstFit strategy
  };

  /// public enum extending int_par_type_S for ScenarioReductionSolver
  enum int_par_type_SRS {
    intAlgorithm = intLastAlgPar,      ///< Algorithm selection (0=Baseline, 1=Dupacova, 2=BestFit, 3=FirstFit)
    intShuffle = intLastAlgPar + 1,    ///< Enable shuffling for FirstFit (0=false, 1=true)
    intRandomSeed = intLastAlgPar + 2, ///< Random seed for shuffling
    intUseWarmstart = intLastAlgPar + 3, ///< Enable warm start for local search (0=false, 1=true)
    intLastParSRS                      ///< First allowed parameter for derived classes
  };

  /// public enum extending dbl_par_type_S for ScenarioReductionSolver
  enum dbl_par_type_SRS {
    dblRho = dblLastAlgPar,            ///< Minimum improvement threshold for local search
    dblLastParSRS = dblLastAlgPar + 1  ///< First allowed parameter for derived classes
  };

  /// public enum extending vint_par_type_S for ScenarioReductionSolver
  enum vint_par_type_SRS {
    vintWarmstartIndices = vintLastAlgPar, ///< Custom warm start indices (empty = use Dupacova)
    vintLastParSRS                         ///< First allowed parameter for derived classes
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioReductionSolver
 *  @{ */

  /// constructs a ScenarioReductionSolver with default parameters
  ScenarioReductionSolver();

  /// destructor, virtual as required for base classes
  ~ScenarioReductionSolver() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods to configure ScenarioReductionSolver
 *  @{ */


  /** @brief Sets the Block to be solved and initializes internal data structures.
   *
   * This method attaches the solver to a CapacitatedFacilityLocationBlock and
   * validates that it meets the requirements for scenario reduction. It also
   * generates the abstract variables needed for solution writing and caches
   * the physical data for efficient access.
   *
   * @param block pointer to a CapacitatedFacilityLocationBlock instance
   *
   * @throws std::invalid_argument if block is not a CapacitatedFacilityLocationBlock
   * @throws std::invalid_argument if any capacity != 1.0 when k > 0
   * @throws std::logic_error if number of facilities != number of customers
   * @throws std::runtime_error if unable to lock the Block
   */
  void set_Block(Block* block) override;

  /** @brief Refreshes cached data from the Block and validates problem structure.
   *
   * This method is called internally by set_Block() and validates that the
   * Block data satisfies all requirements for scenario reduction. It caches
   * frequently accessed data and normalizes scenario probabilities if needed.
   *
   * @param k number of scenarios to select (0 <= k <= total scenarios)
   *
   * @throws std::invalid_argument if k < 0 or k > number of scenarios
   * @throws std::invalid_argument if any capacity != 1.0 when k > 0
   * @throws std::logic_error if Block is not set or not a CFL Block
   * @throws std::logic_error if number of facilities != number of customers
   * @throws std::runtime_error if Block has missing or invalid data
   */
  void refresh_cached_data(int k);

  /** @brief Solves the scenario reduction problem using the selected algorithm.
   *
   * Executes the currently selected scenario reduction algorithm (set via
   * intAlgorithm parameter) and stores the solution in internal data structures.
   * The solution can then be retrieved via get_var_solution().
   *
   * @param changedvars ignored for this solver (required by base class)
   * @return Solver::kOK if successful, Solver::kError otherwise
   */ 
  virtual int compute(bool changedvars = false) override;

  /** @brief Checks if a solution is available.
   * @return true if compute() has been successfully called, false otherwise
   */
  [[nodiscard]] bool has_var_solution() override;

  /** @brief Returns the objective value of the current solution.
   * @return the Wasserstein distance of the reduced scenario set
   * @throws std::logic_error if no solution is available
   */
  [[nodiscard]] OFValue get_var_value() override;

  /** @brief Writes the solution to the Block's y and x variables.
   *
   * Sets the y[i] variables in the Block to 1.0 for selected scenarios
   * (representatives) and 0.0 for unselected ones.
   * 
   * Also sets the x[i][j] variables to represent optimal assignments:
   * x[i][j] = 1.0 if customer j is assigned to facility i (based on 
   * minimum transportation cost), 0.0 otherwise. This assignment information
   * is essential for proper weight aggregation in scenario reduction.
   *
   * @param solc ignored configuration parameter (for base class compatibility)
   * @throws std::logic_error if no Block is set or no solution is available
   */
  void get_var_solution(Configuration* solc = nullptr) override;

  /** @brief Returns the binary solution indicating selected scenarios.
   * @return vector where reduced_atoms[i] = true if scenario i is selected
   */
  const IntSolution& get_reduced_atoms() const { return reduced_atoms; }

  /** @brief Returns the total number of scenarios.
   * @return number of scenarios in the original set
   */
  const Index& get_nb_atoms() const { return nb_atoms; }

  /** @brief Returns the number of scenarios to select.
   * @return k, the size of the reduced scenario set
   */
  const Index& get_nb_reduced() const { return nb_reduced; }

  /** @brief Returns pointer to the scenario probability weights.
   * @return pointer to normalized probability vector (may be nullptr if not set)
   */
  const DVector* get_weights() const { return weights; }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR HANDLING THE PARAMETERS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the ScenarioReductionSolver
 * @{ */

  /** @brief Sets an integer parameter.
   *
   * @param par parameter identifier (intAlgorithm, intShuffle, intRandomSeed,
   *            intUseWarmstart)
   * @param value the parameter value
   * @throws std::invalid_argument if value is out of valid range
   */
  void set_par(idx_type par, int value) override;

  /** @brief Sets a double parameter.
   *
   * @param par parameter identifier (dblRho)
   * @param value the parameter value
   * @throws std::invalid_argument if value is out of valid range
   */
  void set_par(idx_type par, double value) override;

  /** @brief Sets a vector of integers parameter.
   *
   * @param par parameter identifier (vintWarmstartIndices)
   * @param value vector of scenario indices for warm start
   */
  void set_par(idx_type par, std::vector<int>&& value) override;

  /** @brief Returns an integer parameter value.
   * @param par parameter identifier
   * @return the current value of the parameter
   */
  int get_int_par(idx_type par) const override;

  /** @brief Returns a double parameter value.
   * @param par parameter identifier
   * @return the current value of the parameter
   */
  double get_dbl_par(idx_type par) const override;

  /** @brief Returns a vector of integers parameter.
   * @param par parameter identifier
   * @return reference to the parameter vector
   */
  const std::vector<int>& get_vint_par(idx_type par) const override;
  
  /** @brief Returns the default value for an integer parameter.
   * @param par parameter identifier
   * @return the default value
   */
  int get_dflt_int_par(idx_type par) const override;

  /** @brief Returns the default value for a double parameter.
   * @param par parameter identifier
   * @return the default value
   */
  double get_dflt_dbl_par(idx_type par) const override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

  /// pair type for storing (scenario index, distance) in algorithms
  using IndexDistancePair = std::pair<Index, double>;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private data members
 * @{ */

  Algorithm algorithm = Algorithm::Dupacova;  ///< selected reduction algorithm

  double rho = 0.0;            ///< minimum improvement threshold for local search
  bool shuffle = false;        ///< whether to shuffle indices in FirstFit
  bool use_warmstart = false;  ///< whether to use warm start for local search
  std::vector<Index> warmstart_indices;  ///< custom warm start indices
  std::mt19937 rng{std::random_device{}()};  ///< random number generator
  double dist_dupa = std::numeric_limits<double>::infinity();  ///< Dupacova distance cache

  Index nb_atoms;     ///< total number of scenarios
  Index nb_reduced;   ///< number of scenarios to select (k)

  std::vector<Index> indices_to_choose;  ///< scenarios not yet selected
  std::vector<Index> ind_red;            ///< indices in reduced set

  const DVector* weights;                 ///< pointer to scenario probabilities
  const CMatrix* f_transportation_costs;  ///< pointer to distance matrix
  
  std::unique_ptr<DVector> normalized_weights;  ///< normalized probabilities if needed

  IntSolution reduced_atoms;   ///< binary solution (true = selected)
  double f_solution_value;     ///< objective value (Wasserstein distance)

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private methods for scenario reduction algorithms
 * @{ */

  /** @brief Selects next scenario in Dupačová's forward algorithm.
   *
   * Given the current reduced set, selects the scenario that minimizes the
   * Wasserstein distance when added to the set. Uses optimal probability
   * redistribution for distance calculation.
   *
   * @param ind current binary selection vector
   * @param min_cost precomputed minimum distances for efficiency
   * @return tuple of (selected scenario index, index in remaining set)
   */
  std::tuple<int, int> pick_dupacova(const std::vector<bool>& ind, 
    const std::vector<double>& min_cost);

  /** @brief Implements Dupačová's forward selection algorithm.
   *
   * Greedily builds the reduced set by selecting scenarios that minimize
   * the Wasserstein distance at each step.
   *
   * @return Solver::kOK if successful
   */
  int compute_dupacova();

  /** @brief Implements local search algorithms (BestFit/FirstFit).
   *
   * Iteratively improves an initial solution by swapping scenarios until
   * no improvement is found.
   *
   * @return Solver::kOK if successful
   */
  int compute_local_search();

  /** @brief Implements baseline algorithm selecting highest weights.
   *
   * Simply selects the k scenarios with the largest probability weights.
   *
   * @return Solver::kOK if successful
   */
  int compute_baseline();

  /** @brief Initializes starting solution for local search.
   *
   * Creates initial reduced set either from warm start indices or randomly.
   *
   * @return Wasserstein distance of initial solution
   */
  double init_local_search();

  /** @brief Checks if a swap provides sufficient improvement.
   *
   * @param trial_d candidate distance after swap
   * @param curr_d current distance
   * @param dist_dupa reference distance from Dupačová
   * @return true if improvement is sufficient based on rho parameter
   */
  bool improvement_condition(double trial_d, double curr_d, double dist_dupa) const;

  /** @brief Finds best swap using exhaustive search (BestFit).
   *
   * @param curr_d current Wasserstein distance
   * @return tuple of (remove index, add index, new distance)
   */
  std::tuple<Index, Index, double> pick_ij_bestfit(double& curr_d);

  /** @brief Helper for BestFit to evaluate adding a scenario.
   *
   * @param curr_indices current reduced set indices
   * @return pair of (best scenario to add, resulting distance)
   */
  IndexDistancePair bestfit_selection(const std::vector<Index>& curr_indices);

  /** @brief Finds first improving swap (FirstFit).
   *
   * @param curr_d current Wasserstein distance
   * @return tuple of (remove index, add index, new distance)
   */
  std::tuple<Index, Index, double> pick_ij_firstfit(double& curr_d);

  /** @brief Helper for FirstFit to evaluate adding a scenario.
   *
   * @param curr_indices current reduced set indices
   * @param curr_d current distance for early termination
   * @return pair of (first improving scenario, resulting distance)
   */
  IndexDistancePair firstfit_selection(
    const std::vector<Index>& curr_indices, 
    double curr_d);

  /** @brief Swaps a scenario in the reduced set.
   *
   * @param i index of scenario to remove from reduced set
   * @param j index of scenario to add to reduced set
   */
  void swap_indices(Index i, Index j);

  /** @brief Updates binary solution vector from index list.
   *
   * Synchronizes reduced_atoms with ind_red after modifications.
   */
  void update_reduced_atoms();

  /** @brief Validates warm start indices for correctness.
   *
   * @param indices proposed warm start scenario indices
   * @param n total number of scenarios
   * @param m required size of reduced set
   * @throws std::invalid_argument if indices are invalid or duplicated
   */
  void validate_warmstart_indices(const std::vector<Index>& indices, int n, int m);

/** @} ---------------------------------------------------------------------*/

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
