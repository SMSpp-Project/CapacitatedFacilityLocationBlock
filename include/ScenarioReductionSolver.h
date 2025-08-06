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
 * ScenarioReductionSolver interprets the CFL instance as a Discrete 
 * Scenario Reduction Problem and solves it using heuristic methods.
 * 
 * This solver works directly with the physical representation of 
 * the problem and enforces specific requirements on the CFL data:
 * - All capacities must equal 1.0
 * - Demands must represent scenario probabilities (auto-normalized if needed)
 * - Square distance matrix (facilities = customers)
 * 
 * This implementation includes multiple scenario reduction algorithms:
 * - Baseline: Select scenarios with highest probability weights
 * - Dupacova's forward algorithm (default)
 * - Local search with BestFit strategy
 * - Local search with FirstFit strategy
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
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

#include <mutex>
#include <random>
#include <unordered_map>

#include "Solver.h"
#include "CapacitatedFacilityLocationBlock.h" 
#include "ScenarioGenerator.h"  // For ScenarioIndex type

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
 * A purely physical solver for the Discrete Scenario Reduction problem 
 * interpreted as a Capacitated Facility Location (CFL) instance.
 * This solver works directly with the physical data arrays and implements
 * heuristic algorithms.
 * 
 * IMPORTANT REQUIREMENTS for the CapacitatedFacilityLocationBlock:
 * - All facility capacities must be 1.0 (allows sending all probability mass 
 *   to a single facility if needed)
 * - Customer demands represent scenario probabilities and should sum to 1.0
 *   (if they don't, the solver will automatically normalize them)
 * - Number of facilities must equal number of customers (square distance matrix)
 * - Transportation costs represent distances between scenarios in the ell-norm
 * 
 * Implemented heuristic algorithms:
 * - Baseline: Select scenarios with highest probability weights
 * - Dupacova's forward algorithm (default)
 * - Local search with BestFit strategy 
 * - Local search with FirstFit strategy
 * 
 * For exact MILP-based scenario reduction, users should attach a MILPSolver
 * to the CapacitatedFacilityLocationBlock instead.
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
    dblEll = dblLastAlgPar + 1,        ///< Power in ell-Wasserstein distance (default: 2.0)
    dblLastParSRS                      ///< First allowed parameter for derived classes
  };

  /// public enum extending vint_par_type_S for ScenarioReductionSolver
  enum vint_par_type_SRS {
    vintWarmstartIndices = vintLastAlgPar, ///< Custom warm start indices (empty = use Dupacova)
    vintLastParSRS                         ///< First allowed parameter for derived classes
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

  /** Sets the Block that the Solver has to solve and caches physical data.
   *  
   *  The block must be a CapacitatedFacilityLocationBlock with specific requirements:
   *  - All facility capacities must be 1.0
   *  - Customer demands represent scenario probabilities (auto-normalized if needed)
   *  - Number of facilities must equal number of customers 
   *  
   *  @throws std::invalid_argument if block is not CapacitatedFacilityLocationBlock
   *  @throws std::invalid_argument if any capacity != 1.0 and k > 0
   *  @throws std::logic_error if number of facilities != number of customers
   */
  void set_Block(Block* block) override;

  /** Helper method to refresh cached data from the block.
   *  
   *  Validates and caches scenario reduction data:
   *  - Ensures all capacities are 1.0 (when k > 0)
   *  - Normalizes demand probabilities if they don't sum to 1.0
   *  - Caches pointers to transportation costs and normalized weights
   *  
   *  @param k Number of scenarios to select (must be <= number of scenarios)
   *  @throws std::invalid_argument if k < 0 or k > number of scenarios
   *  @throws std::invalid_argument if any capacity != 1.0 and k > 0
   *  @throws std::logic_error if facilities != customers
   */
  void refresh_cached_data(int k);

  /// Solve the scenario reduction problem from a choice of methods 
  virtual int compute(bool changedvars = false) override;

  /// Get the power ell in the objective ell-Wasserstein distance
  [[nodiscard]] float get_ell() const { return ell; }

  /// Returns true if a variable solution is available
  [[nodiscard]] bool has_var_solution() override;

  /// Returns the value of the current solution, if any
  [[nodiscard]] OFValue get_var_value() override;

  /// Writes the current solution in the Block
  void get_var_solution(Configuration* solc = nullptr) override;

  const IntSolution& get_reduced_atoms() const { return reduced_atoms; }

  // Override parameter methods
  void set_par(idx_type par, int value) override;
  void set_par(idx_type par, double value) override;
  void set_par(idx_type par, std::vector<int>&& value) override;
  int get_int_par(idx_type par) const override;
  double get_dbl_par(idx_type par) const override;
  const std::vector<int>& get_vint_par(idx_type par) const override;
  
  // Static default parameter methods
  static int get_dflt_int_par(idx_type par);
  static double get_dflt_dbl_par(idx_type par);

/*--------------------------------------------------------------------------*/

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
  bool use_warmstart = false; // Whether to use warm start for local search
  std::vector<Index> warmstart_indices; // Custom warm start indices (empty = use Dupacova)
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
  
  /// Normalized weights (only allocated if original weights are not normalized)
  std::unique_ptr<DVector> normalized_weights;

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
   * Perform scenario reduction using Baseline algorithm
   * Selects k scenarios with highest probability weights
   * 
   * @return Status code (kOK if successful)
   */
  int compute_baseline();

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

  /**
   * Validate warm start indices
   * 
   * @param indices The indices to validate
   * @param n Total number of scenarios
   * @param m Number of scenarios to select
   * @throws std::invalid_argument if indices are invalid
   */
  void validate_warmstart_indices(const std::vector<Index>& indices, int n, int m);

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