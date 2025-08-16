/*--------------------------------------------------------------------------*/
/*--------------- File ScenarioReductionSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ScenarioReductionSolver class.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

#include "ScenarioReductionSolver.h"

#include <limits>    // std::numeric_limits (not in SMSTypedefs)
#include <numeric>   // std::iota (not in SMSTypedefs)
#include <unordered_map>
#include <unordered_set>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// Register the solver in the factory
SMSpp_insert_in_factory_cpp_1(ScenarioReductionSolver);

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::ScenarioReductionSolver() 
  : nb_atoms(0), nb_reduced(0),
    weights(nullptr), f_transportation_costs(nullptr), f_solution_value(0) 
{}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::~ScenarioReductionSolver() 
{
  // Nothing to clean up
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::set_Block(Block* block)
{
    std::lock_guard<std::recursive_mutex> lock(f_mutex);
    
    // Check if this is the same block we already have
    if (block == f_Block) {
        std::cout << "Block already set, nothing to do" << std::endl;
        return;  // Nothing to do
    }

    // Call base class implementation to set f_Block
    Solver::set_Block(block);

    if (f_Block) {
        // Check if block is correct type
        auto cfl_block = dynamic_cast<CapacitatedFacilityLocationBlock*>(f_Block);
        if (!cfl_block) {
            throw std::invalid_argument("ScenarioReductionSolver only works with CapacitatedFacilityLocationBlock");
        }

        // Lock the block if not owned
        bool owned = f_Block->is_owned_by(f_id);
        if (!owned && !f_Block->lock(f_id)) {
            throw std::runtime_error("ScenarioReductionSolver: unable to lock the Block");
        }

        // Generate abstract variables to write back the solution (y variables)
        // *AND* x variables which are necessary to infer the optimal weights
        f_Block->generate_abstract_variables();

        if (!owned) {
            f_Block->unlock(f_id);
        }

        // Initialize or refresh cached data from the block
        refresh_cached_data(cfl_block->get_NMaxFacilities());

        // Initialize solution structures
        reduced_atoms.resize(nb_atoms, false);
        ind_red.reserve(nb_reduced);
        f_solution_value = 0.0;
    }
}

/*--------------------------------------------------------------------------*/
/*------------------- MAIN COMPUTATION METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::compute(bool changedvars) 
{
  std::lock_guard<std::recursive_mutex> lock(f_mutex);  // Exception-safe RAII lock
  
  // Make sure we have a block to work with
  if (!get_Block()) {
    return kError;
  }
  
  // Clear existing data and solution
  indices_to_choose.clear();
  ind_red.clear();
  std::fill(reduced_atoms.begin(), reduced_atoms.end(), false);
  
  int result;
  // Select the appropriate algorithm
  switch (algorithm) {
    case Algorithm::Baseline:
      result = compute_baseline();
      break;
    case Algorithm::Dupacova:
      result = compute_dupacova();
      break;
    case Algorithm::BestFit:
      [[fallthrough]];  // Both BestFit and FirstFit use compute_local_search()
    case Algorithm::FirstFit:
      result = compute_local_search();
      break;
    default:
      result = kError;
      break;
  }
  
  return result;
  // Lock is automatically released when lock_guard goes out of scope
}

/*--------------------------------------------------------------------------*/
  // Implementation of Dupačová's forward algorithm for scenario reduction
/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::compute_dupacova()
{
  int m = static_cast<int>(nb_reduced);
  int n = static_cast<int>(nb_atoms);
  
  // Initialize the indices to choose from
  indices_to_choose.resize(n);
  std::iota(indices_to_choose.begin(), indices_to_choose.end(), 0);
  
  // For every atom i, save the minimal distance among the current atoms j
  std::vector<double> minimum_d(n, std::numeric_limits<double>::infinity());
  
  for(int k = 0; k < m; ++k) {
    // Find the closest atom to add on a greedy Wasserstein-based criterion
    int j_best, j_tmp;
    std::tie(j_best, j_tmp) = pick_dupacova(reduced_atoms, minimum_d);
    
    // Updates
    for(int i = 0; i < n; i++) {
      minimum_d[i] = std::min(minimum_d[i], (*f_transportation_costs)[i][j_best]);
    }
    reduced_atoms[j_best] = true;
    ind_red.push_back(j_best);

    indices_to_choose.erase(indices_to_choose.begin() + j_tmp);
  }
  
  // Compute final distance
  double dot_product = std::inner_product(minimum_d.begin(), minimum_d.end(), weights->begin(), 0.0);
  dist_dupa = dot_product; // Save the distance from Dupacova for local search
  f_solution_value = dot_product; // Already the ell-th power of Wasserstein distance
  
  return kOK;
}

/*--------------------------------------------------------------------------*/
/*--------------------- BASELINE ALGORITHM IMPLEMENTATION -----------------*/
/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::compute_baseline()
{
  // Create vector of (weight, index) pairs
  std::vector<std::pair<double, Index>> weight_index_pairs;
  weight_index_pairs.reserve(nb_atoms);
  
  for (Index i = 0; i < nb_atoms; ++i) {
    weight_index_pairs.emplace_back((*weights)[i], i);
  }
  
  // Sort by weight in descending order
  std::sort(weight_index_pairs.begin(), weight_index_pairs.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
  
  // Select top k scenarios
  std::fill(reduced_atoms.begin(), reduced_atoms.end(), false);
  ind_red.clear();
  for (Index i = 0; i < nb_reduced && i < nb_atoms; ++i) {
    Index scenario_idx = weight_index_pairs[i].second;
    reduced_atoms[scenario_idx] = true;
    ind_red.push_back(scenario_idx);
  }
  
  // Calculate objective value (ell-th power of Wasserstein distance)
  std::vector<double> min_distances(nb_atoms);
  for (Index i = 0; i < nb_atoms; ++i) {
    min_distances[i] = std::numeric_limits<double>::infinity();
    for (Index j = 0; j < nb_atoms; ++j) {
      if (reduced_atoms[j]) {
        min_distances[i] = std::min(min_distances[i], 
                                   (*f_transportation_costs)[i][j]);
      }
    }
  }
  
  double total_distance = std::inner_product(min_distances.begin(), 
                                           min_distances.end(), 
                                           weights->begin(), 0.0);
  f_solution_value = total_distance; 
  
  return kOK;
}

/*--------------------------------------------------------------------------*/
/*-------------------- LOCAL SEARCH IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::compute_local_search()
{
  // Initialize local search with random indices
  double curr_d = init_local_search();
  
  bool improvement = true;
  while (improvement) {
    int i, j;
    double trial_d;
    
    // Choose the strategy for picking the next pair
    if (algorithm == Algorithm::BestFit) {
      std::tie(i, j, trial_d) = pick_ij_bestfit(curr_d);
    } else {
      std::tie(i, j, trial_d) = pick_ij_firstfit(curr_d);
    }
    
    // Check if we found a valid improvement
    if (i >= 0 && j >= 0 && improvement_condition(trial_d, curr_d, dist_dupa)) {
      curr_d = trial_d;
      swap_indices(i, j);
    } else {
      improvement = false;
    }
  }
  
  // Update reduced_atoms vector from ind_red
  update_reduced_atoms();
  
  // Calculate final ell-th power of Wasserstein distance
  f_solution_value = curr_d;
  
  return kOK;
}

/*--------------------------------------------------------------------------*/

double ScenarioReductionSolver::init_local_search()
{
  int n = static_cast<int>(nb_atoms);
  int m = static_cast<int>(nb_reduced);
  
  // Initialize the indices to choose from
  indices_to_choose.resize(n);
  std::iota(indices_to_choose.begin(), indices_to_choose.end(), 0);
  
  ind_red.clear();
  
  if (use_warmstart) {
    if (!warmstart_indices.empty()) {
      // Use custom warm start indices
      validate_warmstart_indices(warmstart_indices, n, m);
      ind_red = warmstart_indices;
      // Sort to maintain consistency
      std::sort(ind_red.begin(), ind_red.end());
    } else {
      // Use Dupačová as default warm start
      compute_dupacova(); // This will set ind_red and reduced_atoms
    }
    
    // Update indices_to_choose to exclude selected indices
    indices_to_choose.clear();
    for (int i = 0; i < n; ++i) {
      if (std::find(ind_red.begin(), ind_red.end(), i) == ind_red.end()) {
        indices_to_choose.push_back(i);
      }
    }
  } else {
    // Random initialization
    ind_red.clear();
    std::vector<Index> shuffled_indices(n);
    std::iota(shuffled_indices.begin(), shuffled_indices.end(), 0);
    std::shuffle(shuffled_indices.begin(), shuffled_indices.end(), rng);
    
    // Take the first m elements as the initial reduced set
    ind_red.assign(shuffled_indices.begin(), shuffled_indices.begin() + m);
    std::sort(ind_red.begin(), ind_red.end());
    
    // Update indices_to_choose
    indices_to_choose.clear();
    for (int i = 0; i < n; ++i) {
      if (std::find(ind_red.begin(), ind_red.end(), i) == ind_red.end()) {
        indices_to_choose.push_back(i);
      }
    }
  }
  
  // Update reduced_atoms
  update_reduced_atoms();
  
  // Calculate initial Wasserstein distance
  std::vector<double> min_distances(n);
  for (int i = 0; i < n; ++i) {
    min_distances[i] = std::numeric_limits<double>::infinity();
    for (auto j : ind_red) {
      min_distances[i] = std::min(min_distances[i], (*f_transportation_costs)[i][j]);
    }
  }
  
  return std::inner_product(min_distances.begin(), min_distances.end(), weights->begin(), 0.0);
}

/*--------------------------------------------------------------------------*/

bool ScenarioReductionSolver::improvement_condition(
  double trial_d, double curr_d, double dist_dupa) const
{
  if (rho <= 0.0) {
    return trial_d < curr_d;
  } else {
    return trial_d < curr_d - rho * dist_dupa;
  }
}

/*--------------------------------------------------------------------------*/
/*------------------------ BESTFIT IMPLEMENTATION -------------------------*/
/*--------------------------------------------------------------------------*/

std::tuple<ScenarioReductionSolver::Index, ScenarioReductionSolver::Index, double>
ScenarioReductionSolver::pick_ij_bestfit(double& curr_d)
{
  // Holder for distances
  std::vector<double> distances(ind_red.size(), std::numeric_limits<double>::infinity());
  std::unordered_map<Index, Index> j_map;
  
  // Try removing each atom in the reduced set
  for (size_t i = 0; i < ind_red.size(); ++i) {
    // Temporarily remove this atom
    Index atom_to_remove = ind_red[i];
    ind_red.erase(ind_red.begin() + i);
    
    // Find the best atom to add
    auto [best_j, dist] = bestfit_selection(ind_red);
    
    // Remember this combination
    j_map[atom_to_remove] = best_j;
    distances[i] = dist;
    
    // Put the atom back
    ind_red.insert(ind_red.begin() + i, atom_to_remove);
  }
  
  // Find the best combination
  auto min_it = std::min_element(distances.begin(), distances.end());
  if (min_it == distances.end()) {
    return std::make_tuple(-1, -1, std::numeric_limits<double>::infinity());
  }
  
  size_t best_idx = std::distance(distances.begin(), min_it);
  Index best_i = ind_red[best_idx];
  Index best_j = j_map[best_i];
  
  return std::make_tuple(best_i, best_j, distances[best_idx]);
}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::IndexDistancePair
ScenarioReductionSolver::bestfit_selection(const std::vector<Index>& curr_indices)
{
  int n = static_cast<int>(nb_atoms);
  
  // Calculate minimum distances to current reduced set
  std::vector<double> min_on_ind_red(n, std::numeric_limits<double>::infinity());
  for (int i = 0; i < n; ++i) {
    for (auto j : curr_indices) {
      min_on_ind_red[i] = std::min(min_on_ind_red[i], (*f_transportation_costs)[i][j]);
    }
  }
  
  // Find the best atom to add from indices_to_choose
  Index best_j = -1;
  double best_dist = std::numeric_limits<double>::infinity();
  
  for (auto j : indices_to_choose) {
    // Calculate combined minimums
    std::vector<double> combined_min(n);
    for (int i = 0; i < n; ++i) {
      combined_min[i] = std::min(min_on_ind_red[i], (*f_transportation_costs)[i][j]);
    }
    
    // Calculate objective value
    double obj_val = std::inner_product(combined_min.begin(), combined_min.end(), weights->begin(), 0.0);
    
    // Update best if improvement found
    if (obj_val < best_dist) {
      best_dist = obj_val;
      best_j = j;
    }
  }
  
  return {best_j, best_dist};
}

/*--------------------------------------------------------------------------*/
/*------------------------ FIRSTFIT IMPLEMENTATION -------------------------*/
/*--------------------------------------------------------------------------*/

std::tuple<ScenarioReductionSolver::Index, ScenarioReductionSolver::Index, double>
ScenarioReductionSolver::pick_ij_firstfit(double& curr_d)
{
  // If shuffling is enabled, shuffle the reduced set
  if (shuffle) {
    std::shuffle(ind_red.begin(), ind_red.end(), rng);
  }
  
  // Try each atom in the reduced set
  for (auto i : ind_red) {
    // Create a temporary set without i
    std::vector<Index> temp_indices;
    temp_indices.reserve(ind_red.size() - 1);
    for (auto idx : ind_red) {
      if (idx != i) {
        temp_indices.push_back(idx);
      }
    }
    
    // Try to find a suitable replacement
    auto [j, dist] = firstfit_selection(temp_indices, curr_d);
    
    // If we found an improvement, return it
    if (j != static_cast<Index>(-1)) {
      return {i, j, dist};
    }
  }
  
  // If we need to restore the original order of ind_red
  if (shuffle) {
    std::sort(ind_red.begin(), ind_red.end());
  }
  
  // No improvement found
  return {-1, -1, std::numeric_limits<double>::infinity()};
}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::IndexDistancePair
ScenarioReductionSolver::firstfit_selection(
  const std::vector<Index>& curr_indices,
  double curr_d)
{
  int n = static_cast<int>(nb_atoms);
  
  // Calculate minimum distances to current reduced set
  std::vector<double> min_on_ind_red(n, std::numeric_limits<double>::infinity());
  for (int i = 0; i < n; ++i) {
    for (auto j : curr_indices) {
      min_on_ind_red[i] = std::min(min_on_ind_red[i], (*f_transportation_costs)[i][j]);
    }
  }
  
  // Try each atom in indices_to_choose
  for (auto j : indices_to_choose) {
    // Calculate combined minimums
    std::vector<double> combined_min(n);
    for (int i = 0; i < n; ++i) {
      combined_min[i] = std::min(min_on_ind_red[i], (*f_transportation_costs)[i][j]);
    }
    
    // Calculate objective value
    double trial_d = std::inner_product(combined_min.begin(), combined_min.end(), weights->begin(), 0.0);
    
    // Return first improvement found
    if (improvement_condition(trial_d, curr_d, dist_dupa)) {
      return {j, trial_d};
    }
  }
  
  // No improvement found
  return {-1, std::numeric_limits<double>::infinity()};
}

/*--------------------------------------------------------------------------*/
/*------------------------ DUPACOVA IMPLEMENTATION -------------------------*/
/*--------------------------------------------------------------------------*/

std::tuple<int, int> ScenarioReductionSolver::pick_dupacova(
  const std::vector<bool>& red_ind, 
  const std::vector<double>& min_cost) 
{
  // For every atom j in indices_to_choose, compute Wasserstein distance with closed formula
  std::vector<double> inner_min(nb_atoms); // inner_min in closed formula
  std::vector<float> distances(indices_to_choose.size());

  for (size_t idx = 0; idx < indices_to_choose.size(); idx++) {
    auto j = indices_to_choose[idx];
    // compute every component 0\leq i \leq n-1 of inner_min by recursive formula
    for (int i = 0; i < nb_atoms; i++) {
      inner_min[i] = std::min(min_cost[i], (*f_transportation_costs)[i][j]);
    }
    distances[idx] = std::inner_product(inner_min.begin(), inner_min.end(), weights->begin(), 0.0); 
  }

  // compute argmin_j distances[j]
  auto min_it = std::min_element(distances.begin(), distances.end());
  int j_tmp = std::distance(distances.begin(), min_it); // index in indices_to_choose

  return std::make_pair(indices_to_choose[j_tmp], j_tmp);
}

/*--------------------------------------------------------------------------*/
/*-------------------------- HELPER METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::swap_indices(Index i, Index j)
{
  // Find and remove i from ind_red
  auto it_i = std::find(ind_red.begin(), ind_red.end(), i);
  if (it_i != ind_red.end()) {
    ind_red.erase(it_i);
  }
  
  // Find and remove j from indices_to_choose
  auto it_j = std::find(indices_to_choose.begin(), indices_to_choose.end(), j);
  if (it_j != indices_to_choose.end()) {
    indices_to_choose.erase(it_j);
  }
  
  // Add j to ind_red and i to indices_to_choose
  ind_red.push_back(j);
  indices_to_choose.push_back(i);
  
  // Sort the vectors to maintain order
  std::sort(ind_red.begin(), ind_red.end());
  std::sort(indices_to_choose.begin(), indices_to_choose.end());
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::update_reduced_atoms()
{
  // Reset all atoms to false
  std::fill(reduced_atoms.begin(), reduced_atoms.end(), false);
  
  // Set atoms in ind_red to true
  for (auto i : ind_red) {
    if (i >= 0 && i < static_cast<Index>(reduced_atoms.size())) {
      reduced_atoms[i] = true;
    }
  }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::get_var_solution(Configuration* solc)
{
  std::lock_guard<std::recursive_mutex> lock(f_mutex);
  
  if (!f_Block) {
    throw std::logic_error("ScenarioReductionSolver::get_var_solution: no Block set");
  }
  
  if (!has_var_solution()) {
    throw std::logic_error("ScenarioReductionSolver::get_var_solution: no solution available");
  }
  
  auto cfl_block = dynamic_cast<CapacitatedFacilityLocationBlock*>(f_Block);
  if (!cfl_block) {
    throw std::logic_error("ScenarioReductionSolver::get_var_solution: Block is not CapacitatedFacilityLocationBlock");
  }
  
  // Check that abstract representation exists
  if (!cfl_block->get_y(0)) {
    throw std::logic_error("ScenarioReductionSolver::get_var_solution: variables not generated in Block");
  }
  
  // Write solution to Block's y variables
  for (Index i = 0; i < nb_atoms; ++i) {
    ColVariable* y_var = cfl_block->get_y(i);
    if (y_var) {
      // Set value: 1.0 if facility/scenario selected, 0.0 otherwise
      y_var->set_value(reduced_atoms[i] ? 1.0 : 0.0);
    }
  }
  
  // Write solution to Block's x variables (assignments)
  // For each customer (scenario), find the closest open facility (representative)
  // and set x[facility][customer] = 1.0 for that assignment
  if (f_transportation_costs) {
    for (Index customer = 0; customer < nb_atoms; ++customer) {
      // Find the closest open facility for this customer
      Index best_facility = 0;
      double min_cost = std::numeric_limits<double>::infinity();
      
      for (Index facility = 0; facility < nb_atoms; ++facility) {
        // Only consider open facilities (selected representatives)
        if (reduced_atoms[facility]) {
          double cost = (*f_transportation_costs)[customer][facility];
          if (cost < min_cost) {
            min_cost = cost;
            best_facility = facility;
          }
        }
      }
      
      // Set all x variables for this customer
      for (Index facility = 0; facility < nb_atoms; ++facility) {
        ColVariable* x_var = cfl_block->get_x(facility, customer);
        if (x_var) {
          // Set 1.0 for the best assignment, 0.0 for all others
          x_var->set_value((facility == best_facility && reduced_atoms[facility]) ? 1.0 : 0.0);
        }
      }
    }
  }
}

/*--------------------------------------------------------------------------*/

bool ScenarioReductionSolver::has_var_solution() 
{
  // Solution exists if compute() was successful and we have selected scenarios
  return (f_solution_value >= 0) && !reduced_atoms.empty() && 
         std::any_of(reduced_atoms.begin(), reduced_atoms.end(), 
                     [](bool val) { return val; });
}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::OFValue ScenarioReductionSolver::get_var_value() 
{
  return f_solution_value;
}

/*--------------------------------------------------------------------------*/
/*---------------------- PARAMETER METHODS ---------------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::set_par(idx_type par, int value) {
  switch(par) {
    case intAlgorithm:
      if (value < 0 || value > 3) {
        throw std::invalid_argument("Invalid algorithm value");
      }
      algorithm = static_cast<Algorithm>(value);
      break;
    case intShuffle:
      shuffle = (value != 0);
      break;
    case intRandomSeed:
      rng.seed(static_cast<unsigned int>(value));
      break;
    case intUseWarmstart:
      use_warmstart = (value != 0);
      break;
    default:
      Solver::set_par(par, value);
  }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::set_par(idx_type par, double value) {
  switch(par) {
    case dblRho:
      if (value < 0.0) {
        throw std::invalid_argument("rho must be non-negative");
      }
      rho = value;
      break;
    default:
      Solver::set_par(par, value);
  }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::set_par(idx_type par, std::vector<int>&& value) {
  switch(par) {
    case vintWarmstartIndices:
      warmstart_indices.clear();
      warmstart_indices.reserve(value.size());
      for (int idx : value) {
        warmstart_indices.push_back(static_cast<Index>(idx));
      }
      break;
    default:
      Solver::set_par(par, std::move(value));
  }
}

/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::get_int_par(idx_type par) const {
  switch(par) {
    case intAlgorithm:
      return static_cast<int>(algorithm);
    case intShuffle:
      return shuffle ? 1 : 0;
    case intRandomSeed:
      // Note: We can't retrieve the seed from mt19937, so return a default
      return 0;
    case intUseWarmstart:
      return use_warmstart ? 1 : 0;
    default:
      return Solver::get_int_par(par);
  }
}

/*--------------------------------------------------------------------------*/

double ScenarioReductionSolver::get_dbl_par(idx_type par) const {
  switch(par) {
    case dblRho:
      return rho;
    default:
      return Solver::get_dbl_par(par);
  }
}

/*--------------------------------------------------------------------------*/

const std::vector<int>& ScenarioReductionSolver::get_vint_par(idx_type par) const {
  if (par == vintWarmstartIndices) {
    // We need to return a const ref to vector<int>, but we have vector<Index>
    // Create a static thread_local to hold the conversion
    static thread_local std::vector<int> temp_indices;
    temp_indices.clear();
    temp_indices.reserve(warmstart_indices.size());
    for (Index idx : warmstart_indices) {
      temp_indices.push_back(static_cast<int>(idx));
    }
    return temp_indices;
  }
  return Solver::get_vint_par(par);
}

/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::get_dflt_int_par(idx_type par) const {
  switch(par) {
    case intAlgorithm: return 1;    // Dupacova
    case intShuffle: return 0;       // No shuffling
    case intRandomSeed: return 0;    // Default seed
    case intUseWarmstart: return 0;  // No warm start by default
    default: return 0;  // Base class default
  }
}

/*--------------------------------------------------------------------------*/

double ScenarioReductionSolver::get_dflt_dbl_par(idx_type par) const {
  switch(par) {
    case dblRho: return 0.0;
    default: return 0.0;  // Base class default
  }
}

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::refresh_cached_data(int k)
{
  // Check if block exists
  if (!f_Block) {
    throw std::logic_error("ScenarioReductionSolver::refresh_cached_data: no Block set");
  }
  
  // Cast to the specialized type
  auto cfl_block = dynamic_cast<CapacitatedFacilityLocationBlock*>(f_Block);
  if (!cfl_block) {
    throw std::logic_error("ScenarioReductionSolver::refresh_cached_data: Block is not CapacitatedFacilityLocationBlock");
  }
  
  // Validate k parameter
  if (k < 0) {
    throw std::invalid_argument("ScenarioReductionSolver::refresh_cached_data: k must be non-negative");
  }
  
  // Now use the specialized methods with the properly typed pointer
  nb_atoms = static_cast<Index>(cfl_block->get_NCustomers());
  nb_reduced = static_cast<Index>(cfl_block->get_NFacilities());
  
  // Validate that we have a square distance matrix (customers == facilities for scenario reduction)
  if (nb_atoms != nb_reduced) {
    throw std::logic_error("ScenarioReductionSolver::refresh_cached_data: for scenario reduction, number of customers must equal number of facilities");
  }
  
  // Validate k against problem size
  if (k > nb_atoms) {
    throw std::invalid_argument("ScenarioReductionSolver::refresh_cached_data: k cannot exceed number of scenarios");
  }
  
  // Get pointers to data
  const DVector& original_demands = cfl_block->get_Demands();
  const DVector& capacities = cfl_block->get_Capacities();
  f_transportation_costs = &cfl_block->get_Transportation_Costs();
  
  if (original_demands.empty()) {
    throw std::runtime_error("ScenarioReductionSolver::refresh_cached_data: Block has no demand data");
  }
  
  if (capacities.empty()) {
    throw std::runtime_error("ScenarioReductionSolver::refresh_cached_data: Block has no capacity data");
  }
  
  if (!f_transportation_costs) {
    throw std::runtime_error("ScenarioReductionSolver::refresh_cached_data: Block has no transportation cost data");
  }
  
  // Define tolerance for floating point comparisons
  const double tolerance = 1e-6;
  
  // For scenario reduction, capacities should all be 1.0 to allow sending all mass to a single facility
  if (k > 0) {
    for (size_t i = 0; i < capacities.size(); ++i) {
      if (std::abs(capacities[i] - 1.0) > tolerance) {
        throw std::invalid_argument("ScenarioReductionSolver::refresh_cached_data: for scenario reduction, all capacities must be 1.0. Capacity " + 
                                   std::to_string(i) + " is " + std::to_string(capacities[i]));
      }
    }
  }
  
  // Check if demands are already normalized (sum to 1.0 within tolerance)
  double sum = std::accumulate(original_demands.begin(), original_demands.end(), 0.0);
  
  if (k > 0 && std::abs(sum - 1.0) > tolerance) {
    // Normalize the weights and store in normalized_weights
    normalized_weights = std::make_unique<DVector>(original_demands.size());
    for (size_t i = 0; i < original_demands.size(); ++i) {
      (*normalized_weights)[i] = original_demands[i] / sum;
    }
    weights = normalized_weights.get();
  } else {
    // Weights are already normalized or k=0, use original demands
    normalized_weights.reset();  // Free any existing normalized weights
    weights = &original_demands;
  }
  
  // Resize data structures and set parameters
  reduced_atoms.resize(nb_atoms, false);
  ind_red.clear();
  ind_red.reserve(k);
  nb_reduced = k;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::validate_warmstart_indices(
    const std::vector<Index>& indices, int n, int m)
{
  // Check size
  if (indices.size() != static_cast<size_t>(m)) {
    throw std::invalid_argument(
      "ScenarioReductionSolver::validate_warmstart_indices: warmstart_indices size (" + 
      std::to_string(indices.size()) + ") must equal nb_reduced (" + 
      std::to_string(m) + ")");
  }
  
  // Check bounds and uniqueness
  std::unordered_set<Index> seen;
  for (Index idx : indices) {
    if (idx < 0 || idx >= n) {
      throw std::invalid_argument(
        "ScenarioReductionSolver::validate_warmstart_indices: warmstart index " + 
        std::to_string(idx) + " out of bounds [0, " + 
        std::to_string(n-1) + "]");
    }
    if (!seen.insert(idx).second) {
      throw std::invalid_argument(
        "ScenarioReductionSolver::validate_warmstart_indices: duplicate warmstart index " + 
        std::to_string(idx));
    }
  }
}

/*--------------------------------------------------------------------------*/
/*------------------ End File ScenarioReductionSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
