/*--------------------------------------------------------------------------*/
/*--------------- File ScenarioReductionSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ScenarioReductionSolver class.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 */

#include "ScenarioReductionSolver.h"

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
  : f_CFLBlock(nullptr), f_n_facilities(0), f_n_customers(0),
    f_capacities(nullptr), f_fixed_costs(nullptr), f_demands(nullptr),
    f_transportation_costs(nullptr), f_unsplittable(false),
    f_solution_value(0) 
{
  // Nothing else to do in constructor
}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::~ScenarioReductionSolver() 
{
  // No need to delete the data pointers - they're just references to the block's data
}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::set_Block(Block* block) 
{
  // First call the base class implementation to set f_Block
  Solver::set_Block(block);
  
  // Try to cast the block to CapacitatedFacilityLocationBlock
  f_CFLBlock = dynamic_cast<CapacitatedFacilityLocationBlock*>(block);
  if (!f_CFLBlock) {
    throw std::invalid_argument("ScenarioReductionSolver only works with CapacitatedFacilityLocationBlock");
  }
  
  // Cache the problem dimensions and references to the data
  refresh_cached_data();
  
  // Initialize solution structures
  f_facility_solution.resize(f_n_facilities, false);
  f_transportation_solution.resize(f_n_facilities * f_n_customers, 0.0);
  f_solution_value = 0.0;
}

/*--------------------------------------------------------------------------*/

int ScenarioReductionSolver::compute(bool changedvars) 
{
  // Make sure we have a block to work with
  if (!f_CFLBlock) {
    return kError;
  }
  
  // Process any pending modifications
  process_pending_modifications();
  
  // This is where you'll implement your scenario reduction algorithm
  // For now, we'll just create a dummy implementation
  
  // TODO: Replace with your actual algorithm implementation
  // ---------------------------------------------------------
  // Example of how you might use the physical data:
  
  // Clear existing solution
  std::fill(f_facility_solution.begin(), f_facility_solution.end(), false);
  std::fill(f_transportation_solution.begin(), f_transportation_solution.end(), 0.0);
  
  // Simple greedy algorithm (just for illustration)
  // 1. Sort facilities by fixed cost
  std::vector<std::pair<Cost, Index>> sorted_facilities;
  for (Index i = 0; i < f_n_facilities; ++i) {
    sorted_facilities.push_back({(*f_fixed_costs)[i], i});
  }
  std::sort(sorted_facilities.begin(), sorted_facilities.end());
  
  // 2. Open facilities in order of increasing fixed cost until all demand is satisfied
  double total_demand = 0.0;
  for (const auto& demand : *f_demands) {
    total_demand += demand;
  }
  
  double capacity_available = 0.0;
  for (const auto& [cost, i] : sorted_facilities) {
    if (capacity_available >= total_demand) {
      break;
    }
    f_facility_solution[i] = true;
    capacity_available += (*f_capacities)[i];
  }
  
  // 3. Assign customers to facilities (simple greedy allocation)
  for (Index j = 0; j < f_n_customers; ++j) {
    Demand remaining_demand = (*f_demands)[j];
    for (Index i = 0; i < f_n_facilities; ++i) {
      if (!f_facility_solution[i]) continue;
      
      // Find available capacity for this facility
      Demand facility_used = 0.0;
      for (Index j2 = 0; j2 < f_n_customers; ++j2) {
        facility_used += f_transportation_solution[i * f_n_customers + j2] * (*f_demands)[j2];
      }
      Demand available = (*f_capacities)[i] - facility_used;
      
      if (available > 0) {
        Demand assigned = std::min(remaining_demand, available);
        f_transportation_solution[i * f_n_customers + j] = assigned / (*f_demands)[j];
        remaining_demand -= assigned;
        if (remaining_demand <= 0) break;
      }
    }
  }
  
  // 4. Calculate objective value
  f_solution_value = 0.0;
  for (Index i = 0; i < f_n_facilities; ++i) {
    if (f_facility_solution[i]) {
      f_solution_value += (*f_fixed_costs)[i];
      for (Index j = 0; j < f_n_customers; ++j) {
        f_solution_value += f_transportation_solution[i * f_n_customers + j] * (*f_transportation_costs)[i][j];
      }
    }
  }
  // ---------------------------------------------------------
  
  return kOK;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::get_var_solution(Configuration* solc) 
{
  // Make sure we have a block to work with
  if (!f_CFLBlock) {
    throw std::logic_error("No CapacitatedFacilityLocationBlock set");
  }
  
  // For now, do nothing - we won't try to write to the block's variables
  // since we're ignoring the abstract representation
  
  // Later replace this with:
  // if (f_CFLBlock has abstract representation) {
  //   f_CFLBlock->set_facility_solution(f_facility_solution.begin());
  //   f_CFLBlock->set_transportation_solution(f_transportation_solution.begin());
  // }
}

// For get_var_value():
OFValue ScenarioReductionSolver::get_var_value() 
{
  // Simply return your internally calculated objective value
  return f_solution_value;
}

/*--------------------------------------------------------------------------*/

ScenarioReductionSolver::OFValue ScenarioReductionSolver::get_var_value() 
{
  return f_solution_value;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::add_Modification(sp_Mod& mod) 
{
  // First call the base implementation to store the modification in the queue
  Solver::add_Modification(mod);
  
  // For performance, we might want to process modifications immediately
  // but for simplicity, we'll process them all together in compute()
}

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::refresh_cached_data() 
{
  if (!f_CFLBlock) return;
  
  // Cache the problem dimensions
  f_n_facilities = f_CFLBlock->get_NFacilities();
  f_n_customers = f_CFLBlock->get_NCustomers();
  
  // Cache references to the problem data
  // Note: We're storing pointers to the data in the block, not copying it
  f_capacities = &f_CFLBlock->get_Capacities();
  f_fixed_costs = &f_CFLBlock->get_Fixed_Costs();
  f_demands = &f_CFLBlock->get_Demands();
  f_transportation_costs = &f_CFLBlock->get_Transportation_Costs();
  f_unsplittable = f_CFLBlock->get_UnSplittable();
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionSolver::process_pending_modifications() 
{
  // Process all pending modifications
  sp_Mod mod;
  while ((mod = pop()) != nullptr) {
    // Handle different types of modifications
    
    // Check for changes in facility costs
    if (auto cflMod = std::dynamic_pointer_cast<CapacitatedFacilityLocationBlockRngdMod>(mod)) {
      if (cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgFCost ||
          cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgTCost ||
          cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgCap ||
          cflMod->type() == CapacitatedFacilityLocationBlockMod::eChgDem) {
        // Data has changed, refresh our cached references
        refresh_cached_data();
      }
    }
    // Check for subset-based modifications
    else if (auto cflSubsetMod = std::dynamic_pointer_cast<CapacitatedFacilityLocationBlockSbstMod>(mod)) {
      // Similar handling as above
      refresh_cached_data();
    }
    // Check for nuclear option (reload)
    else if (std::dynamic_pointer_cast<NBModification>(mod)) {
      // Complete reload of the problem
      refresh_cached_data();
      
      // Reinitialize solution structures
      f_facility_solution.resize(f_n_facilities, false);
      f_transportation_solution.resize(f_n_facilities * f_n_customers, 0.0);
      f_solution_value = 0.0;
    }
  }
}

/*--------------------------------------------------------------------------*/
/*------------------ End File ScenarioReductionSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/