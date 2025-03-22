/*--------------------------------------------------------------------------*/
/*------------------ File ScenarioReductionSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ScenarioReductionSolver class, which solves
 * the scenario reduction problem formulated as a Capacitated Facility
 * Location Problem.
 * 
 * For now, this solver works directly with the physical representation of 
 * the problem.
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

#include "Solver.h"
#include "CapacitatedFacilityLocationBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioReductionSolver --------------------*/
/*--------------------------------------------------------------------------*/
/**
 * A solver for Capacitated Facility Location problems that works directly
 * with the physical representation of the problem, bypassing the abstract
 * representation. This solver implements scenario reduction techniques.
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

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

  /// Default constructor
  ScenarioReductionSolver();

  /// Destructor
  ~ScenarioReductionSolver() override;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

  /// Sets the Block that the Solver has to solve and caches physical data
  void set_Block(Block* block) override;

  /// Optimizes the problem using scenario reduction techniques
  int compute(bool changedvars = false) override;

  /// Returns the value of the current solution, if any
  OFValue get_var_value() override;

  /// Writes the current solution in the Block
  void get_var_solution(Configuration* solc = nullptr) override;

  /// Reacts to changes in the Block's data
  void add_Modification(sp_Mod& mod) override;

  /// TODO: Add proper :Solution when the abstract representation is done 
  const IntSolution& get_facility_solution() const { return f_facility_solution; }

  /// TODO: Add proper :Solution when the abstract representation is done 
  const CntSolution& get_transportation_solution() const { return f_transportation_solution; }

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

  /// The CapacitatedFacilityLocationBlock we're solving
  CapacitatedFacilityLocationBlock* f_CFLBlock;

  /// Problem dimensions
  Index f_n_facilities;
  Index f_n_customers;

  /// Problem data - these are references to the data in the block
  const DVector* f_capacities;
  const CVector* f_fixed_costs;
  const DVector* f_demands;
  const CMatrix* f_transportation_costs;
  bool f_unsplittable;

  /// Solution data
  IntSolution f_facility_solution;   // Binary values indicating which facilities are open
  CntSolution f_transportation_solution; // Continuous values indicating allocation
  double f_solution_value;           // Objective function value of current solution

  /// Helper method to refresh cached data from the block
  void refresh_cached_data();

  /// Helper method to process modifications
  void process_pending_modifications();

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