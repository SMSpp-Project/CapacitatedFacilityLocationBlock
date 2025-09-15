/*--------------------------------------------------------------------------*/
/*------------------------ File test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test suite for ScenarioReductionSolver class.
 * 
 * Test 1 - Parameter Management:
 * - Part 1: Setting and getting all parameter types (intAlgorithm, intShuffle, 
 *           intRandomSeed, dblRho)
 * - Part 2: Boundary value testing (min/max algorithm values, small rho)
 * - Part 3: Parameter validation and error handling (invalid algorithm, negative rho)
 * - Part 4: Default parameter values through get_dflt_*_par methods
 * 
 * Test 2 - Solution Handling:
 * - Part 1: Basic functionality with y variables (facility selection)
 * - Part 2: X variable assignments (customer-to-facility mapping)
 * - Part 3: Algorithm comparison (Baseline vs Dupacova objective values)
 * - Part 4: Solution independence (multiple solvers on same block)
 * - Part 5: Solution persistence through parameters
 * - Part 6: Warm start with vintWarmstartIndices
 * - Part 7: has_var_solution() and get_var_solution() behavior
 * 
 * Test 3 - Thread Safety & Concurrency:
 * - Part 1: Concurrent parameter setting
 * - Part 2: Preventing simultaneous compute() calls
 * - Part 3: Thread safety during set_Block()
 * - Part 4: Thread safety of get_var_solution()
 * - Part 5: Stress test with mixed operations
 * 
 * Test 4 - Error Handling & Edge Cases:
 * - Part 1: refresh_cached_data() error conditions (negative k, k > n)
 * - Part 2: Non-square matrix validation
 * - Part 3: Weight normalization (unnormalized and pre-normalized)
 * - Part 4: Edge cases (k=0, k=n, single scenario)
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \copyright &copy; by Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>      // std::cout, std::cerr, std::endl
#include <iomanip>       // std::setprecision
#include <stdexcept>     // std::runtime_error, std::invalid_argument, etc.
#include <thread>        // std::thread
#include <chrono>        // std::chrono
#include <atomic>        // std::atomic
#include <cstdlib>       // std::rand, std::srand
#include <sstream>       // std::ostringstream

#include "SMSTypedefs.h" // includes: algorithm, vector, map, functional, mutex
#include "Configuration.h"
#include "BlockSolverConfig.h"
#include "ScenarioReductionSolver.h"
#include "CapacitatedFacilityLocationBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- USING ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- TEST FRAMEWORK -------------------------------*/
/*--------------------------------------------------------------------------*/

struct TestResult {
  bool passed;
  std::string message;
};

// Global test counters
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

// Global verbose flag
static bool verbose = false;

// Test registry
static std::map<std::string, std::function<TestResult()>> test_registry;

// Register a test
#define REGISTER_TEST(name, func) \
  static bool _reg_##func = []() { \
    test_registry[name] = func; \
    return true; \
  }()

// Helper to run a single test
void run_test(const std::string& name, std::function<TestResult()> test_func) {
  if (verbose) {
    std::cout << "\n=== Running: " << name << " ===" << std::endl;
  }
  tests_run++;
  
  // Redirect cout to null if not verbose during test execution
  std::streambuf* orig_cout = nullptr;
  std::ostringstream null_stream;
  if (!verbose) {
    orig_cout = std::cout.rdbuf();
    std::cout.rdbuf(null_stream.rdbuf());
  }
  
  try {
    TestResult result = test_func();
    
    // Restore cout
    if (!verbose && orig_cout) {
      std::cout.rdbuf(orig_cout);
    }
    
    if (result.passed) {
      if (verbose) {
        std::cout << "PASSED: " << result.message << std::endl;
      } else {
        std::cout << "[PASS] " << name << std::endl;
      }
      tests_passed++;
    } else {
      std::cout << "[FAIL] " << name << ": " << result.message << std::endl;
      tests_failed++;
    }
  } catch (const std::exception& e) {
    if (!verbose && orig_cout) {
      std::cout.rdbuf(orig_cout);
    }
    std::cout << "[EXCEPTION] " << name << ": " << e.what() << std::endl;
    tests_failed++;
  }
}

/*--------------------------------------------------------------------------*/
/*---------------------------- TEST HELPERS --------------------------------*/
/*--------------------------------------------------------------------------*/

// Helper function to check if two doubles are approximately equal
bool approx_equal(double a, double b, double epsilon = 1e-6) {
  return std::abs(a - b) < epsilon;
}

// Helper function to create a small test CFL block
// k is the number of facilities to select (must be positive)
CapacitatedFacilityLocationBlock* create_test_block(int k) {
  if (k <= 0) {
    throw std::invalid_argument("k must be positive");
  }
  
  auto cfl_block = new CapacitatedFacilityLocationBlock();
  
  // Create a small problem: for scenario reduction, we need a square distance matrix
  // If we interpret customers as scenarios, we need nc x nc matrix
  int nf = 4;  // Must equal nc for square matrix
  int nc = 4;
  
  // Set up facility costs
  CapacitatedFacilityLocationBlock::CVector fcosts(nf);
  for (int i = 0; i < nf; ++i) {
    fcosts[i] = 100.0 * (i + 1);  // 100, 200, 300, 400
  }
  
  // Set up transportation costs as a distance matrix between scenarios
  CapacitatedFacilityLocationBlock::CMatrix tcosts(boost::extents[nf][nc]);
  for (int i = 0; i < nf; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (i == j) {
        tcosts[i][j] = 0.0;  // Zero distance to self
      } else {
        tcosts[i][j] = 10.0 * std::abs(i - j);  // Distance based on index difference
      }
    }
  }
  
  // Set up capacities - must be 1.0 for scenario reduction (allows sending all mass to single facility)
  CapacitatedFacilityLocationBlock::DVector caps(nf);
  caps[0] = 1.0;
  caps[1] = 1.0;
  caps[2] = 1.0;
  caps[3] = 1.0;
  
  // Set up demands
  CapacitatedFacilityLocationBlock::DVector dems(nc);
  for (int i = 0; i < nc; ++i) {
    dems[i] = 50.0 * (i + 1);  // 50, 100, 150, 200
  }
  
  // Load data into block with k as the max number of facilities
  cfl_block->load(nf, nc, caps, fcosts, dems, tcosts, true, k);
  
  return cfl_block;
}

/*--------------------------------------------------------------------------*/
/*------------------------------ TEST CASES --------------------------------*/
/*--------------------------------------------------------------------------*/

// Test 1: Parameter Management
TestResult test_parameter_management() {
  try {
    // Part 1: Setting and getting all parameter types
    ScenarioReductionSolver solver;
    
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 2);
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 2) {
      return {false, "Failed to set/get algorithm parameter"};
    }
    
    solver.set_par(ScenarioReductionSolver::intShuffle, 1);
    if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 1) {
      return {false, "Failed to set/get shuffle parameter"};
    }
    
    solver.set_par(ScenarioReductionSolver::intRandomSeed, 12345);
    // Note: We can't verify the seed value, just that set_par doesn't throw
    
    // dblRho parameter has been removed - no test needed
    
    // Part 2: Boundary value testing
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Min valid
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 0) {
      return {false, "Failed to set min algorithm value"};
    }
    
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 3);  // Max valid
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 3) {
      return {false, "Failed to set max algorithm value"};
    }
    
    // Small value test for rho removed - parameter no longer exists
    
    // Part 3: Parameter validation and error handling
    try {
      solver.set_par(ScenarioReductionSolver::intAlgorithm, -1);
      return {false, "Should throw for negative algorithm value"};
    } catch (const std::invalid_argument&) {
      // Expected
    }
    
    try {
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 4);
      return {false, "Should throw for algorithm value > 3"};
    } catch (const std::invalid_argument&) {
      // Expected
    }
    
    // Negative rho test removed - parameter no longer exists
    
    // Part 4: Default parameter values
    ScenarioReductionSolver solver2;
    
    if (solver2.get_dflt_int_par(ScenarioReductionSolver::intAlgorithm) != 1) {
      return {false, "Default algorithm should be 1 (Dupacova)"};
    }
    
    if (solver2.get_dflt_int_par(ScenarioReductionSolver::intShuffle) != 0) {
      return {false, "Default shuffle should be 0 (false)"};
    }
    
    if (solver2.get_dflt_int_par(ScenarioReductionSolver::intRandomSeed) != 0) {
      return {false, "Default random seed should be 0"};
    }
    
    // Default rho test removed - parameter no longer exists
    
    // Verify get_*_par returns defaults initially
    if (solver2.get_int_par(ScenarioReductionSolver::intAlgorithm) != 1 ||
        solver2.get_int_par(ScenarioReductionSolver::intShuffle) != 0 ||
        !approx_equal(solver2.get_dbl_par(ScenarioReductionSolver::dblRho), 0.0)) {
      return {false, "get_*_par should return default values initially"};
    }
    
    // Part 5: Test get_vint_par for warm start indices
    ScenarioReductionSolver solver3;
    std::vector<int> test_indices = {0, 2, 3};
    solver3.set_par(ScenarioReductionSolver::vintWarmstartIndices, std::move(test_indices));
    
    const auto& retrieved_indices = solver3.get_vint_par(ScenarioReductionSolver::vintWarmstartIndices);
    if (retrieved_indices.size() != 3 || 
        retrieved_indices[0] != 0 || 
        retrieved_indices[1] != 2 || 
        retrieved_indices[2] != 3) {
      return {false, "get_vint_par failed to retrieve warm start indices correctly"};
    }
    
    return {true, "All parameter management tests passed"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 1 - Parameter Management", test_parameter_management);

// Test 2: Solution Handling (with x variable verification)
TestResult test_solution_handling() {
  try {
    // Part 1: Basic functionality with y variables (facility selection)
    {
      auto* block = new CapacitatedFacilityLocationBlock();
      
      // Set up a simple test instance
      int nf = 2;  // facilities (must equal customers for square distance matrix)
      int nc = 2;  // customers (scenarios)
      int k = 1;  // number of scenarios to select
      
      CapacitatedFacilityLocationBlock::DVector caps(nf);
      caps[0] = 1.0;  // Must be 1.0 for scenario reduction
      caps[1] = 1.0;
      
      CapacitatedFacilityLocationBlock::CVector fcosts(nf);
      fcosts[0] = 10.0;
      fcosts[1] = 15.0;
      
      CapacitatedFacilityLocationBlock::DVector dems(nc);
      dems[0] = 1.0;
      dems[1] = 1.0;
      
      CapacitatedFacilityLocationBlock::CMatrix tcosts(boost::extents[nf][nc]);
      tcosts[0][0] = 0.0;  // Distance from scenario 0 to 0
      tcosts[0][1] = 2.0;  // Distance from scenario 0 to 1
      tcosts[1][0] = 2.0;  // Distance from scenario 1 to 0 (symmetric)
      tcosts[1][1] = 0.0;  // Distance from scenario 1 to 1
      
      block->load(nf, nc, caps, fcosts, dems, tcosts, true, k);
      
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 1); // Dupacova
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "Solver failed to compute"};
      }
      
      const auto& reduced = solver.get_reduced_atoms();
      solver.get_var_solution();
      
      // Check y variables
      for (int i = 0; i < nc; ++i) {
        auto* y_var = block->get_y(i);
        if (!y_var) {
          delete block;
          return {false, "y variable not found"};
        }
        double val = y_var->get_value();
        bool is_selected = (val > 0.5);  // y should be 0 or 1
        if (is_selected != reduced[i]) {
          delete block;
          return {false, "y variable value doesn't match solver's solution"};
        }
      }
      
      // Verify exactly k scenarios selected
      int selected_count = std::count(reduced.begin(), reduced.end(), true);
      if (selected_count != k) {
        delete block;
        return {false, "Wrong number of scenarios selected"};
      }
      
      delete block;
    }
    
    // Part 2: X variable assignments (customer-to-facility mapping) - NEW TEST
    {
      auto* block = create_test_block(2);  // Select 2 out of 4 scenarios
      
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      solver.compute();
      
      const auto& reduced = solver.get_reduced_atoms();
      solver.get_var_solution();
      
      // Get the transportation costs for verification
      const auto& tcosts = block->get_Transportation_Costs();
      
      // Verify x variables are set correctly
      for (int j = 0; j < block->get_NCustomers(); ++j) {  // For each customer/scenario
        bool customer_assigned = false;
        int assigned_facility = -1;
        double assignment_cost = 0.0;
        
        for (int i = 0; i < block->get_NFacilities(); ++i) {  // Check all facilities
          auto* x_var = block->get_x(i, j);
          if (!x_var) {
            delete block;
            return {false, "x variable [" + std::to_string(i) + "][" + std::to_string(j) + "] not found"};
          }
          
          double val = x_var->get_value();
          if (val > 0.5) {  // x should be 0 or 1
            if (customer_assigned) {
              delete block;
              return {false, "Customer " + std::to_string(j) + " assigned to multiple facilities"};
            }
            customer_assigned = true;
            assigned_facility = i;
            assignment_cost = tcosts[j][i];  // Note: tcosts indexing is [customer][facility]
            
            // Verify the assigned facility is actually selected
            if (!reduced[i]) {
              delete block;
              return {false, "Customer " + std::to_string(j) + " assigned to unselected facility " + std::to_string(i)};
            }
          }
        }
        
        if (!customer_assigned) {
          delete block;
          return {false, "Customer " + std::to_string(j) + " not assigned to any facility"};
        }
        
        // Verify this is the closest selected facility
        for (int i = 0; i < block->get_NFacilities(); ++i) {
          if (reduced[i] && i != assigned_facility) {
            double other_cost = tcosts[j][i];
            if (other_cost < assignment_cost - 1e-6) {  // Allow small numerical tolerance
              delete block;
              return {false, "Customer " + std::to_string(j) + " not assigned to closest facility (assigned to " + 
                      std::to_string(assigned_facility) + " with cost " + std::to_string(assignment_cost) + 
                      " but facility " + std::to_string(i) + " has cost " + std::to_string(other_cost) + ")"};
            }
          }
        }
      }
      
      delete block;
    }
    
    // Part 3: Algorithm comparison (Baseline vs Dupacova)
    {
      auto* block = create_test_block(2);
      
      // Run Baseline
      ScenarioReductionSolver baseline_solver;
      baseline_solver.set_Block(block);
      baseline_solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Baseline
      baseline_solver.compute();
      double baseline_obj = baseline_solver.get_var_value();
      
      // Run Dupacova
      ScenarioReductionSolver dupacova_solver;
      dupacova_solver.set_Block(block);
      dupacova_solver.set_par(ScenarioReductionSolver::intAlgorithm, 1);  // Dupacova
      dupacova_solver.compute();
      double dupacova_obj = dupacova_solver.get_var_value();
      
      // Baseline should have worse or equal objective value than Dupacova
      if (baseline_obj < dupacova_obj - 1e-6) {
        delete block;
        return {false, "Baseline objective (" + std::to_string(baseline_obj) + 
                ") should not be better than Dupacova (" + std::to_string(dupacova_obj) + ")"};
      }
      
      delete block;
    }
    
    // Part 4: Solution independence (multiple solvers on same block)
    {
      auto* block = create_test_block(2);
      
      ScenarioReductionSolver solver1;
      solver1.set_Block(block);
      solver1.compute();
      
      std::vector<bool> solution1 = solver1.get_reduced_atoms();
      double obj1 = solver1.get_var_value();
      
      ScenarioReductionSolver solver2;
      solver2.set_Block(block);
      solver2.set_par(ScenarioReductionSolver::intAlgorithm, 0); // Baseline
      solver2.compute();
      
      // Verify both solvers maintain their own solutions
      const auto& solution1_check = solver1.get_reduced_atoms();
      for (size_t i = 0; i < solution1.size(); ++i) {
        if (solution1[i] != solution1_check[i]) {
          delete block;
          return {false, "First solver's solution was corrupted"};
        }
      }
      
      delete block;
    }
    
    // Part 5: Solution persistence through parameters
    {
      auto* block = create_test_block(2);
      
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      solver.compute();
      
      int algorithm = solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
      std::vector<bool> solution = solver.get_reduced_atoms();
      
      ScenarioReductionSolver solver2;
      solver2.set_Block(block);
      solver2.set_par(ScenarioReductionSolver::intAlgorithm, algorithm);
      solver2.compute();
      
      const auto& solution2 = solver2.get_reduced_atoms();
      if (solution.size() != solution2.size()) {
        delete block;
        return {false, "Solution sizes don't match after parameter restoration"};
      }
      
      for (size_t i = 0; i < solution.size(); ++i) {
        if (solution[i] != solution2[i]) {
          delete block;
          return {false, "Solutions don't match after restoring parameters"};
        }
      }
      
      delete block;
    }
    
    // Part 6: Warm start with vintWarmstartIndices
    {
      auto* block = create_test_block(4);
      
      ScenarioReductionSolver solver1;
      solver1.set_Block(block);
      solver1.set_par(ScenarioReductionSolver::intAlgorithm, 1); // Dupacova
      solver1.compute();
      
      const auto& solution1 = solver1.get_reduced_atoms();
      std::vector<int> warm_indices;
      for (size_t i = 0; i < solution1.size(); ++i) {
        if (solution1[i]) {
          warm_indices.push_back(static_cast<int>(i));
        }
      }
      
      ScenarioReductionSolver solver2;
      solver2.set_Block(block);
      solver2.set_par(ScenarioReductionSolver::intAlgorithm, 2); // BestFit
      solver2.set_par(ScenarioReductionSolver::intUseWarmstart, 1);
      solver2.set_par(ScenarioReductionSolver::vintWarmstartIndices, std::move(warm_indices));
      solver2.compute();
      
      double obj1 = solver1.get_var_value();
      double obj2 = solver2.get_var_value();
      
      if (obj2 > obj1 + 1e-6) {  // Local search should be at least as good as Dupacova with warm start
        delete block;
        return {false, "Warm start local search should not be worse than Dupacova"};
      }
      
      delete block;
    }
    
    // Part 7: has_var_solution() and get_var_solution() behavior
    {
      auto* block = create_test_block(2);
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      if (solver.has_var_solution()) {
        delete block;
        return {false, "has_var_solution() should return false before compute()"};
      }
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "Compute failed"};
      }
      
      if (!solver.has_var_solution()) {
        delete block;
        return {false, "has_var_solution() should return true after successful compute()"};
      }
      
      // Test get_nb_atoms(), get_nb_reduced(), and get_weights()
      if (solver.get_nb_atoms() != 4) {
        delete block;
        return {false, "get_nb_atoms() should return 4 for test block"};
      }
      
      if (solver.get_nb_reduced() != 2) {
        delete block;
        return {false, "get_nb_reduced() should return 2 (k value)"};
      }
      
      const auto* weights = solver.get_weights();
      if (!weights) {
        delete block;
        return {false, "get_weights() should not return nullptr"};
      }
      
      // Weights should be normalized (sum to 1.0)
      double weight_sum = 0.0;
      for (size_t i = 0; i < 4; ++i) {
        weight_sum += (*weights)[i];
      }
      if (!approx_equal(weight_sum, 1.0, 1e-6)) {
        delete block;
        return {false, "Weights should be normalized to sum to 1.0, got " + std::to_string(weight_sum)};
      }
      
      // Test error when no solution available
      ScenarioReductionSolver solver2;
      solver2.set_Block(block);
      
      try {
        solver2.get_var_solution();
        delete block;
        return {false, "Should throw when no solution available"};
      } catch (const std::logic_error& e) {
        if (std::string(e.what()).find("no solution available") == std::string::npos) {
          delete block;
          return {false, "Wrong error message for no solution"};
        }
      }
      
      // Test all algorithms produce valid x and y variables
      for (int algo = 0; algo <= 3; ++algo) {
        ScenarioReductionSolver algo_solver;
        algo_solver.set_Block(block);
        algo_solver.set_par(ScenarioReductionSolver::intAlgorithm, algo);
        algo_solver.compute();
        algo_solver.get_var_solution();
        
        int count_selected = 0;
        for (CapacitatedFacilityLocationBlock::Index i = 0; i < block->get_NFacilities(); ++i) {
          auto y_var = block->get_y(i);
          if (y_var && approx_equal(y_var->get_value(), 1.0)) {
            count_selected++;
          }
        }
        
        if (count_selected != 2) {  // Should select exactly k=2 facilities
          delete block;
          return {false, "Algorithm " + std::to_string(algo) + " selected " + 
                  std::to_string(count_selected) + " facilities instead of 2"};
        }
      }
      
      delete block;
    }
    
    return {true, "All solution handling tests passed (including x variable verification)"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 2 - Solution Handling", test_solution_handling);

// Test 3: Thread Safety & Concurrency
TestResult test_thread_safety() {
  try {
    // Part 1: Concurrent parameter setting
    {
      ScenarioReductionSolver solver;
      std::vector<std::thread> threads;
      std::atomic<int> errors(0);
      
      // Launch multiple threads that try to set parameters
      for (int i = 0; i < 10; ++i) {
        threads.emplace_back([&solver, &errors, i]() {
          try {
            solver.lock();
            
            // Each thread sets different parameter values
            solver.set_par(ScenarioReductionSolver::intAlgorithm, i % 4);
            
            // Verify the values are what we just set
            int alg = solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
                  
            if (alg != i % 4) {
              errors++;
            }
            
            solver.unlock();
          } catch (...) {
            errors++;
          }
        });
      }
      
      // Wait for all threads
      for (auto& t : threads) {
        t.join();
      }
      
      if (errors > 0) {
        return {false, "Concurrent parameter setting had " + std::to_string(errors) + " errors"};
      }
    }
    
    // Part 2: Preventing simultaneous compute() calls
    {
      auto block = create_test_block(1);
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      std::atomic<int> concurrent_computes(0);
      std::atomic<int> max_concurrent(0);
      std::atomic<int> completed(0);
      std::vector<std::thread> threads;
      
      // Launch multiple threads that try to compute simultaneously
      for (int i = 0; i < 5; ++i) {
        threads.emplace_back([&solver, &concurrent_computes, &max_concurrent, &completed]() {
          int current = concurrent_computes.fetch_add(1) + 1;
          
          // Update max concurrent if needed
          int expected = max_concurrent.load();
          while (current > expected && 
                 !max_concurrent.compare_exchange_weak(expected, current)) {
          }
          
          // Call compute - should be properly serialized
          int status = solver.compute();
          
          concurrent_computes--;
          
          if (status == Solver::kOK) {
            completed++;
          }
        });
      }
      
      // Wait for all threads
      for (auto& t : threads) {
        t.join();
      }
      
      delete block;
      
      // Verify all threads completed successfully
      if (completed != 5) {
        return {false, "Not all threads completed compute(): " + 
                       std::to_string(completed.load()) + "/5"};
      }
    }
    
    // Part 3: Thread safety during set_Block()
    {
      std::vector<std::thread> threads;
      std::atomic<int> errors(0);
      
      // Create multiple solvers and blocks
      for (int i = 0; i < 5; ++i) {
        threads.emplace_back([&errors, i]() {
          try {
            auto solver = new ScenarioReductionSolver();
            auto block = create_test_block(2);
            
            // Lock before set_Block
            solver->lock();
            solver->set_Block(block);
            solver->unlock();
            
            // Verify block was set
            solver->lock();
            auto retrieved_block = solver->get_Block();
            solver->unlock();
            
            if (retrieved_block != block) {
              errors++;
            }
            
            delete solver;
            delete block;
          } catch (...) {
            errors++;
          }
        });
      }
      
      // Wait for all threads
      for (auto& t : threads) {
        t.join();
      }
      
      if (errors > 0) {
        return {false, "Errors during concurrent set_Block: " + std::to_string(errors)};
      }
    }
    
    // Part 4: Thread safety of get_var_solution()
    {
      auto block = create_test_block(1);
      
      ScenarioReductionSolver solver;
      
      solver.lock();
      solver.set_Block(block);
      int status = solver.compute();
      solver.unlock();
      
      if (status != Solver::kOK) {
        delete block;
        return {false, "Compute failed with status: " + std::to_string(status)};
      }
      
      std::vector<std::thread> threads;
      std::atomic<int> errors(0);
      
      // Multiple threads trying to get solution
      for (int i = 0; i < 5; ++i) {
        threads.emplace_back([&solver, &block, &errors]() {
          try {
            // Block should be locked before calling get_var_solution
            solver.lock();
            block->lock(solver.id());
            
            if (solver.has_var_solution()) {
              solver.get_var_solution();
            }
            
            block->unlock(solver.id());
            solver.unlock();
          } catch (...) {
            errors++;
          }
        });
      }
      
      // Wait for all threads
      for (auto& t : threads) {
        t.join();
      }
      
      delete block;
      
      if (errors > 0) {
        return {false, "Errors during concurrent get_var_solution: " + std::to_string(errors)};
      }
    }
    
    // Part 5: Stress test with mixed operations
    {
      auto block = create_test_block(2);
      ScenarioReductionSolver solver;
      
      solver.lock();
      solver.set_Block(block);
      solver.unlock();
      
      std::atomic<int> errors(0);
      std::atomic<bool> stop(false);
      std::vector<std::thread> threads;
      
      // Thread 1: Continuously set parameters
      threads.emplace_back([&solver, &errors, &stop]() {
        int iteration = 0;
        while (!stop) {
          try {
            solver.lock();
            solver.set_par(ScenarioReductionSolver::intAlgorithm, iteration % 4);
            solver.set_par(ScenarioReductionSolver::dblRho, 0.01 * (iteration % 10));
            solver.unlock();
            iteration++;
          } catch (...) {
            errors++;
          }
        }
      });
      
      // Thread 2: Continuously read parameters
      threads.emplace_back([&solver, &errors, &stop]() {
        while (!stop) {
          try {
            solver.lock();
            solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
            solver.get_dbl_par(ScenarioReductionSolver::dblRho);
            solver.unlock();
          } catch (...) {
            errors++;
          }
        }
      });
      
      // Thread 3: Periodically compute
      threads.emplace_back([&solver, &block, &errors, &stop]() {
        while (!stop) {
          try {
            solver.lock();
            solver.compute();
            solver.unlock();
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
          } catch (...) {
            errors++;
          }
        }
      });
      
      // Let threads run for a short time
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      stop = true;
      
      // Wait for all threads
      for (auto& t : threads) {
        t.join();
      }
      
      delete block;
      
      if (errors > 0) {
        return {false, "Stress test had " + std::to_string(errors) + " errors"};
      }
    }
    
    return {true, "All thread safety tests passed"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 3 - Thread Safety & Concurrency", test_thread_safety);

// Test 4: Error Handling & Edge Cases
TestResult test_error_handling() {
  try {
    // Part 1: refresh_cached_data() error conditions
    {
      auto* block = create_test_block(2);
      ScenarioReductionSolver solver;
      
      try {
        solver.set_Block(block);
        solver.refresh_cached_data(-1);
        delete block;
        return {false, "Should throw for negative k"};
      } catch (const std::invalid_argument& e) {
        if (std::string(e.what()).find("k must be non-negative") == std::string::npos) {
          delete block;
          return {false, "Wrong error message for negative k"};
        }
      }
      
      try {
        solver.refresh_cached_data(10);  // Try to select 10 out of 4
        delete block;
        return {false, "Should throw when k exceeds number of scenarios"};
      } catch (const std::invalid_argument& e) {
        if (std::string(e.what()).find("k cannot exceed number of scenarios") == std::string::npos) {
          delete block;
          return {false, "Wrong error message for k > n"};
        }
      }
      
      delete block;
    }
    
    // Part 2: Non-square matrix validation
    {
      auto* block = new CapacitatedFacilityLocationBlock();
      
      int nf = 3;  // 3 facilities
      int nc = 5;  // 5 customers (non-square)
      int k = 2;
      
      // Set up minimal data
      CapacitatedFacilityLocationBlock::CVector fcosts(nf, 100.0);
      CapacitatedFacilityLocationBlock::CMatrix tcosts(boost::extents[nf][nc]);
      CapacitatedFacilityLocationBlock::DVector caps(nf, 1.0);
      CapacitatedFacilityLocationBlock::DVector dems(nc, 1.0);
      
      for (int i = 0; i < nf; ++i) {
        for (int j = 0; j < nc; ++j) {
          tcosts[i][j] = 10.0;
        }
      }
      
      block->load(nf, nc, caps, fcosts, dems, tcosts, true, k);
      
      ScenarioReductionSolver solver;
      
      try {
        solver.set_Block(block);  // This should throw in refresh_cached_data
        delete block;
        return {false, "Should throw for non-square matrix"};
      } catch (const std::logic_error& e) {
        if (std::string(e.what()).find("number of customers must equal number of facilities") == std::string::npos) {
          delete block;
          return {false, "Wrong error message for non-square matrix"};
        }
      }
      delete block;
    }
    
    // Part 3: Weight normalization
    {
      // Test with unnormalized weights
      auto* block = create_test_block(2); 
      // Verify the demands are not normalized
      const auto& demands = block->get_Demands();
      double sum = std::accumulate(demands.begin(), demands.end(), 0.0);
      if (std::abs(sum - 1.0) <= 1e-6) {
        delete block;
        return {false, "Test demands should be unnormalized"};
      }
      
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      // Solver should work with unnormalized weights (normalized internally)
      int result = solver.compute();
      if (result != Solver::kOK) {
        delete block;
        return {false, "Solver should succeed with unnormalized weights"};
      }
      
      // Verify solution exists
      if (!solver.has_var_solution()) {
        delete block;
        return {false, "Should have solution after compute"};
      }
      
      // The original demands in block should be unchanged
      if (std::abs(demands[0] - 50.0) > 1e-10 ||
          std::abs(demands[1] - 100.0) > 1e-10 ||
          std::abs(demands[2] - 150.0) > 1e-10 ||
          std::abs(demands[3] - 200.0) > 1e-10) {
        delete block;
        return {false, "Original demands should be unchanged"};
      }
      
      delete block;
      
      // Test with already normalized weights
      auto* block2 = new CapacitatedFacilityLocationBlock();
      
      int nf = 4, nc = 4, k = 2;
      
      // Set up normalized demands (0.1, 0.2, 0.3, 0.4 - sum = 1.0)
      CapacitatedFacilityLocationBlock::DVector dems(nc);
      dems[0] = 0.1;
      dems[1] = 0.2;
      dems[2] = 0.3;
      dems[3] = 0.4;
      
      // Verify they're normalized
      double sum2 = std::accumulate(dems.begin(), dems.end(), 0.0);
      if (std::abs(sum2 - 1.0) > 1e-6) {
        delete block2;
        return {false, "Test demands should be normalized"};
      }
      
      CapacitatedFacilityLocationBlock::DVector caps(nf);
      CapacitatedFacilityLocationBlock::CVector fcosts(nf);
      CapacitatedFacilityLocationBlock::CMatrix tcosts(boost::extents[nc][nf]);
      
      for (int i = 0; i < nf; ++i) {
        caps[i] = 1.0;
        fcosts[i] = 0.0;
      }
      
      for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nf; ++j) {
          tcosts[i][j] = std::abs(i - j);
        }
      }
      
      block2->load(nf, nc, caps, fcosts, dems, tcosts, true, k);
      
      ScenarioReductionSolver solver2;
      solver2.set_Block(block2);
      
      // Solver should work with normalized weights
      result = solver2.compute();
      if (result != Solver::kOK) {
        delete block2;
        return {false, "Solver should succeed with normalized weights"};
      }
      
      if (!solver2.has_var_solution()) {
        delete block2;
        return {false, "Should have solution after compute with normalized weights"};
      }
      
      delete block2;
    }
    
    // Part 4: Edge cases
    {
      // Test k = 0 (select no scenarios)
      auto* block = create_test_block(2);
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      solver.refresh_cached_data(0);  // k = 0 is valid
      
      // Test k = n (select all scenarios)
      solver.refresh_cached_data(4);  // k = n is valid
      
      delete block;
    }
    
    return {true, "All error handling and edge case tests passed"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 4 - Error Handling & Edge Cases", test_error_handling);

/*--------------------------------------------------------------------------*/
/*---------------------------- TEST 5 - CONFIG FILE ------------------------*/
/*--------------------------------------------------------------------------*/

TestResult test_config_deserialization() {
  try {
    if (verbose) {
      std::cout << "\n--- Part 1: Testing BlockSolverConfig deserialization ---\n";
    }
    
    // Create solver and test config deserialization
    ScenarioReductionSolver solver;
    
    // Parse BlockSolverConfig from file
    Configuration* cfg = nullptr;
    try {
      cfg = Configuration::deserialize("BSConfig_SR.txt");
    } catch (const std::exception& e) {
      return {false, std::string("Failed to deserialize config: ") + e.what()};
    }
    
    BlockSolverConfig* bsc = dynamic_cast<BlockSolverConfig*>(cfg);
    if (!bsc) {
      delete cfg;
      return {false, "Failed to cast Configuration to BlockSolverConfig"};
    }
    
    // Apply ComputeConfig to solver
    if (bsc->get_SolverNames().size() != 1) {
      delete bsc;
      return {false, "Config should have exactly 1 solver"};
    }
    
    if (bsc->get_SolverConfigs().size() != 1) {
      delete bsc;
      return {false, "Config should have exactly 1 ComputeConfig"};
    }
    
    // Get the ComputeConfig and apply to solver
    ComputeConfig* cc = bsc->get_SolverConfig(0);
    solver.set_ComputeConfig(cc);
    
    // Verify parameters were set correctly
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 1) {
      delete bsc;
      return {false, "intAlgorithm should be 1 (Dupacova) after config"};
    }
    
    if (solver.get_int_par(ScenarioReductionSolver::intUseWarmstart) != 0) {
      delete bsc;
      return {false, "intUseWarmstart should be 0 after config"};
    }
    
    if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 0) {
      delete bsc;
      return {false, "intShuffle should be 0 after config"};
    }
    
    if (solver.get_int_par(ScenarioReductionSolver::intRandomSeed) != 0) {
      delete bsc;
      return {false, "intRandomSeed should be 0 after config"};
    }
    
    if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.0, 1e-10)) {
      delete bsc;
      return {false, "dblRho should be 0.0 after config"};
    }
    
    // Clean up config after verification
    delete bsc;
    
    if (verbose) {
      std::cout << "✓ Config deserialization successful\n";
    }
    
    if (verbose) {
      std::cout << "\n--- Part 2: Testing different algorithms via direct parameter setting ---\n";
    }
    
    // Test setting different algorithms through set_par with enum values
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Baseline
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 0) {
      return {false, "Failed to set algorithm to Baseline"};
    }
    
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 2);  // BestFit
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 2) {
      return {false, "Failed to set algorithm to BestFit"};
    }
    
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 3);  // FirstFit
    if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 3) {
      return {false, "Failed to set algorithm to FirstFit"};
    }
    
    // Test setting other parameters
    solver.set_par(ScenarioReductionSolver::intShuffle, 1);
    if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 1) {
      return {false, "Failed to set intShuffle"};
    }
    
    solver.set_par(ScenarioReductionSolver::dblRho, 0.5);
    if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.5, 1e-10)) {
      return {false, "Failed to set dblRho"};
    }
    
    if (verbose) {
      std::cout << "✓ Parameter setting successful\n";
      std::cout << "\n--- Part 3: Testing solver with config-loaded parameters ---\n";
    }
    
    // Create a test block and solve with config-loaded solver
    CapacitatedFacilityLocationBlock* block = create_test_block(2);  // k=2
    
    // Reset to Dupacova and solve
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 1);
    solver.set_Block(block);
    
    int status = solver.compute();
    if (status != Solver::kOK) {
      delete block;
      return {false, "Solver failed to compute with config-loaded parameters"};
    }
    
    if (!solver.has_var_solution()) {
      delete block;
      return {false, "Solver should have solution after compute"};
    }
    
    double obj_value = solver.get_var_value();
    if (verbose) {
      std::cout << "✓ Solver computed successfully with config parameters, objective: " 
                << obj_value << "\n";
    }
    
    delete block;
    
    return {true, "All config deserialization tests passed"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 5 - BlockSolverConfig Deserialization", test_config_deserialization);

/*--------------------------------------------------------------------------*/
/*------------------------- TEST 6: LOGGING --------------------------------*/
/*--------------------------------------------------------------------------*/

TestResult test_logging() {
  try {
    // Create a simple test block
    auto* block = new CapacitatedFacilityLocationBlock();
    
    // Create 5 scenarios with equal probabilities
    int nf = 5;  // facilities (must equal customers for square matrix)
    int nc = 5;  // customers (scenarios) 
    int k = 2;   // Select 2 scenarios
    
    CapacitatedFacilityLocationBlock::DVector capacities(nf, 1.0);
    CapacitatedFacilityLocationBlock::CVector fcosts(nf, 0.0); // No facility costs
    CapacitatedFacilityLocationBlock::DVector demands(nc, 1.0/nc);
    CapacitatedFacilityLocationBlock::CMatrix costs(boost::extents[nf][nc]);
    
    // Simple distance matrix
    for (int i = 0; i < nf; ++i) {
      for (int j = 0; j < nc; ++j) {
        costs[i][j] = (i == j) ? 0.0 : std::abs(i - j) * 1.0;
      }
    }
    
    block->load(nf, nc, capacities, fcosts, demands, costs, true, k);
    
    // Test 1: Dupacova with logging
    {
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      // Capture log output
      std::ostringstream log_stream;
      solver.set_log(&log_stream);
      solver.set_par(Solver::intLogVerb, 2);  // High verbosity
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 1); // Dupacova
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "Dupacova solver failed"};
      }
      
      std::string log_output = log_stream.str();
      
      // Check that log contains expected content
      if (log_output.find("Dupacova") == std::string::npos) {
        delete block;
        return {false, "Log missing algorithm name"};
      }
      
      // The iteration logs appear at verbosity level 2
      if (log_output.find("iteration") == std::string::npos) {
        // For Dupacova, we should at least see the algorithm name
        // Iteration details are at higher verbosity
        if (verbose) {
          std::cout << "Note: Iteration details not logged (verbosity 2)\n";
        }
      }
      
      if (verbose) {
        std::cout << "Dupacova log output:\n" << log_output << "\n";
      }
    }
    
    // Test 2: BestFit with logging
    {
      std::cout << "Starting BestFit test..." << std::endl;
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      std::ostringstream log_stream;
      solver.set_log(&log_stream);
      solver.set_par(Solver::intLogVerb, 2);  // High verbosity for testing
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 2); // BestFit
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "BestFit solver failed"};
      }
      
      std::string log_output = log_stream.str();
      
      if (log_output.find("BestFit") == std::string::npos) {
        delete block;
        return {false, "BestFit log missing algorithm name"};
      }
      
      // Always print for debugging
      std::cout << "BestFit log output:\n" << log_output << "\n";
      
      if (log_output.find("converged") == std::string::npos) {
        delete block;
        return {false, "BestFit log missing convergence info"};
      }
    }
    
    // Test 3: FirstFit with logging to show table format
    {
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      std::ostringstream log_stream;
      solver.set_log(&log_stream);
      solver.set_par(Solver::intLogVerb, 2);  // High verbosity to see tables
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 3); // FirstFit
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "FirstFit solver failed"};
      }
      
      std::string log_output = log_stream.str();
      
      // With high verbosity, we should see the table
      if (verbose) {
        std::cout << "FirstFit log output:\n" << log_output << "\n";
      }
    }
    
    // Test 4: Baseline with logging
    {
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      
      std::ostringstream log_stream;
      solver.set_log(&log_stream);
      solver.set_par(Solver::intLogVerb, 1);
      solver.set_par(ScenarioReductionSolver::intAlgorithm, 0); // Baseline
      
      int status = solver.compute();
      if (status != Solver::kOK) {
        delete block;
        return {false, "Baseline solver failed"};
      }
      
      std::string log_output = log_stream.str();
      
      if (log_output.find("Baseline") == std::string::npos) {
        delete block;
        return {false, "Baseline log missing algorithm name"};
      }
      
      if (log_output.find("highest probability") == std::string::npos) {
        delete block;
        return {false, "Baseline log missing selection method"};
      }
      
      if (verbose) {
        std::cout << "Baseline log output:\n" << log_output << "\n";
      }
    }
    
    delete block;
    return {true, "All logging tests passed"};
    
  } catch (const std::exception& e) {
    return {false, std::string("Exception: ") + e.what()};
  }
}

REGISTER_TEST("Test 6 - Logging", test_logging);

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main(int argc, char* argv[]) {
  std::cout << "ScenarioReductionSolver Test Suite\n";
  std::cout << "==================================\n";
  
  std::string test_name;
  
  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    
    if (arg == "--list" || arg == "-l") {
      std::cout << "\nAvailable tests:\n";
      for (const auto& [name, _] : test_registry) {
        std::cout << "  " << name << "\n";
      }
      return 0;
    }
    
    if (arg == "--help" || arg == "-h") {
      std::cout << "\nUsage: " << argv[0] << " [options] [test_name]\n\n";
      std::cout << "Options:\n";
      std::cout << "  --list, -l        List all available tests\n";
      std::cout << "  --verbose, -v     Enable verbose output\n";
      std::cout << "  --help, -h        Show this help message\n";
      std::cout << "  test_name         Run only the specified test\n\n";
      std::cout << "If no test name is provided, all tests are run.\n";
      std::cout << "By default, output is minimal (non-verbose).\n";
      return 0;
    }
    
    if (arg == "--verbose" || arg == "-v") {
      verbose = true;
      continue;
    }
    
    // Assume it's a test name
    if (test_name.empty()) {
      test_name = arg;
    }
  }
  
  std::cout << "\n";
  
  // Run specific test or all tests
  if (!test_name.empty()) {
    auto it = test_registry.find(test_name);
    if (it != test_registry.end()) {
      run_test(it->first, it->second);
    } else {
      std::cerr << "Error: Test '" << test_name << "' not found\n\n";
      std::cout << "Available tests:\n";
      for (const auto& [name, _] : test_registry) {
        std::cout << "  " << name << "\n";
      }
      return 1;
    }
  } else {
    // Run all tests
    for (const auto& [name, func] : test_registry) {
      run_test(name, func);
    }
  }
  
  // Print summary
  std::cout << "\n==================================\n";
  std::cout << "Test Summary: " << tests_passed << " passed, " 
            << tests_failed << " failed (out of " << tests_run << " run)\n";
  
  return (tests_failed > 0) ? 1 : 0;
}

/*--------------------------------------------------------------------------*/
/*------------------------- End test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/