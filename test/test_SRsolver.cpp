/*--------------------------------------------------------------------------*/
/*------------------------ File test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tests for ScenarioReductionSolver class.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * Copyright &copy; by Benoît Tran
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- OVERVIEW ------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
 * TEST SUITE OVERVIEW
 * ===================
 * 
 * This test suite contains 15 tests organized by theme:
 * 
 * PARAMETER MANAGEMENT (3 tests)
 * - parameter_setting: Tests setting/getting all parameter types and boundary values
 * - parameter_validation: Tests parameter validation and error handling
 * - default_values: Tests default parameter values and get_dflt_*_par methods
 * 
 * SOLUTION HANDLING (6 tests)
 * - solution_basic_functionality: Tests writing solutions to Block's y variables and algorithm correctness
 * - solution_independence: Tests multiple solvers working independently
 * - solution_persistence: Tests saving/restoring solver state through parameters
 * - solution_warm_start: Tests warm start functionality with vintWarmstartIndices
 * - has_var_solution_test: Tests has_var_solution() method behavior
 * - get_var_solution_test: Tests get_var_solution() method and error handling
 * 
 * THREAD SAFETY & CONCURRENCY (5 tests)
 * - concurrent_parameter_setting: Tests thread-safe parameter access
 * - simultaneous_compute_prevention: Tests preventing concurrent compute() calls
 * - lock_during_set_block: Tests thread safety during set_Block()
 * - get_var_solution_locking: Tests thread safety of get_var_solution()
 * - stress_test_mixed_operations: Stress test with multiple concurrent operations
 * 
 * ERROR HANDLING & EDGE CASES (1 test)
 * - refresh_cached_data_error_handling: Tests error handling in refresh_cached_data()
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <thread>
#include <chrono>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <algorithm>

#include "SMSTypedefs.h"
#include "Configuration.h"
#include "ScenarioReductionSolver.h"
#include "CapacitatedFacilityLocationBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- USING ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- TEST REGISTRY --------------------------------*/
/*--------------------------------------------------------------------------*/

class TestRegistry {
public:
  using TestFunction = std::function<void()>;
  
  static TestRegistry& instance() {
    static TestRegistry registry;
    return registry;
  }
  
  void register_test(const std::string& name, TestFunction test) {
    tests_[name] = test;
  }
  
  void run_all() {
    std::cout << "Running " << tests_.size() << " tests...\n\n";
    int passed = 0;
    int failed = 0;
    
    for (const auto& [name, test] : tests_) {
      std::cout << "Running test: " << name << " ... ";
      std::cout.flush();
      
      try {
        test();
        std::cout << "PASSED\n";
        ++passed;
      } catch (const std::exception& e) {
        std::cout << "FAILED\n";
        std::cerr << "  Error: " << e.what() << "\n";
        ++failed;
      }
    }
    
    std::cout << "\nTest summary: " << passed << " passed, " 
              << failed << " failed\n";
    
    if (failed > 0) {
      exit(1);
    }
  }
  
  void run_test(const std::string& name) {
    auto it = tests_.find(name);
    if (it == tests_.end()) {
      std::cerr << "Test '" << name << "' not found\n";
      list_tests();
      exit(1);
    }
    
    std::cout << "Running test: " << name << " ... ";
    std::cout.flush();
    
    try {
      it->second();
      std::cout << "PASSED\n";
    } catch (const std::exception& e) {
      std::cout << "FAILED\n";
      std::cerr << "  Error: " << e.what() << "\n";
      exit(1);
    }
  }
  
  void list_tests() {
    std::cout << "Available tests:\n";
    for (const auto& [name, _] : tests_) {
      std::cout << "  " << name << "\n";
    }
  }
  
private:
  std::map<std::string, TestFunction> tests_;
};

#define REGISTER_TEST(name) \
  static void test_##name(); \
  static bool register_##name = []() { \
    TestRegistry::instance().register_test(#name, test_##name); \
    return true; \
  }(); \
  static void test_##name()

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
  
  auto block = new CapacitatedFacilityLocationBlock();
  
  // Create a small problem: for scenario reduction, we need a square distance matrix
  // If we interpret customers as scenarios, we need nc x nc matrix
  int nf = 4;  // Must equal nc for square matrix
  int nc = 4;
  
  // Set up facility costs
  CapacitatedFacilityLocationBlock::CVector fcosts(nf);
  for (int i = 0; i < nf; ++i) {
    fcosts[i] = 100.0 * (i + 1);  // 100, 200, 300
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
  
  // Set up capacities (interpreted as probabilities)
  CapacitatedFacilityLocationBlock::DVector caps(nf);
  caps[0] = 0.25;
  caps[1] = 0.25;
  caps[2] = 0.25;
  caps[3] = 0.25;
  
  // Set up demands
  CapacitatedFacilityLocationBlock::DVector dems(nc);
  for (int i = 0; i < nc; ++i) {
    dems[i] = 50.0 * (i + 1);  // 50, 100, 150, 200
  }
  
  // Load data into block with k as the max number of facilities
  block->load(nf, nc, caps, fcosts, dems, tcosts, false, k);
  
  return block;
}

/*--------------------------------------------------------------------------*/
/*------------------------------ TEST CASES --------------------------------*/
/*--------------------------------------------------------------------------*/

REGISTER_TEST(parameter_setting) {
  ScenarioReductionSolver solver;
  
  // Test setting and getting integer parameters
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 2);
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 2) {
    throw std::runtime_error("Failed to set algorithm parameter");
  }
  
  solver.set_par(ScenarioReductionSolver::intShuffle, 1);
  if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 1) {
    throw std::runtime_error("Failed to set shuffle parameter");
  }
  
  // Note: We can't verify the random seed value since mt19937 doesn't expose it
  // But we can test that set_par doesn't throw
  solver.set_par(ScenarioReductionSolver::intRandomSeed, 12345);
  
  // Test setting and getting double parameters
  solver.set_par(ScenarioReductionSolver::dblEll, 0.5);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 0.5)) {
    throw std::runtime_error("Failed to set ell parameter");
  }
  
  solver.set_par(ScenarioReductionSolver::dblRho, 0.01);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.01)) {
    throw std::runtime_error("Failed to set rho parameter");
  }
  
  // Test boundary values for algorithm
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Min valid
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 0) {
    throw std::runtime_error("Failed to set min algorithm value");
  }
  
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 3);  // Max valid
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 3) {
    throw std::runtime_error("Failed to set max algorithm value");
  }
  
  // Test small positive values for double parameters
  solver.set_par(ScenarioReductionSolver::dblEll, 1e-10);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 1e-10)) {
    throw std::runtime_error("Failed to set small ell value");
  }
  
  solver.set_par(ScenarioReductionSolver::dblRho, 1e-10);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 1e-10)) {
    throw std::runtime_error("Failed to set small rho value");
  }
}

REGISTER_TEST(parameter_validation) {
  ScenarioReductionSolver solver;
  
  // Test invalid parameter values
  try {
    solver.set_par(ScenarioReductionSolver::intAlgorithm, -1);
    throw std::runtime_error("Should have thrown for invalid algorithm value");
  } catch (const std::invalid_argument&) {
    // Expected
  }
  
  try {
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 4);
    throw std::runtime_error("Should have thrown for invalid algorithm value");
  } catch (const std::invalid_argument&) {
    // Expected
  }
  
  // Note: shuffle parameter accepts any integer (non-zero = true, zero = false)
  // So there are no invalid values to test
  
  try {
    solver.set_par(ScenarioReductionSolver::dblEll, -1.0);
    throw std::runtime_error("Should have thrown for invalid ell value");
  } catch (const std::invalid_argument&) {
    // Expected
  }
  
  try {
    solver.set_par(ScenarioReductionSolver::dblRho, -1.0);
    throw std::runtime_error("Should have thrown for invalid rho value");
  } catch (const std::invalid_argument&) {
    // Expected
  }
}

REGISTER_TEST(default_values) {
  ScenarioReductionSolver solver;
  
  // Test getting default values through get_dflt_*_par
  if (solver.get_dflt_int_par(ScenarioReductionSolver::intAlgorithm) != 1) {
    throw std::runtime_error("Default algorithm should be 1");
  }
  
  if (solver.get_dflt_int_par(ScenarioReductionSolver::intShuffle) != 0) {
    throw std::runtime_error("Default shuffle should be 0");
  }
  
  if (solver.get_dflt_int_par(ScenarioReductionSolver::intRandomSeed) != 0) {
    throw std::runtime_error("Default random seed should be 0");
  }
  
  if (!approx_equal(solver.get_dflt_dbl_par(ScenarioReductionSolver::dblEll), 2.0)) {
    throw std::runtime_error("Default ell should be 2.0");
  }
  
  if (!approx_equal(solver.get_dflt_dbl_par(ScenarioReductionSolver::dblRho), 0.0)) {
    throw std::runtime_error("Default rho should be 0.0");
  }
  
  // Also verify that get_*_par returns the same defaults initially
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 1 ||
      solver.get_int_par(ScenarioReductionSolver::intShuffle) != 0 ||
      !approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 2.0) ||
      !approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.0)) {
    throw std::runtime_error("get_*_par should return default values initially");
  }
}


/*--------------------------------------------------------------------------*/
/*------------------------ THREAD SAFETY TESTS -----------------------------*/
/*--------------------------------------------------------------------------*/

REGISTER_TEST(concurrent_parameter_setting) {
  ScenarioReductionSolver solver;
  std::vector<std::thread> threads;
  std::atomic<int> errors(0);
  
  // Launch multiple threads that try to set parameters
  for (int i = 0; i < 10; ++i) {
    threads.emplace_back([&solver, &errors, i]() {
      try {
        solver.lock();
        
        // Each thread sets different parameter values
        solver.set_par(ScenarioReductionSolver::intAlgorithm, i % 4);  // Only 0-3 are valid
        solver.set_par(ScenarioReductionSolver::dblRho, 0.1 * i);
        
        // Verify the values are what we just set
        int alg = solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
        double rho = solver.get_dbl_par(ScenarioReductionSolver::dblRho);
        
        if (alg != i % 4) {
          errors++;
        }
        if (std::abs(rho - 0.1 * i) > 1e-6) {
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
    throw std::runtime_error("Concurrent parameter setting had " + 
                            std::to_string(errors) + " errors");
  }
}

REGISTER_TEST(simultaneous_compute_prevention) {
  // Test that compute() calls are properly serialized by the mutex
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
      // Increment concurrent counter
      int current = concurrent_computes.fetch_add(1) + 1;
      
      // Update max concurrent if needed
      int expected = max_concurrent.load();
      while (current > expected && 
             !max_concurrent.compare_exchange_weak(expected, current)) {
      }
      
      // Call compute - this should block on the mutex if another compute is running
      int status = solver.compute();
      
      // Decrement concurrent counter
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
  
  // Clean up
  delete block;
  
  // Verify that compute calls were serialized (max 1 concurrent)
  if (max_concurrent > 1) {
    throw std::runtime_error("Multiple threads were in compute() simultaneously: " + 
                           std::to_string(max_concurrent.load()));
  }
  
  // Verify all threads completed
  if (completed != 5) {
    throw std::runtime_error("Not all threads completed compute(): " + 
                           std::to_string(completed.load()) + "/5");
  }
}

REGISTER_TEST(lock_during_set_block) {
  std::vector<std::thread> threads;
  std::atomic<int> errors(0);
  
  // Create multiple solvers and blocks
  for (int i = 0; i < 5; ++i) {
    threads.emplace_back([&errors, i]() {
      try {
        auto solver = new ScenarioReductionSolver();
        auto block = create_test_block(2);  // Test with 2 out of 3 facilities
        
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
    throw std::runtime_error("Errors during concurrent set_Block: " + 
                            std::to_string(errors));
  }
}

REGISTER_TEST(get_var_solution_locking) {
  auto block = create_test_block(1);  // Reduce to 1 facility
  
  ScenarioReductionSolver solver;
  
  solver.lock();
  solver.set_Block(block);
  int status = solver.compute();
  solver.unlock();
  
  if (status != Solver::kOK) {
    delete block;
    throw std::runtime_error("Compute failed with status: " + 
                            std::to_string(status));
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
    throw std::runtime_error("Errors during concurrent get_var_solution: " + 
                            std::to_string(errors));
  }
}

REGISTER_TEST(stress_test_mixed_operations) {
  auto block = create_test_block(2);  // Use 2 facilities for stress test
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
        solver.set_par(ScenarioReductionSolver::intAlgorithm, iteration % 4);  // Only 0-3 are valid
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
    throw std::runtime_error("Stress test had " + 
                            std::to_string(errors) + " errors");
  }
}

/*--------------------------------------------------------------------------*/
/*---------------------- SOLUTION SUPPORT TESTS ----------------------------*/
/*--------------------------------------------------------------------------*/

REGISTER_TEST(solution_basic_functionality) {
  // Test that ScenarioReductionSolver correctly writes solution to Block's y variables
  // and that different algorithms produce valid results
  
  // Test 1: Basic Dupacova algorithm functionality
  {
    auto* block = new CapacitatedFacilityLocationBlock();
    
    // Set up a simple test instance - need square matrix for scenario reduction
    int nf = 2;  // facilities (must equal customers for square distance matrix)
    int nc = 2;  // customers (scenarios)
    int k = 1;  // number of scenarios to select
    
    CapacitatedFacilityLocationBlock::DVector caps(nf);
    caps[0] = 0.5;
    caps[1] = 0.5;
    
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
    
    block->load(nf, nc, caps, fcosts, dems, tcosts, false, k);
    
    // Create and configure solver
    ScenarioReductionSolver solver;
    solver.set_Block(block);
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 1); // Dupacova
    
    // Compute solution
    int status = solver.compute();
    if (status != Solver::kOK) {
      delete block;
      throw std::runtime_error("Solver failed to compute");
    }
    
    // Get the solver's solution
    const auto& reduced = solver.get_reduced_atoms();
    
    // Verify the solver wrote solution to the Block's y variables
    solver.get_var_solution();
    
    // Check that y variables have the correct values
    for (int i = 0; i < nc; ++i) {
      auto* y_var = block->get_y(i);
      if (!y_var) {
        delete block;
        throw std::runtime_error("y variable not found");
      }
      double val = y_var->get_value();
      bool is_selected = (val > 0.5);  // y should be 0 or 1
      if (is_selected != reduced[i]) {
        delete block;
        throw std::runtime_error("y variable value doesn't match solver's solution");
      }
    }
    
    // Verify the solver computed successfully
    int selected_count = std::count(reduced.begin(), reduced.end(), true);
    if (selected_count != k) {
      delete block;
      throw std::runtime_error("Wrong number of scenarios selected");
    }
    
    delete block;
  }
  
  // Test 2: Baseline algorithm selects scenarios with highest weights
  {
    auto* block = create_test_block(2);  // Select 2 out of 4 scenarios
    ScenarioReductionSolver solver;
    solver.set_Block(block);
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Baseline
    
    int result = solver.compute();
    if (result != Solver::kOK) {
      delete block;
      throw std::runtime_error("Baseline algorithm compute() failed");
    }
    
    // Get the selected scenarios
    const auto& selected = solver.get_reduced_atoms();
    
    // Get the demands (weights)
    const auto& demands = block->get_Demands();
    
    // Verify that selected scenarios have the highest weights
    std::vector<double> selected_weights, unselected_weights;
    for (CapacitatedFacilityLocationBlock::Index i = 0; i < block->get_NFacilities(); ++i) {
      if (selected[i]) {
        selected_weights.push_back(demands[i]);
      } else {
        unselected_weights.push_back(demands[i]);
      }
    }
    
    // The demands were set as 50, 100, 150, 200
    // So for k=2, we should select indices 2 and 3 (with weights 150 and 200)
    if (selected_weights.size() != 2) {
      delete block;
      throw std::runtime_error("Should select exactly 2 scenarios");
    }
    
    // All selected weights should be >= all unselected weights
    double min_selected = *std::min_element(selected_weights.begin(), selected_weights.end());
    double max_unselected = unselected_weights.empty() ? 0 : 
                            *std::max_element(unselected_weights.begin(), unselected_weights.end());
    
    if (min_selected < max_unselected) {
      delete block;
      throw std::runtime_error("Baseline didn't select scenarios with highest weights");
    }
    
    // Check that the solution can be written
    solver.get_var_solution();
    
    delete block;
  }
  
  // Test 3: Compare Baseline vs Dupacova algorithms
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
      throw std::runtime_error(
        "Baseline objective (" + std::to_string(baseline_obj) + 
        ") should not be better than Dupacova (" + std::to_string(dupacova_obj) + ")"
      );
    }
    
    delete block;
  }
}

REGISTER_TEST(solution_independence) {
  // Test that multiple solvers can work independently on the same block
  auto* block = create_test_block(2);
  
  ScenarioReductionSolver solver1;
  solver1.set_Block(block);
  solver1.compute();
  
  // Save first solver's solution
  std::vector<bool> solution1 = solver1.get_reduced_atoms();
  double obj1 = solver1.get_var_value();
  
  // Create second solver with different algorithm
  ScenarioReductionSolver solver2;
  solver2.set_Block(block);
  solver2.set_par(ScenarioReductionSolver::intAlgorithm, 0); // Baseline
  solver2.compute();
  
  // Verify both solvers maintain their own solutions
  const auto& solution1_check = solver1.get_reduced_atoms();
  for (size_t i = 0; i < solution1.size(); ++i) {
    if (solution1[i] != solution1_check[i]) {
      delete block;
      throw std::runtime_error("First solver's solution was corrupted");
    }
  }
  
  delete block;
}

REGISTER_TEST(solution_persistence) {
  // Test that solver state can be saved and restored through parameters
  auto* block = create_test_block(2);
  
  ScenarioReductionSolver solver;
  solver.set_Block(block);
  solver.compute();
  
  // Save solver state
  int algorithm = solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
  double rho = solver.get_dbl_par(ScenarioReductionSolver::dblRho);
  double ell = solver.get_dbl_par(ScenarioReductionSolver::dblEll);
  std::vector<bool> solution = solver.get_reduced_atoms();
  double obj_value = solver.get_var_value();
  
  // Create new solver and restore state
  ScenarioReductionSolver solver2;
  solver2.set_Block(block);
  solver2.set_par(ScenarioReductionSolver::intAlgorithm, algorithm);
  solver2.set_par(ScenarioReductionSolver::dblRho, rho);
  solver2.set_par(ScenarioReductionSolver::dblEll, ell);
  solver2.compute();
  
  // Verify same solution is obtained
  const auto& solution2 = solver2.get_reduced_atoms();
  if (solution.size() != solution2.size()) {
    delete block;
    throw std::runtime_error("Solution sizes don't match");
  }
  
  // With same parameters, should get same solution
  for (size_t i = 0; i < solution.size(); ++i) {
    if (solution[i] != solution2[i]) {
      delete block;
      throw std::runtime_error("Solutions don't match after restoring parameters");
    }
  }
  
  delete block;
}

REGISTER_TEST(solution_warm_start) {
  // Test warm start functionality through vintWarmstartIndices parameter
  auto* block = create_test_block(4);
  
  // First solve to get a reference solution
  ScenarioReductionSolver solver1;
  solver1.set_Block(block);
  solver1.set_par(ScenarioReductionSolver::intAlgorithm, 1); // Dupacova
  solver1.compute();
  
  // Get the solution indices
  const auto& solution1 = solver1.get_reduced_atoms();
  std::vector<int> warm_indices;
  for (size_t i = 0; i < solution1.size(); ++i) {
    if (solution1[i]) {
      warm_indices.push_back(static_cast<int>(i));
    }
  }
  
  // Create second solver with warm start
  ScenarioReductionSolver solver2;
  solver2.set_Block(block);
  solver2.set_par(ScenarioReductionSolver::intAlgorithm, 2); // Local search
  solver2.set_par(ScenarioReductionSolver::intUseWarmstart, 1);
  solver2.set_par(ScenarioReductionSolver::vintWarmstartIndices, std::move(warm_indices));
  solver2.compute();
  
  // Verify warm start was used (local search with good warm start should give good result)
  double obj1 = solver1.get_var_value();
  double obj2 = solver2.get_var_value();
  
  if (obj2 > obj1 * 1.1) { // Local search shouldn't be much worse than Dupacova with warm start
    delete block;
    throw std::runtime_error("Warm start doesn't seem to be working properly");
  }
  
  delete block;
}

REGISTER_TEST(has_var_solution_test) {
  // Test has_var_solution() method
  
  // Test 1: has_var_solution() returns false before compute()
  auto* block1 = create_test_block(2);
  ScenarioReductionSolver solver1;
  solver1.set_Block(block1);
  
  if (solver1.has_var_solution()) {
    delete block1;
    throw std::runtime_error("has_var_solution() should return false before compute()");
  }
  
  // Test 2: has_var_solution() returns true after successful compute()
  int status = solver1.compute();
  if (status != Solver::kOK) {
    delete block1;
    throw std::runtime_error("Compute failed");
  }
  
  if (!solver1.has_var_solution()) {
    delete block1;
    throw std::runtime_error("has_var_solution() should return true after successful compute()");
  }
  
  delete block1;
  
  // Test 3: has_var_solution() returns false when no Block is set
  ScenarioReductionSolver solver2;
  if (solver2.has_var_solution()) {
    throw std::runtime_error("has_var_solution() should return false when no Block is set");
  }
  
  // Test 4: Test with different algorithms
  auto* block2 = create_test_block(2);
  ScenarioReductionSolver solver3;
  solver3.set_Block(block2);
  
  // Test Dupacova (default)
  solver3.compute();
  if (!solver3.has_var_solution()) {
    delete block2;
    throw std::runtime_error("has_var_solution() should return true after Dupacova");
  }
  
  // Test BestFit
  solver3.set_par(ScenarioReductionSolver::intAlgorithm, 2);
  solver3.compute();
  if (!solver3.has_var_solution()) {
    delete block2;
    throw std::runtime_error("has_var_solution() should return true after BestFit");
  }
  
  // Test FirstFit
  solver3.set_par(ScenarioReductionSolver::intAlgorithm, 3);
  solver3.compute();
  if (!solver3.has_var_solution()) {
    delete block2;
    throw std::runtime_error("has_var_solution() should return true after FirstFit");
  }
  
  delete block2;
  
  // Test 5: Edge case - k = n (select all scenarios)
  auto* block3 = create_test_block(4);  // 4 facilities, select all 4
  ScenarioReductionSolver solver4;
  solver4.set_Block(block3);
  solver4.compute();
  
  if (!solver4.has_var_solution()) {
    delete block3;
    throw std::runtime_error("has_var_solution() should return true when k = n");
  }
  
  delete block3;
  
  std::cout << "✓ has_var_solution() tests passed\n";
}

REGISTER_TEST(get_var_solution_test) {
  // Test get_var_solution() method
  
  // Test 1: Error when no Block is set
  {
    ScenarioReductionSolver solver;
    try {
      solver.get_var_solution();
      throw std::runtime_error("Should throw when no Block is set");
    } catch (const std::logic_error& e) {
      // Expected
      if (std::string(e.what()).find("no Block set") == std::string::npos) {
        throw std::runtime_error("Wrong error message");
      }
    }
  }
  
  // Test 2: Error when no solution is available
  {
    auto* block = create_test_block(2);
    ScenarioReductionSolver solver;
    solver.set_Block(block);
    
    try {
      solver.get_var_solution();
      throw std::runtime_error("Should throw when no solution available");
    } catch (const std::logic_error& e) {
      // Expected
      if (std::string(e.what()).find("no solution available") == std::string::npos) {
        delete block;
        throw std::runtime_error("Wrong error message");
      }
    }
    delete block;
  }
  
  // Test 3: Variables are now automatically generated by set_Block
  // This test is no longer valid since we follow the MILPSolver pattern
  // where set_Block automatically generates variables
  
  // Test 4: Successful solution writing
  {
    auto* block = create_test_block(2);
    
    ScenarioReductionSolver solver;
    solver.set_Block(block);
    solver.compute();
    
    // Get the solution before writing
    const auto& reduced_atoms = solver.get_reduced_atoms();
    
    // Write solution to block
    solver.get_var_solution();
    
    // Verify values were written correctly
    for (CapacitatedFacilityLocationBlock::Index i = 0; i < block->get_NFacilities(); ++i) {
      auto y_var = block->get_y(i);
      if (!y_var) {
        delete block;
        throw std::runtime_error("Variable y[" + std::to_string(i) + "] is null");
      }
      
      double expected = reduced_atoms[i] ? 1.0 : 0.0;
      double actual = y_var->get_value();
      
      if (!approx_equal(actual, expected)) {
        delete block;
        throw std::runtime_error(
          "Variable y[" + std::to_string(i) + "] has wrong value: " +
          std::to_string(actual) + " != " + std::to_string(expected)
        );
      }
    }
    
    delete block;
  }
  
  // Test 5: Test with all algorithms
  {
    for (int algo = 0; algo <= 3; ++algo) {  // All algorithms including Baseline
      auto* block = create_test_block(3);
      
      ScenarioReductionSolver solver;
      solver.set_Block(block);
      solver.set_par(ScenarioReductionSolver::intAlgorithm, algo);
      solver.compute();
      solver.get_var_solution();
      
      // Just verify it doesn't crash and some variables are set
      int count_selected = 0;
      for (CapacitatedFacilityLocationBlock::Index i = 0; i < block->get_NFacilities(); ++i) {
        auto y_var = block->get_y(i);
        if (y_var && approx_equal(y_var->get_value(), 1.0)) {
          count_selected++;
        }
      }
      
      if (count_selected != 3) {  // Should select exactly k=3 facilities
        delete block;
        throw std::runtime_error(
          "Algorithm " + std::to_string(algo) + " selected " + 
          std::to_string(count_selected) + " facilities instead of 3"
        );
      }
      
      delete block;
    }
  }
  
  std::cout << "✓ get_var_solution() tests passed\n";
}

REGISTER_TEST(refresh_cached_data_error_handling) {
  // Test error handling in refresh_cached_data() method
  
  // Test 1: Negative k value
  {
    auto* block = create_test_block(2);
    ScenarioReductionSolver solver;
    
    try {
      solver.set_Block(block);
      solver.refresh_cached_data(-1);
      delete block;
      throw std::runtime_error("Should throw for negative k");
    } catch (const std::invalid_argument& e) {
      // Expected
      if (std::string(e.what()).find("k must be non-negative") == std::string::npos) {
        delete block;
        throw std::runtime_error("Wrong error message for negative k");
      }
    }
    delete block;
  }
  
  // Test 2: k exceeds number of scenarios
  {
    auto* block = create_test_block(2);  // Block has 4 scenarios
    ScenarioReductionSolver solver;
    
    try {
      solver.set_Block(block);
      solver.refresh_cached_data(10);  // Try to select 10 out of 4
      delete block;
      throw std::runtime_error("Should throw when k exceeds number of scenarios");
    } catch (const std::invalid_argument& e) {
      // Expected
      if (std::string(e.what()).find("k cannot exceed number of scenarios") == std::string::npos) {
        delete block;
        throw std::runtime_error("Wrong error message for k > n");
      }
    }
    delete block;
  }
  
  // Test 3: Non-square matrix (customers != facilities)
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
    
    block->load(nf, nc, caps, fcosts, dems, tcosts, false, k);
    
    ScenarioReductionSolver solver;
    
    try {
      solver.set_Block(block);  // This should throw in refresh_cached_data
      delete block;
      throw std::runtime_error("Should throw for non-square matrix");
    } catch (const std::logic_error& e) {
      // Expected
      if (std::string(e.what()).find("number of customers must equal number of facilities") == std::string::npos) {
        delete block;
        throw std::runtime_error("Wrong error message for non-square matrix");
      }
    }
    delete block;
  }
  
  // Test 4: Valid cases - ensure no exceptions
  {
    auto* block = create_test_block(2);
    ScenarioReductionSolver solver;
    
    // Should not throw
    solver.set_Block(block);
    solver.refresh_cached_data(0);  // k = 0 is valid
    solver.refresh_cached_data(2);  // k = 2 is valid
    solver.refresh_cached_data(4);  // k = n is valid
    
    delete block;
  }
  
  std::cout << "✓ refresh_cached_data() error handling tests passed\n";
}

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main(int argc, char* argv[]) {
  std::cout << "ScenarioReductionSolver Test Suite\n";
  std::cout << "==================================\n\n";
  
  if (argc > 1) {
    std::string arg(argv[1]);
    
    if (arg == "--list" || arg == "-l") {
      TestRegistry::instance().list_tests();
      return 0;
    }
    
    if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " [options] [test_name]\n\n";
      std::cout << "Options:\n";
      std::cout << "  --list, -l     List all available tests\n";
      std::cout << "  --help, -h     Show this help message\n";
      std::cout << "  test_name      Run only the specified test\n\n";
      std::cout << "If no arguments are provided, all tests are run.\n";
      return 0;
    }
    
    // Run specific test
    TestRegistry::instance().run_test(arg);
  } else {
    // Run all tests
    TestRegistry::instance().run_all();
  }
  
  return 0;
}

/*--------------------------------------------------------------------------*/
/*------------------------- End test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/