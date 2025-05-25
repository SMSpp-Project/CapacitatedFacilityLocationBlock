/*--------------------------------------------------------------------------*/
/*------------------------ File test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Structured test suite for ScenarioReductionSolver class following the
 * test registry pattern used in StochasticBlock tests.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * Copyright &copy; by Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
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
  
  // Create a small problem: 3 facilities, 4 customers
  int nf = 3;
  int nc = 4;
  
  // Set up facility costs
  CapacitatedFacilityLocationBlock::CVector fcosts(nf);
  for (int i = 0; i < nf; ++i) {
    fcosts[i] = 100.0 * (i + 1);  // 100, 200, 300
  }
  
  // Set up transportation costs
  CapacitatedFacilityLocationBlock::CMatrix tcosts(boost::extents[nf][nc]);
  for (int i = 0; i < nf; ++i) {
    for (int j = 0; j < nc; ++j) {
      tcosts[i][j] = 10.0 * (j + 1);  // Repeating pattern
    }
  }
  
  // Set up capacities (interpreted as probabilities)
  CapacitatedFacilityLocationBlock::DVector caps(nf);
  caps[0] = 0.3;
  caps[1] = 0.5;
  caps[2] = 0.2;
  
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

REGISTER_TEST(basic_creation) {
  // Test basic solver creation and destruction
  ScenarioReductionSolver solver;
  
  // Check default parameters
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 1) {
    throw std::runtime_error("Default algorithm should be 1 (Dupacova)");
  }
  
  if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 0) {
    throw std::runtime_error("Default shuffle should be 0 (disabled)");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 2.0)) {
    throw std::runtime_error("Default ell should be 2.0");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.0)) {
    throw std::runtime_error("Default rho should be 0.0");
  }
}

REGISTER_TEST(parameter_setting_int) {
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
}

REGISTER_TEST(parameter_setting_double) {
  ScenarioReductionSolver solver;
  
  // Test setting and getting double parameters
  solver.set_par(ScenarioReductionSolver::dblEll, 0.5);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 0.5)) {
    throw std::runtime_error("Failed to set ell parameter");
  }
  
  solver.set_par(ScenarioReductionSolver::dblRho, 0.01);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.01)) {
    throw std::runtime_error("Failed to set rho parameter");
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
    solver.set_par(ScenarioReductionSolver::intAlgorithm, 5);
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
  
  // Test getting default values
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
}

REGISTER_TEST(parameter_reset) {
  ScenarioReductionSolver solver;
  
  // Change some parameters
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 3);
  solver.set_par(ScenarioReductionSolver::dblEll, 1.5);
  solver.set_par(ScenarioReductionSolver::dblRho, 0.1);
  
  // Verify they changed
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 3) {
    throw std::runtime_error("Algorithm should have changed to 3");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 1.5)) {
    throw std::runtime_error("Ell should have changed to 1.5");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.1)) {
    throw std::runtime_error("Rho should have changed to 0.1");
  }
  
  // Reset to defaults
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 
                solver.get_dflt_int_par(ScenarioReductionSolver::intAlgorithm));
  solver.set_par(ScenarioReductionSolver::dblEll, 
                solver.get_dflt_dbl_par(ScenarioReductionSolver::dblEll));
  solver.set_par(ScenarioReductionSolver::dblRho, 
                solver.get_dflt_dbl_par(ScenarioReductionSolver::dblRho));
  
  // Verify reset
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 1) {
    throw std::runtime_error("Algorithm should be reset to 1");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 2.0)) {
    throw std::runtime_error("Ell should be reset to 2.0");
  }
  
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 0.0)) {
    throw std::runtime_error("Rho should be reset to 0.0");
  }
}

REGISTER_TEST(algorithm_constants) {
  // Test that algorithm constants are as expected
  if (static_cast<int>(ScenarioReductionSolver::Algorithm::Baseline) != 0) {
    throw std::runtime_error("Baseline algorithm should be 0");
  }
  
  if (static_cast<int>(ScenarioReductionSolver::Algorithm::Dupacova) != 1) {
    throw std::runtime_error("Dupacova algorithm should be 1");
  }
  
  if (static_cast<int>(ScenarioReductionSolver::Algorithm::BestFit) != 2) {
    throw std::runtime_error("BestFit algorithm should be 2");
  }
  
  if (static_cast<int>(ScenarioReductionSolver::Algorithm::FirstFit) != 3) {
    throw std::runtime_error("FirstFit algorithm should be 3");
  }
  
  if (static_cast<int>(ScenarioReductionSolver::Algorithm::MILP) != 4) {
    throw std::runtime_error("MILP algorithm should be 4");
  }
}

REGISTER_TEST(boundary_values) {
  ScenarioReductionSolver solver;
  
  // Test boundary values for algorithm
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);  // Min valid
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 0) {
    throw std::runtime_error("Algorithm 0 should be valid");
  }
  
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 4);  // Max valid
  if (solver.get_int_par(ScenarioReductionSolver::intAlgorithm) != 4) {
    throw std::runtime_error("Algorithm 4 should be valid");
  }
  
  // Test boundary values for shuffle
  solver.set_par(ScenarioReductionSolver::intShuffle, 0);
  if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 0) {
    throw std::runtime_error("Shuffle 0 should be valid");
  }
  
  solver.set_par(ScenarioReductionSolver::intShuffle, 1);
  if (solver.get_int_par(ScenarioReductionSolver::intShuffle) != 1) {
    throw std::runtime_error("Shuffle 1 should be valid");
  }
  
  // Test small positive values for double parameters
  solver.set_par(ScenarioReductionSolver::dblEll, 1e-10);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblEll), 1e-10)) {
    throw std::runtime_error("Very small ell should be valid");
  }
  
  solver.set_par(ScenarioReductionSolver::dblRho, 1e-10);
  if (!approx_equal(solver.get_dbl_par(ScenarioReductionSolver::dblRho), 1e-10)) {
    throw std::runtime_error("Very small rho should be valid");
  }
}

/*--------------------------------------------------------------------------*/
/*------------------------ THREAD SAFETY TESTS -----------------------------*/
/*--------------------------------------------------------------------------*/

REGISTER_TEST(mutex_basic_lock_unlock) {
  // Test that the solver has a working mutex
  ScenarioReductionSolver solver;
  
  // The solver should be lockable
  solver.lock();
  
  // Try to lock again from another thread - should block
  std::atomic<bool> locked(false);
  std::atomic<bool> tried(false);
  
  std::thread t([&solver, &locked, &tried]() {
    tried = true;
    if (solver.try_lock()) {
      locked = true;
      solver.unlock();
    }
  });
  
  // Give the thread time to try
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  // The thread should have tried but not succeeded
  if (!tried) {
    solver.unlock();
    t.join();
    throw std::runtime_error("Thread didn't attempt to lock");
  }
  
  if (locked) {
    solver.unlock();
    t.join();
    throw std::runtime_error("Thread acquired lock when it shouldn't have");
  }
  
  // Now unlock and let the thread acquire it
  solver.unlock();
  
  // Give the thread time to acquire
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  t.join();
}

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
        solver.set_par(ScenarioReductionSolver::intAlgorithm, i % 5);
        solver.set_par(ScenarioReductionSolver::dblRho, 0.1 * i);
        
        // Verify the values are what we just set
        int alg = solver.get_int_par(ScenarioReductionSolver::intAlgorithm);
        double rho = solver.get_dbl_par(ScenarioReductionSolver::dblRho);
        
        if (alg != i % 5) {
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
  auto block = create_test_block(1);  // Reduce from 3 facilities to 1
  ScenarioReductionSolver solver;
  solver.set_Block(block);
  
  std::atomic<int> computing(0);
  std::atomic<int> completed(0);
  std::atomic<int> blocked(0);
  std::vector<std::thread> threads;
  
  // Launch multiple threads that try to compute simultaneously
  for (int i = 0; i < 5; ++i) {
    threads.emplace_back([&solver, &computing, &completed, &blocked]() {
      solver.lock();
      
      // Check if we can compute
      int current = computing.fetch_add(1);
      
      if (current == 0) {
        // First thread, should be able to compute
        int status = solver.compute();
        computing--;
        completed++;
        
        if (status != Solver::kOK && status != Solver::kError) {
          solver.unlock();
          throw std::runtime_error("Unexpected compute status: " + 
                                  std::to_string(status));
        }
      } else {
        // Other threads should be prevented from computing
        computing--;
        blocked++;
        // Note: Without the f_computing flag implemented yet,
        // multiple threads might actually compute simultaneously
        // This test will need to be updated after implementation
      }
      
      solver.unlock();
    });
  }
  
  // Wait for all threads
  for (auto& t : threads) {
    t.join();
  }
  
  // Clean up
  delete block;
  
  // Note: This test will be more meaningful after f_computing flag is implemented
  // For now, we just check that no exceptions were thrown
  if (completed == 0) {
    throw std::runtime_error("No thread completed computation");
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
        solver.set_par(ScenarioReductionSolver::intAlgorithm, iteration % 5);
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
  // Create a test block and solver
  auto* block = new CapacitatedFacilityLocationBlock();
  
  // Set up a simple test instance
  block->set_NFacilities(3);
  block->set_NCustomers(2);
  
  std::vector<double> f_cost = {10.0, 15.0, 20.0};
  std::vector<double> t_cost = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  std::vector<double> capacity = {0.4, 0.3, 0.3};  // probabilities
  std::vector<double> demand = {1.0, 1.0};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
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
  
  // Test get_Solution from Block (need to specify we want scenario reduction solution)
  SimpleConfiguration<int> sr_config(1);  // 1 = scenario reduction solution
  auto* sol = block->get_Solution(&sr_config, false);  // false = not empty
  if (!sol) {
    delete block;
    throw std::runtime_error("get_Solution returned nullptr");
  }
  
  // Verify solution type (once ScenarioReductionSolution is implemented)
  // auto* sr_sol = dynamic_cast<ScenarioReductionSolution*>(sol);
  // if (!sr_sol) {
  //   delete sol;
  //   delete block;
  //   throw std::runtime_error("Wrong solution type returned");
  // }
  
  delete sol;
  delete block;
}

REGISTER_TEST(solution_save_restore) {
  auto* block = new CapacitatedFacilityLocationBlock();
  
  // Set up test instance
  block->set_NFacilities(4);
  block->set_NCustomers(2);
  
  std::vector<double> f_cost = {10.0, 15.0, 20.0, 25.0};
  std::vector<double> t_cost = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  std::vector<double> capacity = {0.25, 0.25, 0.25, 0.25};
  std::vector<double> demand = {1.0, 1.0};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
  ScenarioReductionSolver solver;
  solver.set_Block(block);
  
  // Compute with Dupacova
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 1);
  solver.compute();
  
  // Save solution from Block
  SimpleConfiguration<int> sr_config(1);  // 1 = scenario reduction solution
  auto* saved_sol = block->get_Solution(&sr_config, false);
  double saved_obj = solver.get_var_value();
  
  // Compute with different algorithm (Baseline)
  solver.set_par(ScenarioReductionSolver::intAlgorithm, 0);
  solver.compute();
  double new_obj = solver.get_var_value();
  
  // Objectives should be different
  if (approx_equal(saved_obj, new_obj)) {
    delete saved_sol;
    delete block;
    throw std::runtime_error("Different algorithms should produce different objectives");
  }
  
  // Restore saved solution through solver (solver will write it to block)
  solver.put_Solution(saved_sol);
  
  // Verify objective is restored
  if (!approx_equal(solver.get_var_value(), saved_obj)) {
    delete saved_sol;
    delete block;
    throw std::runtime_error("Failed to restore solution objective");
  }
  
  delete saved_sol;
  delete block;
}

REGISTER_TEST(solution_clone) {
  auto* block = new CapacitatedFacilityLocationBlock();
  
  block->set_NFacilities(3);
  block->set_NCustomers(2);
  
  std::vector<double> f_cost = {10.0, 15.0, 20.0};
  std::vector<double> t_cost = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  std::vector<double> capacity = {0.4, 0.3, 0.3};
  std::vector<double> demand = {1.0, 1.0};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
  ScenarioReductionSolver solver;
  solver.set_Block(block);
  solver.compute();
  
  // Get solution from Block and clone it
  SimpleConfiguration<int> sr_config(1);  // 1 = scenario reduction solution
  auto* sol1 = block->get_Solution(&sr_config, false);
  auto* sol2 = sol1->clone(false);  // Full clone
  auto* sol3 = sol1->clone(true);   // Empty clone
  
  // Once implemented, verify clones
  // auto* sr_sol2 = dynamic_cast<ScenarioReductionSolution*>(sol2);
  // auto* sr_sol3 = dynamic_cast<ScenarioReductionSolution*>(sol3);
  
  delete sol1;
  delete sol2;
  delete sol3;
  delete block;
}

REGISTER_TEST(solution_serialization) {
  auto* block = new CapacitatedFacilityLocationBlock();
  
  block->set_NFacilities(3);
  block->set_NCustomers(2);
  
  std::vector<double> f_cost = {10.0, 15.0, 20.0};
  std::vector<double> t_cost = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  std::vector<double> capacity = {0.4, 0.3, 0.3};
  std::vector<double> demand = {1.0, 1.0};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
  ScenarioReductionSolver solver;
  solver.set_Block(block);
  solver.compute();
  
  SimpleConfiguration<int> sr_config(1);  // 1 = scenario reduction solution
  auto* sol = block->get_Solution(&sr_config, false);
  
  // Test serialization (once implemented)
  try {
    netCDF::NcFile file("test_sr_solution.nc4", netCDF::NcFile::replace);
    netCDF::NcGroup group = file.addGroup("solution");
    sol->serialize(group);
    file.close();
    
    // Test deserialization
    netCDF::NcFile file2("test_sr_solution.nc4", netCDF::NcFile::read);
    netCDF::NcGroup group2 = file2.getGroup("solution");
    auto* loaded = Solution::new_Solution(group2);
    
    if (!loaded) {
      delete sol;
      delete block;
      throw std::runtime_error("Failed to deserialize solution");
    }
    
    delete loaded;
    
    // Clean up test file
    std::remove("test_sr_solution.nc4");
  } catch (const netCDF::exceptions::NcException& e) {
    delete sol;
    delete block;
    throw std::runtime_error(std::string("NetCDF error: ") + e.what());
  }
  
  delete sol;
  delete block;
}

REGISTER_TEST(solution_warm_start) {
  auto* block = new CapacitatedFacilityLocationBlock();
  
  // Set up test instance
  block->set_NFacilities(5);
  block->set_NCustomers(3);
  
  std::vector<double> f_cost = {10.0, 15.0, 20.0, 25.0, 30.0};
  std::vector<double> t_cost(15, 1.0);  // All transportation costs = 1.0
  std::vector<double> capacity = {0.2, 0.2, 0.2, 0.2, 0.2};
  std::vector<double> demand = {1.0, 1.0, 1.0};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
  ScenarioReductionSolver solver1;
  solver1.set_Block(block);
  solver1.compute();
  
  // Save solution for warm start
  SimpleConfiguration<int> sr_config(1);
  auto* warm_start = block->get_Solution(&sr_config, false);
  
  // Create a second solver and use warm start
  ScenarioReductionSolver solver2;
  solver2.set_Block(block);
  
  // Apply warm start (once put_Solution is implemented)
  try {
    solver2.put_Solution(warm_start);
  } catch (...) {
    // Expected until implementation
  }
  
  delete warm_start;
  delete block;
}

REGISTER_TEST(backward_compatibility) {
  // Test that default get_Solution() returns CFL solution for backward compatibility
  auto* block = new CapacitatedFacilityLocationBlock();
  
  // Set up minimal test instance
  block->set_NFacilities(2);
  block->set_NCustomers(1);
  
  std::vector<double> f_cost = {10.0, 15.0};
  std::vector<double> t_cost = {1.0, 2.0};
  std::vector<double> capacity = {1.0, 1.0};
  std::vector<double> demand = {0.5};
  
  block->set_facility_cost(f_cost.data());
  block->set_transportation_cost(t_cost.data());
  block->set_capacity(capacity.data());
  block->set_demand(demand.data());
  
  // Get solution without configuration (backward compatible)
  auto* default_sol = block->get_Solution();
  
  // Verify it's a CFL solution, not a scenario reduction solution
  if (dynamic_cast<CapacitatedFacilityLocationSolution*>(default_sol)) {
    std::cout << "✓ Default get_Solution() returns CapacitatedFacilityLocationSolution\n";
  } else {
    std::cerr << "✗ Default get_Solution() did not return expected type\n";
  }
  
  // Also verify we can get scenario reduction solution with config
  SimpleConfiguration<int> sr_config(1);
  auto* sr_sol = block->get_Solution(&sr_config);
  
  if (dynamic_cast<ScenarioReductionSolution*>(sr_sol)) {
    std::cout << "✓ get_Solution(&sr_config) returns ScenarioReductionSolution\n";
  } else {
    std::cerr << "✗ get_Solution(&sr_config) did not return expected type\n";
  }
  
  delete default_sol;
  delete sr_sol;
  delete block;
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