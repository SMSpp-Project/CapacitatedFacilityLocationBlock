/*--------------------------------------------------------------------------*/
/*----------------------- File test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test file for the ScenarioReductionSolver class.
 * This test evaluates different scenario reduction algorithms:
 * - Dupacova's forward selection algorithm
 * - Local search with BestFit strategy
 * - Local search with FirstFit strategy (with and without shuffling)
 * - Comparison with optimal solution from MILPSolver
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"

#include "ScenarioReductionSolver.h"
#include "CapacitatedFacilityLocationBlock.h"
#include "DiscreteScenarioSet.h"

#include "BlockSolverConfig.h"
#include "Configuration.h"

#include <format> // C++20 string formatting
#include <iostream>
#include <iomanip>  // For std::setw
#include <vector>
#include <boost/multi_array.hpp>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace netCDF;

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

using Exponent = float;
using ScenarioIndex = ScenarioGenerator::ScenarioIndex;
using Scenario = ScenarioGenerator::Scenario;
using ScenarioSize = ScenarioGenerator::ScenarioSize;

/*--------------------------------------------------------------------------*/
/*-------------------------- GLOBAL VARIABLES ------------------------------*/
/*--------------------------------------------------------------------------*/

std::mt19937 rng;

/*--------------------------------------------------------------------------*/
/*---------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * Print the solution from a ScenarioReductionSolver
 * 
 * @param solver Pointer to the ScenarioReductionSolver with a solution
 * @param block Pointer to the CapacitatedFacilityLocationBlock with problem data
 * @param algorithm_name Name of the algorithm used (for display purposes)
 */
void print_solution(ScenarioReductionSolver* solver, CapacitatedFacilityLocationBlock* block, 
                    const std::string& algorithm_name = "algorithm") {
    std::cout << "Solution from " << algorithm_name << ":" << std::endl;
    if (!solver || !block) {
        std::cerr << "Error: Null pointer!" << std::endl;
        return;
    }
    
    // Calculate and display Wasserstein distance
    float ell = solver->get_ell();
    double distance_value = pow(solver->get_var_value(), ell);
    std::cout << std::format("Wasserstein distance power {} between original and reduced distrib: {} ", 
                            ell, distance_value) << std::endl;
    
    // Display the selected scenarios (facilities in the CFL formulation)
    const auto& selected_scenarios = solver->get_reduced_atoms();
    std::cout << "   Reduced scenario set: ";
    for (int i = 0; i < block->get_NFacilities(); i++) {
        if (selected_scenarios[i]) {
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
}

/**
 * Print solution based on the MILP solver's output
 * 
 * @param block Pointer to the CapacitatedFacilityLocationBlock containing solution
 */
void print_cfl_solution(CapacitatedFacilityLocationBlock* block) {
    // Check if block is valid
    if (!block) {
        std::cerr << "Error: Null block pointer" << std::endl;
        return;
    }
    
    // Get registered solvers
    auto solvers = block->get_registered_solvers();
    if (solvers.empty()) {
        std::cerr << "Error: No solver registered to block" << std::endl;
        return;
    }
    
    // First try to find a MILP solver with solution
    Solver* solver = nullptr;
    for (auto s : solvers) {
        if (dynamic_cast<ScenarioReductionSolver*>(s) == nullptr) {
            solver = s;  // This is the MILP solver
            if (solver->has_var_solution()) {
                break;  // Found MILP solver with solution
            }
        }
    }
    
    // If no MILP solver with solution, try any solver
    if (!solver || !solver->has_var_solution()) {
        for (auto s : solvers) {
            if (s->has_var_solution()) {
                solver = s;
                break;
            }
        }
    }
    
    // Check if we have a solution
    if (!solver || !solver->has_var_solution()) {
        std::cerr << "Error: No solution available from any solver" << std::endl;
        return;
    }
    
    std::cout << "Using solver: " << solver->classname() << std::endl;
    
    // Get the solution
    solver->get_var_solution();
    
    // Get facility decisions (which scenarios are selected)
    CapacitatedFacilityLocationBlock::CntSolution y(block->get_NFacilities());
    block->get_facility_solution(y.begin());
    
    // Print objective values
    double raw_value = solver->get_lb();
    double normalized_value = raw_value / block->get_NFacilities();
    std::cout << "Raw objective value: " << raw_value << std::endl;
    std::cout << "Normalized objective value: " << normalized_value << std::endl;
    std::cout << "This value should match the Wasserstein distance²" << std::endl;

    // Print selected scenarios
    std::cout << "Selected scenarios: ";
    for (size_t i = 0; i < y.size(); i++) {
        if (y[i] > 0) {  
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
}

/*--------------------------------------------------------------------------*/
/**
 * Convert a scenario reduction problem into a capacitated facility location problem
 * 
 * This function formulates the scenario reduction problem as a CFL problem by:
 * - Treating each scenario as both a potential facility and a customer
 * - Setting customer demands to scenario probability weights
 * - Setting all facility capacities to 1.0
 * - Setting all facility fixed costs to 0.0
 * - Computing transportation costs as ell-power of p-norm distances between scenarios
 *
 * @param block The CapacitatedFacilityLocationBlock to load data into
 * @param scenarios A DiscreteScenarioPool holding the scenarios
 * @param weights Probability weights associated with each scenario
 * @param k The maximum number of scenarios to select in the reduced set
 * @param ell The ell-Wasserstein distance power to be minimized (default: 2)
 * @param p The p-norm to use for computing distances (default: 2)
 * @throws std::logic_error if block is nullptr
 * @throws std::invalid_argument if input parameters are invalid
 */
void load_scenario_reduction_problem(
    CapacitatedFacilityLocationBlock* block,
    const DiscreteScenarioSet::DiscreteScenarioPool& scenarios,
    const std::vector<double>& weights,
    CapacitatedFacilityLocationBlock::Index k,
    double ell = 2.0,
    double p = 2.0)
{
    // Validate input parameters
    if (!block) {
        throw std::logic_error("Invalid CapacitatedFacilityLocationBlock pointer provided");
    }
    
    ScenarioIndex n_scenarios = scenarios.shape()[0];
    ScenarioSize scenario_size = scenarios.shape()[1];
    
    if (weights.size() != n_scenarios) {
        throw std::invalid_argument("Number of weights must match number of scenarios");
    }
    if (k > n_scenarios) {
        throw std::invalid_argument("k cannot exceed the number of scenarios");
    }
    if (ell <= 0 || p <= 0) {
        throw std::invalid_argument("ell and p must be positive");
    }
    if (std::any_of(weights.begin(), weights.end(), [](double w) { return w < 0; })) {
        throw std::invalid_argument("All weights must be non-negative");
    }
    
    // Create CFL problem parameters
    CapacitatedFacilityLocationBlock::DVector capacities(n_scenarios);
    CapacitatedFacilityLocationBlock::CVector fixed_costs(n_scenarios);
    CapacitatedFacilityLocationBlock::DVector demands(n_scenarios);
    CapacitatedFacilityLocationBlock::CMatrix transport_costs(boost::extents[n_scenarios][n_scenarios]);
    
    // Map scenarios to Eigen matrix for easier distance calculations
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
        all_scenarios(scenarios.data(), n_scenarios, scenario_size);
    
    // Set up the basic CFL parameters
    for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
        capacities[i] = 1.0;      // Each facility can serve at most one unit
        fixed_costs[i] = 0.0;     // No fixed cost in scenario reduction
        demands[i] = weights[i];  // Demand equals probability weight

        std::cout << "Scenario " << i << ": " << all_scenarios.row(i) << std::endl;
    }
    
    // Compute the distance/transportation cost matrix
    for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
        for (ScenarioIndex j = 0; j < n_scenarios; ++j) {
            if (i == j) {
                transport_costs[i][j] = 0.0;
                continue;
            }
            
            // Calculate p-norm distance between scenarios
            Eigen::VectorXd diff = all_scenarios.row(i) - all_scenarios.row(j);
            double norm;
            
            if (p == 2) {
                norm = diff.norm();  // Euclidean norm for p=2 (more efficient)
            } else {
                norm = std::pow((diff.array().abs().pow(p)).sum(), 1.0/p);
            }
            
            // Transportation cost = ell-power of the norm
            transport_costs[i][j] = std::pow(norm, ell);
        }
    }

    // Load the CFL problem into the provided block
    block->load(
        n_scenarios,       // Number of facilities
        n_scenarios,       // Number of customers
        std::move(capacities),
        std::move(fixed_costs),
        std::move(demands),
        std::move(transport_costs),
        false,             // Not a balanced problem
        k                  // Maximum number of facilities
    );
}

/*--------------------------------------------------------------------------*/

/**
 * Generate a set of distinct scenarios with uniform weights
 *
 * This function creates random scenarios with normally distributed values,
 * ensuring each scenario is unique. The values are truncated to [-20, 20].
 *
 * @param scenarios Reference to a DiscreteScenarioPool to store scenarios
 * @param weights Reference to a vector to store probability weights
 * @param n_scenar Number of scenarios to generate
 * @param scen_size Size of each scenario
 * @param saveToFile Whether to save scenarios and weights to a file
 * @param filename Name of the file to save to (if saveToFile is true)
 */
void generateDistinctScenarios(
    DiscreteScenarioSet::DiscreteScenarioPool& scenarios,
    std::vector<double>& weights,
    size_t n_scenar = 5,
    size_t scen_size = 10,
    bool saveToFile = false,
    const std::string& filename = "scenarios.txt") 
{
    // Create normally distributed values (mean=5.0, stddev=5.0)
    std::normal_distribution<double> normalDist(5.0, 5.0);

    // Set uniform weights
    std::fill(weights.begin(), weights.end(), 1.0 / n_scenar);

    // Generate distinct scenarios
    std::set<std::vector<double>> scenarioSet;
    size_t count = 0;

    while (count < n_scenar) {
        // Generate a candidate scenario
        std::vector<double> scenario(scen_size);
        for (size_t j = 0; j < scen_size; ++j) {
            // Generate truncated normal values between -20 and 20
            double value = normalDist(rng);
            while (value < -20.0 || value > 20.0) {
                value = normalDist(rng);
            }
            scenario[j] = value;
        }

        // Only add if this scenario is distinct from previously generated ones
        if (scenarioSet.find(scenario) == scenarioSet.end()) {
            // Store the scenario in the multi_array
            for (size_t j = 0; j < scen_size; ++j) {
                scenarios[count][j] = scenario[j];
            }
            scenarioSet.insert(scenario);
            count++;
        }
    }

    // Save to file if requested
    if (saveToFile) {
        std::ofstream outFile(filename);
        if (outFile.is_open()) {
            // Write weights first
            outFile << "Weights:\n";
            for (size_t i = 0; i < weights.size(); ++i) {
                outFile << weights[i];
                if (i < weights.size() - 1) {
                    outFile << " ";
                }
            }
            outFile << "\n\n";

            // Write scenarios
            outFile << "Scenarios:\n";
            for (size_t i = 0; i < n_scenar; ++i) {
                outFile << "Scenario " << i + 1 << ": ";
                for (size_t j = 0; j < scen_size; ++j) {
                    outFile << scenarios[i][j];
                    if (j < scen_size - 1) {
                        outFile << " ";
                    }
                }
                outFile << "\n";
            }

            outFile.close();
            std::cout << "Scenarios and weights saved to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file " << filename << " for writing" << std::endl;
        }
    }
}

/**
 * Test a specific scenario reduction algorithm with given parameters
 * 
 * @param solver Pointer to the ScenarioReductionSolver
 * @param block Pointer to the CapacitatedFacilityLocationBlock
 * @param algorithm Algorithm to test
 * @param algorithm_name Name of the algorithm for display
 * @param rho Minimum improvement threshold parameter (for local search)
 * @param shuffle Whether to enable shuffling (for FirstFit)
 */
void test_algorithm(
    ScenarioReductionSolver* solver, 
    CapacitatedFacilityLocationBlock* block, 
    ScenarioReductionSolver::Algorithm algorithm, 
    const std::string& algorithm_name,
    double rho = 0.0, 
    bool shuffle = false) 
{
    std::cout << "\n=== Testing " << algorithm_name << " algorithm ===" << std::endl;
    
    // Configure the solver with the specified algorithm and parameters
    solver->set_algorithm(algorithm);
    solver->set_rho(rho);
    solver->set_shuffle(shuffle);
    
    // Display configuration
    std::cout << "Configuration:" << std::endl;
    std::cout << "  - Algorithm: " << algorithm_name << std::endl;
    std::cout << "  - rho: " << solver->get_rho() << std::endl;
    std::cout << "  - shuffle: " << (solver->get_shuffle() ? "enabled" : "disabled") << std::endl;
    
    // Solve with the selected algorithm
    std::cout << "Solving with " << algorithm_name << " algorithm..." << std::endl;
    int result = solver->compute();
    
    // Display results
    std::cout << "Solver returned: " << result << std::endl;
    if (result == Solver::kOK) {
        print_solution(solver, block, algorithm_name);
    } else {
        std::cout << "Failed to find a solution with " << algorithm_name << " algorithm." << std::endl;
    }
}

/**
 * Compare results from multiple algorithm runs
 * 
 * @param results Vector of tuples containing (algorithm name, distance, selected scenarios)
 */
void compare_results(const std::vector<std::tuple<std::string, double, std::vector<int>>>& results) {
    std::cout << "\n=== Algorithm Comparison ===" << std::endl;
    
    // Print results in a formatted table
    std::cout << std::left;
    std::cout << std::setw(20) << "Algorithm" << std::setw(25) << "Wasserstein Distance²" 
              << "Selected Scenarios" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (const auto& [algo_name, distance, scenarios] : results) {
        std::cout << std::setw(20) << algo_name << std::setw(25) << distance << ": ";
        for (int scenario : scenarios) {
            std::cout << scenario << " ";
        }
        std::cout << std::endl;
    }
    
    // Find the best algorithm based on minimum distance
    auto best_algo = std::min_element(results.begin(), results.end(),
        [](const auto& a, const auto& b) { return std::get<1>(a) < std::get<1>(b); });
    
    std::cout << "\nBest performing algorithm: " << std::get<0>(*best_algo) 
              << " with distance² = " << std::get<1>(*best_algo) << std::endl;
    
    // Theoretical validation: MILP solver should provide the optimal solution
    auto milp_it = std::find_if(results.begin(), results.end(),
        [](const auto& res) { return std::get<0>(res) == "HiGHSMILPSolver"; });
        
    if (milp_it != results.end()) {
        // Allow for small numerical precision differences (e.g., 1e-6)
        const double EPSILON = 1e-6;
        bool milp_is_best = false;
        
        if (std::abs(std::get<1>(*milp_it) - std::get<1>(*best_algo)) < EPSILON) {
            // MILP is within epsilon of the best, consider it optimal
            milp_is_best = true;
        } else if (std::get<1>(*milp_it) <= std::get<1>(*best_algo)) {
            // MILP has lower or equal distance
            milp_is_best = true;
        }
        
        if (!milp_is_best) {
            std::cout << "\nWARNING: The MILP solver did not yield the best solution, which contradicts theory.\n"
                      << "This suggests a possible implementation issue or numerical precision error.\n"
                      << "MILP distance: " << std::get<1>(*milp_it) << ", Best distance: " << std::get<1>(*best_algo) << std::endl;
        } else {
            std::cout << "\nTheoretical check passed: The MILP solver yielded the best solution as expected." << std::endl;
        }
    }
}

/*--------------------------------------------------------------------------*/
/*---------------------------- MAIN -------------------------------------*/
/*--------------------------------------------------------------------------*/

int main() {
    try {
        // Part 1: Setup the scenario reduction problem
        std::cout << "\n=== PART 1: Setting up a scenario reduction problem ===" << std::endl;
        
        // Create a toy scenarios with uniform weights
        const size_t num_scenarios = 5;
        const size_t scenario_size = 10;
        std::cout << std::format("Creating toy scenario set with {} scenarios of size {}...", 
            num_scenarios, scenario_size) << std::endl;

        DiscreteScenarioSet::DiscreteScenarioPool scenarios(
            boost::extents[num_scenarios][scenario_size]);
        std::vector<double> weights(num_scenarios);

        generateDistinctScenarios(scenarios, weights, num_scenarios, scenario_size, true);
        
        // Create a new CFLBlock and load the scenario reduction problem
        const int k = 3;  // Number of scenarios to select
        std::cout << std::format("Loading scenario reduction problem (k={})...", k) << std::endl;
        auto sr_block = new CapacitatedFacilityLocationBlock();
        load_scenario_reduction_problem(sr_block, scenarios, weights, k);

        // Verify the physical variables
        std::cout << "Verifying physical model:" << std::endl;
        std::cout << std::format("  Number of facilities (scenarios): {}", sr_block->get_NFacilities()) << std::endl;
        std::cout << std::format("  Number of customers (scenarios): {}", sr_block->get_NCustomers()) << std::endl;
        std::cout << std::format("  Maximum facilities to select (k): {}", sr_block->get_NMaxFacilities()) << std::endl;
        std::cout << "  All capacities are 1.0: "
                  << std::boolalpha
                  << std::all_of(sr_block->get_Capacities().begin(), sr_block->get_Capacities().end(),
                                [](double val) { return val == 1.0; }) << std::endl;
        
        // Generating and sanity check of the abstract representation in the Block
        std::cout << "Verifying abstract model:" << std::endl;

        auto cfg = Configuration::deserialize( "BPar.txt" );
        if( BlockConfig * bc = dynamic_cast< BlockConfig * >( cfg ) )
            bc->apply( sr_block );
        else {
            std::cerr << "Error: BPar.txt does not contain a BlockConfig" << std::endl;
            return 1;
        }
        sr_block->generate_abstract_variables();
        sr_block->generate_abstract_constraints();
        
        // Check static variables
        auto static_vars = sr_block->get_static_variables();
        std::cout << "  Number of static variable groups: " << static_vars.size() << std::endl;
        if (static_vars.size() == 2) {
            std::cout << "    Expected: 2 groups (y variables and x variables)" << std::endl;
            
            // Check facility variables (y) 
            if (sr_block->get_y(0)) {
                std::cout << "    Found y variables (facility variables): " << std::endl;
            }
            
            // Check transportation variables (x)
            if (sr_block->get_x(0, 0)) {
                std::cout << "    Found x variables (transportation variables): " << std::endl;
            }
        } else {
            std::cerr << "    WARNING: Expected 2 variable groups" << std::endl;
        }
        
        // Check static constraints
        auto static_constr = sr_block->get_static_constraints();
        std::cout << "  Number of static constraint groups: " << static_constr.size() << std::endl;
        if (static_constr.size() == 3) {
            std::cout << "    Expected: 3 groups (satisfaction + capacity + max facilities)" << std::endl;
        } else {
            std::cerr << "    WARNING: Expected 3 constraint groups" << std::endl;
        }

        // Part 2-3: Register and run both solvers
        std::cout << "\n=== PART 2-3: Setting up and running both solvers ===" << std::endl;
        
        // Load the BlockSolverConfig for both solvers
        BlockSolverConfig* bsc = nullptr;
        {
            auto c = Configuration::deserialize("BSPar.txt");
            bsc = dynamic_cast<BlockSolverConfig*>(c);
            
            if (!bsc) {
                std::cerr << "Error: BSPar.txt does not contain a BlockSolverConfig" << std::endl;
                delete(c);
                return 1;
            }
            
            // Apply the configuration to attach both solvers to the block
            bsc->apply(sr_block);
            
            // Verify solver registration
            auto registered_solvers = sr_block->get_registered_solvers();
            if (registered_solvers.empty()) {
                std::cout << "Error: No solver registered to the block!" << std::endl;
                return 1;
            }
            
            std::cout << "Number of registered solvers: " << registered_solvers.size() << std::endl;
            for (auto solver : registered_solvers) {
                std::cout << "  - " << solver->classname() << std::endl;
            }
        }
        
        // Find the ScenarioReductionSolver in registered solvers
        ScenarioReductionSolver* sr_solver = nullptr;
        for (auto solver : sr_block->get_registered_solvers()) {
            sr_solver = dynamic_cast<ScenarioReductionSolver*>(solver);
            if (sr_solver) break;
        }
        
        if (!sr_solver) {
            std::cerr << "Error: ScenarioReductionSolver not found in registered solvers!" << std::endl;
            return 1;
        }

        // Find the MILPSolver in registered solvers
        Solver* milp_solver = nullptr;
        for (auto solver : sr_block->get_registered_solvers()) {
            // Check if this is not the ScenarioReductionSolver
            if (dynamic_cast<ScenarioReductionSolver*>(solver) == nullptr) {
                milp_solver = solver;
                break;
            }
        }
        
        if (!milp_solver) {
            std::cerr << "Error: MILPSolver not found in registered solvers!" << std::endl;
            return 1;
        }
        
        // Part 4: Test all algorithms
        std::cout << "\n=== PART 4: Testing all scenario reduction algorithms ===" << std::endl;
        
        // Create a vector to store results for comparison
        std::vector<std::tuple<std::string, double, std::vector<int>>> results;
        
        // Test Dupacova's algorithm
        test_algorithm(sr_solver, sr_block, ScenarioReductionSolver::Algorithm::Dupacova, "Dupacova");
        
        // Store Dupacova's results
        double dupacova_distance = pow(sr_solver->get_var_value(), sr_solver->get_ell());
        std::vector<int> dupacova_scenarios;
        for (int i = 0; i < sr_block->get_NFacilities(); i++) {
            if (sr_solver->get_reduced_atoms()[i]) {
                dupacova_scenarios.push_back(i);
            }
        }
        results.push_back({"Dupacova", dupacova_distance, dupacova_scenarios});
        
        // Test BestFit algorithm with rho=0
        test_algorithm(sr_solver, sr_block, ScenarioReductionSolver::Algorithm::BestFit, "BestFit (rho=0)", 0.0);
        
        // Store BestFit results
        double bestfit_distance = pow(sr_solver->get_var_value(), sr_solver->get_ell());
        std::vector<int> bestfit_scenarios;
        for (int i = 0; i < sr_block->get_NFacilities(); i++) {
            if (sr_solver->get_reduced_atoms()[i]) {
                bestfit_scenarios.push_back(i);
            }
        }
        results.push_back({"BestFit (rho=0)", bestfit_distance, bestfit_scenarios});
        
        // Test BestFit algorithm with rho=0.01
        test_algorithm(sr_solver, sr_block, ScenarioReductionSolver::Algorithm::BestFit, "BestFit (rho=0.01)", 0.01);
        
        // Store BestFit results with rho=0.01
        double bestfit_rho_distance = pow(sr_solver->get_var_value(), sr_solver->get_ell());
        std::vector<int> bestfit_rho_scenarios;
        for (int i = 0; i < sr_block->get_NFacilities(); i++) {
            if (sr_solver->get_reduced_atoms()[i]) {
                bestfit_rho_scenarios.push_back(i);
            }
        }
        results.push_back({"BestFit (rho=0.01)", bestfit_rho_distance, bestfit_rho_scenarios});
        
        // Test FirstFit algorithm without shuffle
        test_algorithm(sr_solver, sr_block, ScenarioReductionSolver::Algorithm::FirstFit, "FirstFit (no shuffle)", 0.0, false);
        
        // Store FirstFit results
        double firstfit_distance = pow(sr_solver->get_var_value(), sr_solver->get_ell());
        std::vector<int> firstfit_scenarios;
        for (int i = 0; i < sr_block->get_NFacilities(); i++) {
            if (sr_solver->get_reduced_atoms()[i]) {
                firstfit_scenarios.push_back(i);
            }
        }
        results.push_back({"FirstFit (no shuffle)", firstfit_distance, firstfit_scenarios});
        
        // Test FirstFit algorithm with shuffle
        sr_solver->set_random_seed(42);  // Use fixed seed for reproducibility
        test_algorithm(sr_solver, sr_block, ScenarioReductionSolver::Algorithm::FirstFit, "FirstFit (shuffle)", 0.0, true);
        
        // Store FirstFit results with shuffle
        double firstfit_shuffle_distance = pow(sr_solver->get_var_value(), sr_solver->get_ell());
        std::vector<int> firstfit_shuffle_scenarios;
        for (int i = 0; i < sr_block->get_NFacilities(); i++) {
            if (sr_solver->get_reduced_atoms()[i]) {
                firstfit_shuffle_scenarios.push_back(i);
            }
        }
        results.push_back({"FirstFit (shuffle)", firstfit_shuffle_distance, firstfit_shuffle_scenarios});
        
        // Part 5: Solve with MILPSolver as a reference
        std::cout << "\n=== PART 5: Solving with HiGHSMILPSolver as reference ===" << std::endl;
        
        // Solve with the MILP solver
        std::cout << "Solving scenario reduction problem with MILP solver..." << std::endl;
        int result = milp_solver->compute();
        
        // Print results
        std::cout << "HiGHSMILPSolver returned: " << result << std::endl;
        std::vector<int> milp_scenarios;
        double milp_distance = 0.0;
        
        if (result == Solver::kOK) {
            print_cfl_solution(sr_block);
            
            // Get the MILP solution
            milp_solver->get_var_solution();
            CapacitatedFacilityLocationBlock::CntSolution y(sr_block->get_NFacilities());
            sr_block->get_facility_solution(y.begin());
            
            // Store MILP results
            // The MILP solver reports the total cost, we need to normalize it for fair comparison
            // Dividing by number of facilities to match the way we report other distances
            milp_distance = milp_solver->get_lb() / sr_block->get_NFacilities();
            for (size_t i = 0; i < y.size(); i++) {
                if (y[i] > 0) {
                    milp_scenarios.push_back(i);
                }
            }
            results.push_back({"HiGHSMILPSolver", milp_distance, milp_scenarios});
        } else {
            std::cout << "Solver failed to find a solution." << std::endl;
        }
        
        // Part 6: Compare all algorithm results
        compare_results(results);

        // Clean up
        bsc->clear();
        delete bsc;
        delete sr_block;
        // Note: No need to delete sr_solver individually since it was registered with the block
        // and will be cleaned up when the block is deleted
        
        std::cout << "\nTest completed successfully!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}