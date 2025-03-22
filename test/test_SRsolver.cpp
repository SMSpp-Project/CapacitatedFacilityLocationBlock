/*--------------------------------------------------------------------------*/
/*----------------------- File test_SRsolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test file for the ScenarioReductionSolver class.
 * This is a basic test to verify that the solver compiles and runs properly.
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <vector>

#include "ScenarioReductionSolver.h"
#include "CapacitatedFacilityLocationBlock.h"

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// Function to create a simple CFL problem instance
CapacitatedFacilityLocationBlock* create_test_instance() {
    // Create a small test instance with 3 facilities and 5 customers
    const int n_facilities = 3;
    const int n_customers = 5;
    
    // Create capacity vector for facilities
    std::vector<double> capacities(n_facilities);
    capacities[0] = 100.0;
    capacities[1] = 150.0;
    capacities[2] = 120.0;
    
    // Create fixed costs vector for facilities
    std::vector<double> fixed_costs(n_facilities);
    fixed_costs[0] = 500.0;
    fixed_costs[1] = 600.0;
    fixed_costs[2] = 450.0;
    
    // Create demand vector for customers
    std::vector<double> demands(n_customers);
    demands[0] = 30.0;
    demands[1] = 40.0;
    demands[2] = 25.0;
    demands[3] = 35.0;
    demands[4] = 45.0;
    
    // Create transportation cost matrix
    boost::multi_array<double, 2> transport_costs(boost::extents[n_facilities][n_customers]);
    
    // Fill in transportation costs
    for (int i = 0; i < n_facilities; i++) {
        for (int j = 0; j < n_customers; j++) {
            // Simple costs based on facility and customer indices
            transport_costs[i][j] = 10.0 + 5.0 * i + 3.0 * j;
        }
    }
    
    // Create and load the block
    auto block = new CapacitatedFacilityLocationBlock();
    block->load(n_facilities, n_customers, 
                std::move(capacities), 
                std::move(fixed_costs), 
                std::move(demands), 
                std::move(transport_costs));
    
    return block;
}

// Function to print the solution
void print_solution(ScenarioReductionSolver* solver, CapacitatedFacilityLocationBlock* block) {
    std::cout << "Solution:" << std::endl;
    if (!solver || !block) {
        std::cerr << "Error null ptr!" << std::endl;
        return;
    }
    
    // Get the solution directly from the solver
    const auto& facility_solution = solver->get_facility_solution();
    const auto& transport_solution = solver->get_transportation_solution();
    
    // Print which facilities are open
    std::cout << "Open facilities: ";
    for (int i = 0; i < block->get_NFacilities(); i++) {
        if (facility_solution[i]) {
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
    
    // Print transportation assignments
    std::cout << "Transportation assignments:" << std::endl;
    for (int i = 0; i < block->get_NFacilities(); i++) {
        for (int j = 0; j < block->get_NCustomers(); j++) {
            double value = transport_solution[i * block->get_NCustomers() + j];
            if (value > 0.001) { // Only print non-zero assignments
                std::cout << "  Facility " << i << " -> Customer " << j 
                          << ": " << value * 100 << "%" << std::endl;
            }
        }
    }
}

/*--------------------------------------------------------------------------*/
/*------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main() {
    try {
        std::cout << "Creating test instance..." << std::endl;
        auto block = create_test_instance();
        
        std::cout << "Creating solver..." << std::endl;
        auto solver = new ScenarioReductionSolver();
        
        std::cout << "Registering solver with block..." << std::endl;
        block->register_Solver(solver);
        
        std::cout << "Solving problem..." << std::endl;
        int result = solver->compute();
        
        std::cout << "Solver returned: " << result << std::endl;
        if (result == Solver::kOK) {
            std::cout << "Objective value: " << solver->get_var_value() << std::endl;
            print_solution(solver, block);
        } else {
            std::cout << "Solver failed to find a solution." << std::endl;
        }
        
        // Test a modification
        std::cout << "\nTesting modification - changing facility cost..." << std::endl;
        block->chg_facility_cost(700.0, 0);  // Change cost of facility 0
        
        std::cout << "Resolving problem..." << std::endl;
        result = solver->compute();
        
        std::cout << "Solver returned: " << result << std::endl;
        if (result == Solver::kOK) {
            std::cout << "Objective value: " << solver->get_var_value() << std::endl;
            print_solution(solver, block);
        } else {
            std::cout << "Solver failed to find a solution." << std::endl;
        }
        
        
        // Clean up
        delete solver;
        delete block;
        
        std::cout << "Test completed successfully!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_SRsolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/