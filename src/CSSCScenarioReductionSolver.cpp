/*--------------------------------------------------------------------------*/
/*-------------- File CSSCScenarioReductionSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of CSSCScenarioReductionSolver.
 *
 * \author Minh Duc Pham \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

#include "CSSCScenarioReductionSolver.h"
#include "AbstractBlock.h"
#include "ColVariable.h"
#include "DiscreteScenarioSet.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "StochasticBlock.h"

#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <unordered_set>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && (LogVerb > l)) *f_log << x

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------- CSSCComputeConfig METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void CSSCComputeConfig::serialize( netCDF::NcGroup & group ) const
{
 // 1. Let base class handle f_extra_Configuration (BlockSolverConfig)
 //    and all int/dbl/str/vint parameters
 ComputeConfig::serialize( group );

 // 2. Serialize DiscreteScenarioSet into sub-group "ScenarioSet"
 //    Only if we have one to serialize
 if( f_scenario_set ) {
  auto sg = group.addGroup( "ScenarioSet" );
  f_scenario_set->serialize( sg );
  }
}

/*--------------------------------------------------------------------------*/

void CSSCComputeConfig::deserialize( const netCDF::NcGroup & group )
{
 // 1. Let base class read f_extra_Configuration + all parameters
 ComputeConfig::deserialize( group );

 // 2. Clear any previously owned DSS
 delete f_owned_dss;
 f_owned_dss    = nullptr;
 f_scenario_set = nullptr;

 // 3. Try to read DiscreteScenarioSet from sub-group "ScenarioSet"
 auto sg = group.getGroup( "ScenarioSet" );
 if( ! sg.isNull() ) {
  f_owned_dss = new DiscreteScenarioSet();
  f_owned_dss->deserialize( sg );
  f_scenario_set = f_owned_dss;
  }
}

/*--------------------------------------------------------------------------*/

void CSSCComputeConfig::load( std::istream & input )
{
 // 1. Let base class read f_extra_Configuration + all parameters from txt
 ComputeConfig::load( input );

 // 2. Clear any previously owned DSS
 delete f_owned_dss;
 f_owned_dss    = nullptr;
 f_scenario_set = nullptr;

 // 3. Try to read DiscreteScenarioSet from stream
 //    Format: peek at next token, if it is an integer, assume DSS follows
 //    The DSS load() format starts with "N D" (two positive integers)
 std::streampos pos = input.tellg();
 unsigned int N = 0;
 if( ( input >> N ) && N > 0 ) {
  // Put back and let DSS load() read from here
  input.seekg( pos );
  f_owned_dss = new DiscreteScenarioSet();
  f_owned_dss->load( input );
  f_scenario_set = f_owned_dss;
  }
 else {
  // No DSS data in stream, clear error flags and leave f_scenario_set null
  input.clear();
  input.seekg( pos );
  }
}

/*--------------------------------------------------------------------------*/

CSSCComputeConfig * CSSCComputeConfig::clone( void ) const
{
 auto * c = new CSSCComputeConfig();
 // Copy base class fields
 c->f_diff  = f_diff;
 c->f_relax = f_relax;
 c->int_pars  = int_pars;
 c->dbl_pars  = dbl_pars;
 c->str_pars  = str_pars;
 c->vint_pars = vint_pars;
 c->vdbl_pars = vdbl_pars;
 c->vstr_pars = vstr_pars;
 c->f_extra_Configuration =
   f_extra_Configuration ? f_extra_Configuration->clone() : nullptr;
 // Copy DSS: if we own it, reconstruct via load_from_memory()
 // (DSS has unique_ptr field so copy constructor is deleted)
 // if caller-provided (not owned), just copy the pointer
 if( f_owned_dss ) {
  const auto N = f_owned_dss->get_nbScenarios();
  const auto D = f_owned_dss->get_scenario_size();
  std::vector< std::vector< double > > scenarios( N ,
                                                   std::vector< double >( D ) );
  for( unsigned int i = 0 ; i < N ; ++i )
   for( unsigned int d = 0 ; d < D ; ++d )
    scenarios[ i ][ d ] = f_owned_dss->get_scenario_value( i , d );
  std::vector< double > weights( f_owned_dss->get_set_weights().begin() ,
                                 f_owned_dss->get_set_weights().end() );
  c->f_owned_dss = new DiscreteScenarioSet();
  c->f_owned_dss->load_from_memory( scenarios , weights );
  c->f_scenario_set = c->f_owned_dss;
  }
 else {
  c->f_owned_dss    = nullptr;
  c->f_scenario_set = f_scenario_set;  // non-owning copy
  }
 return c;
}

/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( CSSCScenarioReductionSolver );

/*--------------------------------------------------------------------------*/
/*----------------------- CONFIGURATION METHOD -----------------------------*/
/*--------------------------------------------------------------------------*/

void CSSCScenarioReductionSolver::set_ComputeConfig( const ComputeConfig * cfg ) {
 // 1. Let base class handle standard int/dbl/vint parameters
 ScenarioReductionSolver::set_ComputeConfig( cfg );

 // 2. Clear previously stored config and scenario pointer.
 delete f_milp_config;
 f_milp_config  = nullptr;
 f_scenario_set = nullptr;

 if( ! cfg ) return;

 // 3. Cast to CSSCComputeConfig to access CSSC-specific fields
 auto * cssc_cfg = dynamic_cast< const CSSCComputeConfig * >( cfg );
 if( ! cssc_cfg )
  throw std::invalid_argument(
    "CSSCScenarioReductionSolver::set_ComputeConfig: expected a "
    "CSSCComputeConfig. Use CSSCComputeConfig to pass BlockSolverConfig "
    "and DiscreteScenarioSet to this solver." );

 // 4. Clone the BlockSolverConfig from f_extra_Configuration (owned by us)
 if( cssc_cfg->f_extra_Configuration ) {
  auto * bsc = dynamic_cast< BlockSolverConfig * >(
    cssc_cfg->f_extra_Configuration );
  if( ! bsc )
   throw std::invalid_argument(
     "CSSCScenarioReductionSolver::set_ComputeConfig: "
     "f_extra_Configuration must be a BlockSolverConfig." );
  f_milp_config = static_cast< BlockSolverConfig * >( bsc->clone() );
 }

 // 5. Store scenario set pointer (not owned, caller retains ownership)
 f_scenario_set = cssc_cfg->f_scenario_set;
}

/*--------------------------------------------------------------------------*/
/*------------------------ MAIN COMPUTATION METHOD -------------------------*/
/*--------------------------------------------------------------------------*/

int CSSCScenarioReductionSolver::compute( bool changedvars ) {
 std::lock_guard< std::recursive_mutex > lock( f_mutex );

 if( ! get_Block() ) return kError;

 mod_clear();

 // clear previous solution state.
 indices_to_choose.clear();
 ind_red.clear();
 std::fill( reduced_atoms.begin() , reduced_atoms.end() , false );

 BLOG( 0 , "\nCSSCScenarioReductionSolver: reducing "
           << nb_atoms << " scenarios to " << nb_reduced
           << " using CSSC algorithm" << std::endl );

 // Guard 1: MILP solver config is required for both Step 1 sub-problems
 //          and the Step 2 partitioning MILP
 if( ! f_milp_config )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute: no MILP solver config set. "
    "Call set_ComputeConfig() with a valid extra configuration." );

 // Guard 2: scenario set is required for Step 1 (building V matrix)
 if( ! f_scenario_set )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute: no scenario set provided. "
    "Call set_ComputeConfig() with a valid extra configuration." );

 // Step 1: build the N×N opportunity-cost matrix V
 const auto V = compute_V_matrix();

 // Step 2: solve the MILP and populate ind_red / reduced_atoms /
 //         f_solution_value
 solve_cssc_milp( V );

 return kOK;
}

/*--------------------------------------------------------------------------*/
/*-------------------- STEP 1: compute_V_matrix() --------------------------*/
/*--------------------------------------------------------------------------*/
// V[i][j] = F( x*_i , xi_j )                                    eq.(21)
//
// x*_i = argmin_{x in X} F(x, xi_i)                              eq.(22)
//
// procedure for each scenario i:
//   (a) inject xi_i via StochasticBlock::set_data() -> solve CFL -> x*_i
//   (b) for each j!=i: fix y=x*_i, inject xi_j, solve -> V[i][j]
//   (c) restore y to kBinary before next iteration
/*--------------------------------------------------------------------------*/

std::vector< std::vector< double > >
CSSCScenarioReductionSolver::compute_V_matrix() {

 const int n = static_cast< int >( nb_atoms );

 // The CFL block must be the inner block of a StochasticBlock so that
 // set_data() can inject individual scenarios into its uncertain parameters
 auto * stoch_block = dynamic_cast< StochasticBlock * >(
   f_Block->get_f_Block() );
 if( ! stoch_block )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute_V_matrix: f_Block has no "
    "StochasticBlock parent. CSSC Step 1 requires a StochasticBlock wrapper "
    "so that individual scenarios can be injected via set_data()." );

 auto * cfl_block = dynamic_cast< CapacitatedFacilityLocationBlock * >(
   f_Block );

 // Attach a cloned MILP solver to f_Block for the sub-problem solves
 // I clone so the original f_milp_config remains available for Step 2
 // After Step 1 we call clear() to properly detach the solver before
 // the BlockSolverConfig goes out of scope, as required by SMS++ design
 BlockSolverConfig * sub_cfg =
   static_cast< BlockSolverConfig * >( f_milp_config->clone() );
 sub_cfg->apply( f_Block );

 Solver * sub_solver = f_Block->get_registered_solvers().empty()
   ? nullptr : f_Block->get_registered_solvers().front();
 if( ! sub_solver ) {
  sub_cfg->clear();
  delete sub_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::compute_V_matrix: failed to attach "
    "sub-problem solver to the CFL block." );
 }

 const auto scenario_size = f_scenario_set->get_scenario_size();

 std::vector< std::vector< double > > V( n , std::vector< double >( n , 0.0 ));
 std::vector< std::vector< double > > x_star( n , std::vector< double >( n , 0.0 ));

 BLOG( 1 , "\n  CSSC Step 1: computing " << n << "x" << n
           << " V matrix (" << (n + n*(n-1)) << " CFL solves)" << std::endl );

 for( int i = 0 ; i < n ; ++i ) {

  // (a) load scenario i and solve
  {
   std::vector< double > data( scenario_size );
   for( std::size_t d = 0 ; d < scenario_size ; ++d )
    data[ d ] = f_scenario_set->get_scenario_value(
      static_cast< DiscreteScenarioSet::ScenarioIndex >( i ) ,
      static_cast< DiscreteScenarioSet::ScenarioSize   >( d ) );
   stoch_block->set_data( data , eNoBlck , eNoBlck );
  }
  {
   const int st = sub_solver->compute();
   if( st != Solver::kOK && st != Solver::kLowPrecision ) {
    sub_cfg->clear();  delete sub_cfg;
    throw std::runtime_error(
      "CSSC Step 1: sub-problem failed for scenario " + std::to_string( i ) +
      " (status " + std::to_string( st ) + ")" );
   }
   sub_solver->get_var_solution();
  }
  for( int f = 0 ; f < n ; ++f ) {
   ColVariable * yv = cfl_block->get_y( static_cast< Index >( f ) );
   x_star[ i ][ f ] = yv ? yv->get_value() : 0.0;
  }
  V[ i ][ i ] = sub_solver->get_var_value();

  // (b) fix y = x*_i, evaluate under each other scenario j
  for( int j = 0 ; j < n ; ++j ) {
   if( j == i ) continue;

   for( int f = 0 ; f < n ; ++f ) {
    ColVariable * yv = cfl_block->get_y( static_cast< Index >( f ) );
    if( yv ) {
     const double val = x_star[ i ][ f ];
     yv->set_value( val );
     yv->set_type( val > 0.5 ? ColVariable::kZeroIntU
                             : ColVariable::kZeroReal , eNoMod );
    }
   }
   {
    std::vector< double > data_j( scenario_size );
    for( std::size_t d = 0 ; d < scenario_size ; ++d )
     data_j[ d ] = f_scenario_set->get_scenario_value(
       static_cast< DiscreteScenarioSet::ScenarioIndex >( j ) ,
       static_cast< DiscreteScenarioSet::ScenarioSize   >( d ) );
    stoch_block->set_data( data_j , eNoBlck , eNoBlck );
   }
   {
    const int st = sub_solver->compute();
    if( st != Solver::kOK && st != Solver::kLowPrecision ) {
     sub_cfg->clear();  delete sub_cfg;
     throw std::runtime_error(
       "CSSC Step 1: evaluation failed for i=" + std::to_string( i ) +
       " j=" + std::to_string( j ) );
    }
    sub_solver->get_var_solution();
   }
   V[ i ][ j ] = sub_solver->get_var_value();

   // (c) restore y to binary for next iteration
   for( int f = 0 ; f < n ; ++f ) {
    ColVariable * yv = cfl_block->get_y( static_cast< Index >( f ) );
    if( yv ) yv->set_type( ColVariable::kBinary , eNoMod );
   }
  }

  BLOG( 1 , "  Step 1: scenario " << (i+1) << "/" << n
            << "  V[i][i]=" << V[i][i] << std::endl );
 }

 // properly detach the solver from f_Block before deleting the config
 sub_cfg->clear();
 delete sub_cfg;

 BLOG( 1 , "  CSSC Step 1: V matrix complete." << std::endl );
 return V;
}

/*--------------------------------------------------------------------------*/
/*------------------ STEP 2: solve_cssc_milp(V) ----------------------------*/
/*--------------------------------------------------------------------------*/
// MILP partitioning problem (eqs 24-29, Keutchayan et al. 2023):
//
//   min   (1/N) sum_j t_j
//   s.t.  t_j >= sum_j_i x_ij*(V[j][i] - V[j][j])   for all j   (25)
//         t_j >= sum_j_i x_ij*(V[j][j] - V[j][i])   for all j   (26)
//         x_ij <= u_j                               for all i,j (27a)
//         x_jj  = u_j                               for all j   (27b)
//         sum_j_j x_ij = 1                          for all i   (28a)
//         sum_j_j u_j  = K                                      (28b)
/*--------------------------------------------------------------------------*/

void CSSCScenarioReductionSolver::solve_cssc_milp(
  const std::vector< std::vector< double > > & V ) {

 const int n = static_cast< int >( nb_atoms );
 const int K = static_cast< int >( nb_reduced );

 BLOG( 1 , "  CSSC Step 2: building MILP (N=" << n << ", K=" << K
           << ")" << std::endl );

 auto milp_block = std::make_unique< AbstractBlock >();

 // AbstractBlock::~AbstractBlock can safely
 // delete them (delete on stack objects would segfault)
 auto * x_vars_p = new std::vector< ColVariable >( n * n );
 auto * u_vars_p = new std::vector< ColVariable >( n );
 auto * t_vars_p = new std::vector< ColVariable >( n );
 auto & x_vars = *x_vars_p;
 auto & u_vars = *u_vars_p;
 auto & t_vars = *t_vars_p;

 for( int i = 0 ; i < n * n ; ++i ) x_vars[ i ].set_type( ColVariable::kBinary );
 for( int j = 0 ; j < n    ; ++j ) u_vars[ j ].set_type( ColVariable::kBinary );
 for( int j = 0 ; j < n    ; ++j ) t_vars[ j ].set_type( ColVariable::kNonNegative );

 milp_block->add_static_variable( x_vars );
 milp_block->add_static_variable( u_vars );
 milp_block->add_static_variable( t_vars );

 // eq.(25): sum_i (V[j][i]-V[j][j])*x_ij - t_j <= 0
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  for( int i = 0 ; i < n ; ++i ) {
   const double c = V[j][i] - V[j][j];
   if( std::abs(c) > 1e-15 ) terms.emplace_back( &x_vars[i*n+j] , c );
  }
  terms.emplace_back( &t_vars[j] , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( -Inf<double>() );  con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(26): sum_i (V[j][j]-V[j][i])*x_ij - t_j <= 0
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  for( int i = 0 ; i < n ; ++i ) {
   const double c = V[j][j] - V[j][i];
   if( std::abs(c) > 1e-15 ) terms.emplace_back( &x_vars[i*n+j] , c );
  }
  terms.emplace_back( &t_vars[j] , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( -Inf<double>() );  con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(27a): x_ij - u_j <= 0
 for( int i = 0 ; i < n ; ++i )
  for( int j = 0 ; j < n ; ++j ) {
   LinearFunction::v_coeff_pair terms;
   terms.emplace_back( &x_vars[i*n+j] ,  1.0 );
   terms.emplace_back( &u_vars[j]     , -1.0 );
   auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
   auto con = std::make_unique< FRowConstraint >();
   con->set_lhs( -Inf<double>() );  con->set_rhs( 0.0 );
   con->set_function( lf.release() , eNoMod );
   milp_block->add_static_constraint( *con.release() );
  }

 // eq.(27b): x_jj - u_j = 0
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  terms.emplace_back( &x_vars[j*n+j] ,  1.0 );
  terms.emplace_back( &u_vars[j]     , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( 0.0 );  con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(28a): sum_j x_ij = 1
 for( int i = 0 ; i < n ; ++i ) {
  LinearFunction::v_coeff_pair terms;
  for( int j = 0 ; j < n ; ++j )
   terms.emplace_back( &x_vars[i*n+j] , 1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( 1.0 );  con->set_rhs( 1.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(28b): sum_j u_j = K
 {
  LinearFunction::v_coeff_pair terms;
  for( int j = 0 ; j < n ; ++j )
   terms.emplace_back( &u_vars[j] , 1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( static_cast<double>(K) );
  con->set_rhs( static_cast<double>(K) );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // objective: min (1/N) sum_j t_j
 {
  const double c = 1.0 / static_cast<double>(n);
  LinearFunction::v_coeff_pair terms;
  for( int j = 0 ; j < n ; ++j )
   terms.emplace_back( &t_vars[j] , c );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto obj = std::make_unique< FRealObjective >( milp_block.get() , lf.release() );
  obj->set_sense( Objective::eMin );
  milp_block->set_objective( obj.release() );
 }

 BLOG( 1 , "  CSSC Step 2: solving partitioning MILP" << std::endl );

 milp_block->generate_abstract_variables();

 // Attach MILP solver, solve, then detach via clear() before
 // the config goes out of scope, required by BlockSolverConfig design
 BlockSolverConfig * milp_cfg =
   static_cast< BlockSolverConfig * >( f_milp_config->clone() );
 milp_cfg->apply( milp_block.get() );

 Solver * milp_solver = milp_block->get_registered_solvers().empty()
   ? nullptr : milp_block->get_registered_solvers().front();
 if( ! milp_solver ) {
  milp_cfg->clear();  delete milp_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::solve_cssc_milp: no solver attached." );
 }

 const int status = milp_solver->compute();
 if( status != Solver::kOK && status != Solver::kLowPrecision ) {
  milp_cfg->clear();  delete milp_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::solve_cssc_milp: MILP failed (status " +
    std::to_string(status) + ")" );
 }

 milp_solver->get_var_solution();

 BLOG( 0 , "\n  CSSC MILP solved, discrepancy = "
           << std::fixed << std::setprecision(6)
           << milp_solver->get_var_value() << std::endl );

 // detach solver before config is deleted
 milp_cfg->clear();
 delete milp_cfg;

 // extract solution into shared ScenarioReductionSolver output fields
 ind_red.clear();
 for( int j = 0 ; j < n ; ++j )
  if( u_vars[j].get_value() > 0.5 )
   ind_red.push_back( static_cast<Index>(j) );

 // numerical safety: if MILP returned wrong count, take top-K by u value
 if( static_cast<int>(ind_red.size()) != K ) {
  std::vector<int> all(n);
  std::iota( all.begin() , all.end() , 0 );
  std::sort( all.begin() , all.end() , [&](int a, int b){
   return u_vars[a].get_value() > u_vars[b].get_value();
  });
  ind_red.clear();
  for( int k = 0 ; k < K ; ++k )
   ind_red.push_back( static_cast<Index>(all[k]) );
  std::sort( ind_red.begin() , ind_red.end() );
 }

 std::unordered_set<Index> rep_set( ind_red.begin() , ind_red.end() );
 indices_to_choose.clear();
 for( int i = 0 ; i < n ; ++i )
  if( ! rep_set.count(static_cast<Index>(i)) )
   indices_to_choose.push_back( static_cast<Index>(i) );

 update_reduced_atoms();

 // Wasserstein distance (same formula as Dupacova/BestFit)
 double wasserstein = 0.0;
 for( int i = 0 ; i < n ; ++i ) {
  double min_c = std::numeric_limits<double>::infinity();
  for( auto r : ind_red )
   min_c = std::min( min_c , (*f_transportation_costs)[i][r] );
  wasserstein += (*weights)[i] * min_c;
 }
 f_solution_value = wasserstein;

 BLOG( 0 , "  CSSC: Wasserstein distance = "
           << std::fixed << std::setprecision(6)
           << f_solution_value << std::endl );
}

/*--------------------------------------------------------------------------*/
/*---------- End File CSSCScenarioReductionSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/