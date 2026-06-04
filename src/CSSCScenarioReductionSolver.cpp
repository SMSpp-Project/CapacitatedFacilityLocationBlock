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
/*--------------------------------------------------------------------------*/

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
#include <unordered_set>

/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && (LogVerb > l)) *f_log << x

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( CSSCScenarioReductionSolver );

/*--------------------------------------------------------------------------*/
/*---------------- CSSCComputeConfig METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void CSSCComputeConfig::serialize( netCDF::NcGroup & group ) const
{
 ComputeConfig::serialize( group );
 if( f_scenario_set ) {
  auto sg = group.addGroup( "ScenarioSet" );
  f_scenario_set->serialize( sg );
 }
}

/*--------------------------------------------------------------------------*/

void CSSCComputeConfig::deserialize( const netCDF::NcGroup & group )
{
 ComputeConfig::deserialize( group );
 delete f_owned_dss;
 f_owned_dss    = nullptr;
 f_scenario_set = nullptr;
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
 ComputeConfig::load( input );
 delete f_owned_dss;
 f_owned_dss    = nullptr;
 f_scenario_set = nullptr;
 std::streampos pos = input.tellg();
 unsigned int N = 0;
 if( ( input >> N ) && N > 0 ) {
  input.seekg( pos );
  f_owned_dss = new DiscreteScenarioSet();
  f_owned_dss->load( input );
  f_scenario_set = f_owned_dss;
 } else {
  input.clear();
  input.seekg( pos );
 }
}

/*--------------------------------------------------------------------------*/

CSSCComputeConfig * CSSCComputeConfig::clone() const
{
 auto * c = new CSSCComputeConfig();
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
 if( f_owned_dss ) {
  const auto Ns = f_owned_dss->get_nbScenarios();
  const auto D  = f_owned_dss->get_scenario_size();
  std::vector< std::vector< double > > sc( Ns , std::vector< double >( D ) );
  for( unsigned int i = 0 ; i < Ns ; ++i )
   for( unsigned int d = 0 ; d < D  ; ++d )
    sc[ i ][ d ] = f_owned_dss->get_scenario_value( i , d );
  std::vector< double > wt( f_owned_dss->get_set_weights().begin() ,
                             f_owned_dss->get_set_weights().end() );
  c->f_owned_dss = new DiscreteScenarioSet();
  c->f_owned_dss->load_from_memory( sc , wt );
  c->f_scenario_set = c->f_owned_dss;
 } else {
  c->f_owned_dss    = nullptr;
  c->f_scenario_set = f_scenario_set;
 }
 return c;
}

/*--------------------------------------------------------------------------*/
/*----------------------- CONFIGURATION METHOD -----------------------------*/
/*--------------------------------------------------------------------------*/

void CSSCScenarioReductionSolver::set_ComputeConfig( const ComputeConfig * cfg )
{
 ScenarioReductionSolver::set_ComputeConfig( cfg );

 delete f_milp_config;       f_milp_config       = nullptr;
 delete f_milp_config_step2; f_milp_config_step2 = nullptr;
 f_scenario_set = nullptr;

 if( ! cfg ) return;

 auto * cssc_cfg = dynamic_cast< const CSSCComputeConfig * >( cfg );
 if( ! cssc_cfg )
  throw std::invalid_argument(
    "CSSCScenarioReductionSolver::set_ComputeConfig: expected a "
    "CSSCComputeConfig." );

 // f_extra_Configuration -> Step 1 sub-problem solver (LP relaxation)
 if( cssc_cfg->f_extra_Configuration ) {
  auto * bsc = dynamic_cast< BlockSolverConfig * >(
    cssc_cfg->f_extra_Configuration );
  if( ! bsc )
   throw std::invalid_argument(
     "CSSCScenarioReductionSolver::set_ComputeConfig: "
     "f_extra_Configuration must be a BlockSolverConfig." );
  f_milp_config = static_cast< BlockSolverConfig * >( bsc->clone() );
 }

 // f_milp_config -> Step 2 partitioning MILP solver (optional, MIP)
 // If not provided, Step 2 reuses f_milp_config (Step 1 config)
 if( cssc_cfg->f_milp_config ) {
  f_milp_config_step2 =
    static_cast< BlockSolverConfig * >( cssc_cfg->f_milp_config->clone() );
 }

 f_scenario_set = cssc_cfg->f_scenario_set;
}

/*--------------------------------------------------------------------------*/
/*------------------------ MAIN COMPUTATION METHOD -------------------------*/
/*--------------------------------------------------------------------------*/

int CSSCScenarioReductionSolver::compute( bool changedvars )
{
 std::lock_guard< std::recursive_mutex > lock( f_mutex );

 if( ! get_Block() ) return kError;

 mod_clear();

 indices_to_choose.clear();
 ind_red.clear();
 std::fill( reduced_atoms.begin() , reduced_atoms.end() , false );

 if( ! f_milp_config )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute: no MILP solver config. "
    "Call set_ComputeConfig() with a CSSCComputeConfig." );

 if( ! f_scenario_set )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute: no scenario set. "
    "Call set_ComputeConfig() with a CSSCComputeConfig." );

 if( ! f_sub_block )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute: no sub-problem block. "
    "Call set_sub_problem_block() with the real CFL block." );

 BLOG( 0 , "\nCSSCScenarioReductionSolver: reducing "
           << nb_atoms << " scenarios to " << nb_reduced
           << " using CSSC algorithm" << std::endl );

 // Step 1: build the N×N opportunity-cost matrix V
 const auto V = compute_V_matrix();

 // Step 2: solve the MILP partitioning problem
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
// Procedure for each scenario i:
//   (a) inject xi_i via StochasticBlock::set_data() -> solve -> x*_i
//   (b) fix y = x*_i once before the j-loop
//   (c) for each j != i: inject xi_j -> solve -> V[i][j]
//   (d) restore y to kBinary once AFTER the j-loop
//
// f_sub_block is the real CFL block. Its parent (get_f_Block()) must
// be a StochasticBlock for scenario injection.
/*--------------------------------------------------------------------------*/

std::vector< std::vector< double > >
CSSCScenarioReductionSolver::compute_V_matrix()
{
 const int n  = static_cast< int >( nb_atoms );
 const int nf = static_cast< int >( f_sub_block->get_NFacilities() );

 // f_sub_block's parent must be a StochasticBlock
 auto * stoch = dynamic_cast< StochasticBlock * >(
   f_sub_block->get_f_Block() );
 if( ! stoch )
  throw std::logic_error(
    "CSSCScenarioReductionSolver::compute_V_matrix: "
    "f_sub_block has no StochasticBlock parent. "
    "Call base_cfl->set_f_Block(stoch) before compute()." );

 // Attach a cloned MILP solver to f_sub_block.
 // Clone so f_milp_config remains available for Step 2.
 // After Step 1 call clear() to detach before config is deleted.
 BlockSolverConfig * sub_cfg =
   static_cast< BlockSolverConfig * >( f_milp_config->clone() );
 sub_cfg->apply( f_sub_block );

 Solver * sub_solver = f_sub_block->get_registered_solvers().empty()
   ? nullptr : f_sub_block->get_registered_solvers().front();
 if( ! sub_solver ) {
  sub_cfg->clear(); delete sub_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::compute_V_matrix: "
    "failed to attach sub-solver to f_sub_block." );
 }

 const auto scenario_size = f_scenario_set->get_scenario_size();

 std::vector< std::vector< double > > V( n , std::vector< double >( n , 0.0 ) );
 std::vector< std::vector< double > > x_star( n , std::vector< double >( nf , 0.0 ) );

 BLOG( 1 , "\n  CSSC Step 1: computing " << n << "x" << n
           << " V matrix (" << (n + n*(n-1)) << " sub-problem solves)"
           << std::endl );

 for( int i = 0 ; i < n ; ++i ) {

  // (a) inject scenario i and solve -> x*_i
  //     y variables are relaxed to continuous [0,1] for speed.
  //     Relax once before solve, restore to binary after reading x*_i.
  {
   std::vector< double > data( scenario_size );
   for( std::size_t d = 0 ; d < scenario_size ; ++d )
    data[ d ] = f_scenario_set->get_scenario_value(
      static_cast< DiscreteScenarioSet::ScenarioIndex >( i ) ,
      static_cast< DiscreteScenarioSet::ScenarioSize   >( d ) );
   stoch->set_data( data , eNoBlck , eNoBlck );
  }
  // Relax y to continuous [0,1]: LP solve instead of MIP
  for( int f = 0 ; f < nf ; ++f ) {
   ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
   if( yv ) yv->set_type( ColVariable::kContinuous , eModBlck );
  }
  {
   const int st = sub_solver->compute();
   if( st != Solver::kOK && st != Solver::kLowPrecision ) {
    // restore before throwing
    for( int f = 0 ; f < nf ; ++f ) {
     ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
     if( yv ) yv->set_type( ColVariable::kBinary , eModBlck );
    }
    sub_cfg->clear(); delete sub_cfg;
    throw std::runtime_error(
      "CSSC Step 1: sub-problem failed for scenario " +
      std::to_string( i ) + " (status " + std::to_string( st ) + ")" );
   }
   sub_solver->get_var_solution();
  }
  // Restore y to binary after reading x*_i
  for( int f = 0 ; f < nf ; ++f ) {
   ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
   if( yv ) yv->set_type( ColVariable::kBinary , eModBlck );
  }

  // Read x*_i from f_sub_block y variables
  for( int f = 0 ; f < nf ; ++f ) {
   ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
   x_star[ i ][ f ] = yv ? yv->get_value() : 0.0;
  }
  V[ i ][ i ] = sub_solver->get_var_value();

  // (b) fix y = x*_i ONCE before the j-loop.
  //     y is already restored to kBinary above.
  //     Now fix: set_value() + is_fixed(true) + relax to continuous
  //     so solver treats this as LP (fixed y = constants, solve x LP).
  for( int f = 0 ; f < nf ; ++f ) {
   ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
   if( yv ) {
    yv->set_value( x_star[ i ][ f ] );
    yv->is_fixed( true , eModBlck );
   }
  }

  // (c) evaluate x*_i under each other scenario j
  for( int j = 0 ; j < n ; ++j ) {
   if( j == i ) continue;

   {
    std::vector< double > data_j( scenario_size );
    for( std::size_t d = 0 ; d < scenario_size ; ++d )
     data_j[ d ] = f_scenario_set->get_scenario_value(
       static_cast< DiscreteScenarioSet::ScenarioIndex >( j ) ,
       static_cast< DiscreteScenarioSet::ScenarioSize   >( d ) );
    stoch->set_data( data_j , eNoBlck , eNoBlck );
   }
   {
    const int st = sub_solver->compute();
    if( st != Solver::kOK && st != Solver::kLowPrecision ) {
     // restore before throwing so the block is left clean
     for( int f = 0 ; f < nf ; ++f ) {
      ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
      if( yv ) { yv->is_fixed( false , eModBlck );
                 yv->set_type( ColVariable::kBinary , eModBlck ); }
     }
     sub_cfg->clear(); delete sub_cfg;
     throw std::runtime_error(
       "CSSC Step 1: evaluation failed for i=" + std::to_string( i ) +
       " j=" + std::to_string( j ) );
    }
    sub_solver->get_var_solution();
   }
   V[ i ][ j ] = sub_solver->get_var_value();
  }

  // (d) restore y to kBinary ONCE after the j-loop
  for( int f = 0 ; f < nf ; ++f ) {
   ColVariable * yv = f_sub_block->get_y( static_cast< Index >( f ) );
   if( yv ) {
    yv->is_fixed( false , eModBlck );
    yv->set_type( ColVariable::kBinary , eModBlck );
   }
  }

  BLOG( 1 , "  Step 1: scenario " << (i+1) << "/" << n
            << "  V[i][i]=" << V[i][i] << std::endl );
 }

 sub_cfg->clear();
 delete sub_cfg;

 BLOG( 1 , "  CSSC Step 1: V matrix complete." << std::endl );
 return V;
}

/*--------------------------------------------------------------------------*/
/*------------------ STEP 2: solve_cssc_milp(V) ----------------------------*/
/*--------------------------------------------------------------------------*/
// MILP partitioning problem (Keutchayan et al. 2023, eqs 24-29):
//
//   min   (1/N) sum_j t_j
//   s.t.  t_j >= sum_i (V[j][i]-V[j][j]) * x_ij   for all j   (25)
//         t_j >= sum_i (V[j][j]-V[j][i]) * x_ij   for all j   (26)
//         x_ij <= u_j                               for all i,j (27a)
//         x_jj  = u_j                               for all j   (27b)
//         sum_j x_ij = 1                            for all i   (28a)
//         sum_j u_j  = K                                        (28b)
/*--------------------------------------------------------------------------*/

void CSSCScenarioReductionSolver::solve_cssc_milp(
  const std::vector< std::vector< double > > & V )
{
 const int n = static_cast< int >( nb_atoms );
 const int K = static_cast< int >( nb_reduced );

 BLOG( 1 , "  CSSC Step 2: building MILP (N=" << n << ", K=" << K
           << ")" << std::endl );

 auto milp_block = std::make_unique< AbstractBlock >();

 // Heap-allocated so AbstractBlock::~AbstractBlock can safely delete them
 auto * x_vars_p = new std::vector< ColVariable >( n * n );
 auto * u_vars_p = new std::vector< ColVariable >( n );
 auto * t_vars_p = new std::vector< ColVariable >( n );
 auto & x_vars = *x_vars_p;
 auto & u_vars = *u_vars_p;
 auto & t_vars = *t_vars_p;

 for( int i = 0 ; i < n * n ; ++i ) x_vars[i].set_type( ColVariable::kBinary );
 for( int j = 0 ; j < n    ; ++j ) u_vars[j].set_type( ColVariable::kBinary );
 for( int j = 0 ; j < n    ; ++j ) t_vars[j].set_type( ColVariable::kNonNegative );

 milp_block->add_static_variable( x_vars );
 milp_block->add_static_variable( u_vars );
 milp_block->add_static_variable( t_vars );

 // eq.(25)
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  for( int i = 0 ; i < n ; ++i ) {
   const double c = V[j][i] - V[j][j];
   if( std::abs(c) > 1e-15 ) terms.emplace_back( &x_vars[i*n+j] , c );
  }
  terms.emplace_back( &t_vars[j] , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( -Inf<double>() ); con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(26)
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  for( int i = 0 ; i < n ; ++i ) {
   const double c = V[j][j] - V[j][i];
   if( std::abs(c) > 1e-15 ) terms.emplace_back( &x_vars[i*n+j] , c );
  }
  terms.emplace_back( &t_vars[j] , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( -Inf<double>() ); con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(27a)
 for( int i = 0 ; i < n ; ++i )
  for( int j = 0 ; j < n ; ++j ) {
   LinearFunction::v_coeff_pair terms;
   terms.emplace_back( &x_vars[i*n+j] ,  1.0 );
   terms.emplace_back( &u_vars[j]     , -1.0 );
   auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
   auto con = std::make_unique< FRowConstraint >();
   con->set_lhs( -Inf<double>() ); con->set_rhs( 0.0 );
   con->set_function( lf.release() , eNoMod );
   milp_block->add_static_constraint( *con.release() );
  }

 // eq.(27b)
 for( int j = 0 ; j < n ; ++j ) {
  LinearFunction::v_coeff_pair terms;
  terms.emplace_back( &x_vars[j*n+j] ,  1.0 );
  terms.emplace_back( &u_vars[j]     , -1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( 0.0 ); con->set_rhs( 0.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(28a)
 for( int i = 0 ; i < n ; ++i ) {
  LinearFunction::v_coeff_pair terms;
  for( int j = 0 ; j < n ; ++j )
   terms.emplace_back( &x_vars[i*n+j] , 1.0 );
  auto lf  = std::make_unique< LinearFunction >( std::move(terms) );
  auto con = std::make_unique< FRowConstraint >();
  con->set_lhs( 1.0 ); con->set_rhs( 1.0 );
  con->set_function( lf.release() , eNoMod );
  milp_block->add_static_constraint( *con.release() );
 }

 // eq.(28b)
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

 // Use dedicated Step 2 config if provided, otherwise fall back to Step 1 config
 BlockSolverConfig * cfg_to_use = f_milp_config_step2
   ? f_milp_config_step2 : f_milp_config;

 BlockSolverConfig * milp_cfg =
   static_cast< BlockSolverConfig * >( cfg_to_use->clone() );
 milp_cfg->apply( milp_block.get() );

 Solver * milp_solver = milp_block->get_registered_solvers().empty()
   ? nullptr : milp_block->get_registered_solvers().front();
 if( ! milp_solver ) {
  milp_cfg->clear(); delete milp_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::solve_cssc_milp: no solver attached." );
 }

 const int status = milp_solver->compute();
 if( status != Solver::kOK && status != Solver::kLowPrecision ) {
  milp_cfg->clear(); delete milp_cfg;
  throw std::runtime_error(
    "CSSCScenarioReductionSolver::solve_cssc_milp: MILP failed (status " +
    std::to_string(status) + ")" );
 }

 milp_solver->get_var_solution();

 BLOG( 0 , "\n  CSSC MILP solved, discrepancy = "
           << std::fixed << std::setprecision(6)
           << milp_solver->get_var_value() << std::endl );

 milp_cfg->clear();
 delete milp_cfg;

 // Extract u_j -> representative indices
 ind_red.clear();
 for( int j = 0 ; j < n ; ++j )
  if( u_vars[j].get_value() > 0.5 )
   ind_red.push_back( static_cast< Index >(j) );

 // Numerical safety: if MILP returned wrong count, take top-K by u value
 if( static_cast<int>(ind_red.size()) != K ) {
  std::vector<int> all(n);
  std::iota( all.begin() , all.end() , 0 );
  std::sort( all.begin() , all.end() , [&](int a, int b){
   return u_vars[a].get_value() > u_vars[b].get_value();
  });
  ind_red.clear();
  for( int k = 0 ; k < K ; ++k )
   ind_red.push_back( static_cast< Index >(all[k]) );
  std::sort( ind_red.begin() , ind_red.end() );
 }

 // Extract x_ij -> cluster assignment before milp_block is destroyed.
 // f_scenario_assignment[i] = index of the representative that scenario i
 // is assigned to, as determined by the Step 2 MILP.
 f_scenario_assignment.assign( n , ind_red[0] );
 for( int i = 0 ; i < n ; ++i )
  for( int j = 0 ; j < n ; ++j )
   if( x_vars[i*n+j].get_value() > 0.5 ) {
    f_scenario_assignment[ i ] = static_cast< Index >( j );
    break;
   }

 std::unordered_set< Index > rep_set( ind_red.begin() , ind_red.end() );
 indices_to_choose.clear();
 for( int i = 0 ; i < n ; ++i )
  if( ! rep_set.count( static_cast< Index >(i) ) )
   indices_to_choose.push_back( static_cast< Index >(i) );

 update_reduced_atoms();

 // Wasserstein distance
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