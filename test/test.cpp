/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of the feasibility check of the CapacitatedFacilityLocationBlock,
 * on an instance built in memory and with no Solver attached.
 *
 * A point set by hand in the Variable is checked by is_feasible() using the
 * abstract and the physical representation, in the natural and in the
 * knapsack formulations: a violated row makes the point infeasible, and it
 * stops doing so in the abstract representation once it is relaxed, while
 * the physical one, which reads the data, is not affected by relaxing it.
 * The tolerance is checked in its relative and absolute forms, and so is
 * the maximum number of open facilities, built as a row and not. In the
 * "Benders friendly" formulation, which has no X, the open capacity is
 * checked against the demand, and the Solution holds y alone. The rows are
 * relaxed with eNoMod, the Block not handling the Modification of that.
 * Each physical Modification is finally checked for what it says it
 * changes [see Modification::changes()], against what the change does to
 * the model, written down in the tester.
 *
 * The checks hold under NDEBUG as well, not being assert(). The exit code
 * is the number of failed checks.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CapacitatedFacilityLocationBlock.h"

#include "FakeSolver.h"

#include "FRowConstraint.h"

#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using CFLB = CapacitatedFacilityLocationBlock;
using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;  // number of failed checks

/* The instance: 2 facilities of capacity 6 and 10, 2 customers of demand 4
 * each; at most k facilities can be opened. */

static CFLB * make_instance( Index k = CFLB::iInf )
{
 auto b = new CFLB();
 CFLB::CMatrix C( boost::extents[ 2 ][ 2 ] );
 C[ 0 ][ 0 ] = 1; C[ 0 ][ 1 ] = 2;
 C[ 1 ][ 0 ] = 2; C[ 1 ][ 1 ] = 1;
 b->load( 2 , 2 , CFLB::DVector{ 6 , 10 } , CFLB::CVector{ 3 , 3 } ,
	  CFLB::DVector{ 4 , 4 } , std::move( C ) , false , k );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/
/// writes y and x (x[ i ][ j ] row-wise) in the Variable

static void set_point( CFLB * b , std::vector< double > y ,
		       std::vector< double > x )
{
 for( Index i = 0 ; i < 2 ; ++i ) {
  b->get_y( i )->set_value( y[ i ] );
  for( Index j = 0 ; j < 2 ; ++j )
   b->get_x( i , j )->set_value( x[ 2 * i + j ] );
  }
 }

/*--------------------------------------------------------------------------*/
/// the abstract and the physical check, with the given Configuration

static bool abs_ok( CFLB * b , Configuration * c = nullptr )
{
 return( b->is_feasible( true , c ) );
 }

static bool phys_ok( CFLB * b , Configuration * c = nullptr )
{
 return( b->is_feasible( false , c ) );
 }

/*--------------------------------------------------------------------------*/

static void test_natural( void )
{
 auto b = make_instance();
 b->generate_abstract_variables();
 b->generate_abstract_constraints();

 // both facilities open, each customer served by the nearest one
 set_point( b , { 1 , 1 } , { 1 , 0 , 0 , 1 } );
 check( abs_ok( b ) && phys_ok( b ) , "natural: a feasible point fails" );

 // facility 0 serves both customers: 8 > 6
 set_point( b , { 1 , 0 } , { 1 , 1 , 0 , 0 } );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"natural: capacity violated, the point passes" );

 auto cap = b->get_static_constraint_v< FRowConstraint >( "cap" );
 check( cap && ( cap->size() == 2 ) , "natural: no capacity rows" );
 if( ! cap )
  return;

 ( *cap )[ 0 ].relax( true , eNoMod );
 check( abs_ok( b ) , "natural: the violated capacity row is relaxed, yet "
	"the abstract check fails" );
 check( ! phys_ok( b ) , "natural: the physical check reads the relaxed "
	"row" );

 // relaxing a row that is satisfied does not hide the one violated
 ( *cap )[ 0 ].relax( false , eNoMod );
 ( *cap )[ 1 ].relax( true , eNoMod );
 check( ! abs_ok( b ) , "natural: a violated row that is not relaxed "
	"passes" );
 ( *cap )[ 1 ].relax( false , eNoMod );

 // customer 1 half served: the satisfaction row of customer 1 is violated
 set_point( b , { 1 , 1 } , { 1 , 0.5 , 0 , 0 } );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"natural: satisfaction violated, the point passes" );

 auto sat = b->get_static_constraint_v< FRowConstraint >( "sat" );
 check( sat && ( sat->size() == 2 ) , "natural: no satisfaction rows" );
 if( sat ) {
  ( *sat )[ 1 ].relax( true , eNoMod );
  check( abs_ok( b ) , "natural: the violated satisfaction row is relaxed, "
	 "yet the abstract check fails" );
  check( ! phys_ok( b ) , "natural: the physical check reads the relaxed "
	 "satisfaction row" );
  ( *sat )[ 1 ].relax( false , eNoMod );
  }

 // facility 0 loaded with 6.4 > 6: a violation of 0.4
 set_point( b , { 1 , 1 } , { 1 , 0.6 , 0 , 0.4 } );
 SimpleConfiguration< double > tight( 1e-6 );
 check( ( ! abs_ok( b , & tight ) ) && ( ! phys_ok( b , & tight ) ) ,
	"natural: a violation of 0.4 within 1e-6" );

 SimpleConfiguration< std::pair< double , int > > abs05( { 0.5 , 0 } );
 check( abs_ok( b , & abs05 ) && phys_ok( b , & abs05 ) ,
	"natural: a violation of 0.4 beyond an absolute 0.5" );

 SimpleConfiguration< std::pair< double , int > > abs01( { 0.1 , 0 } );
 check( ( ! abs_ok( b , & abs01 ) ) && ( ! phys_ok( b , & abs01 ) ) ,
	"natural: a violation of 0.4 within an absolute 0.1" );

 delete b;
 }

/*--------------------------------------------------------------------------*/

static void test_knapsack( void )
{
 auto b = make_instance();
 SimpleConfiguration< int > ksk( 1 );
 b->generate_abstract_variables( & ksk );
 b->generate_abstract_constraints();

 set_point( b , { 1 , 1 } , { 1 , 0 , 0 , 1 } );
 check( abs_ok( b ) && phys_ok( b ) , "knapsack: a feasible point fails" );

 set_point( b , { 1 , 0 } , { 1 , 1 , 0 , 0 } );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"knapsack: capacity violated, the point passes" );

 // the capacity row of facility 0 is the knapsack one of its sub-Block
 auto ksk0 = b->get_nested_Block( 0 )->
                                 get_static_constraint< FRowConstraint >( 0 );
 check( ksk0 , "knapsack: no knapsack row" );
 if( ksk0 ) {
  ksk0->relax( true , eNoMod );
  check( abs_ok( b ) , "knapsack: the violated knapsack row is relaxed, yet "
	 "the abstract check fails" );
  check( ! phys_ok( b ) , "knapsack: the physical check reads the relaxed "
	 "row" );
  ksk0->relax( false , eNoMod );
  check( ! abs_ok( b ) , "knapsack: the knapsack row is no longer relaxed, "
	 "yet the point passes" );
  }

 delete b;
 }

/*--------------------------------------------------------------------------*/

static void test_max_facilities( void )
{
 // at most one facility: the row is not built, the data are read
 auto b = make_instance( 1 );
 b->generate_abstract_variables();
 b->generate_abstract_constraints();

 set_point( b , { 1 , 1 } , { 1 , 0 , 0 , 1 } );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"max facilities, no row: two open facilities pass" );
 delete b;

 // the row is built
 b = make_instance( 1 );
 b->generate_abstract_variables();
 SimpleConfiguration< int > wc( 7 );
 b->generate_abstract_constraints( & wc );

 set_point( b , { 1 , 1 } , { 1 , 0 , 0 , 1 } );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"max facilities, row: two open facilities pass" );

 auto maxF = b->get_static_constraint< FRowConstraint >( "maxF" );
 check( maxF , "max facilities: no row" );
 if( maxF ) {
  maxF->relax( true , eNoMod );
  check( abs_ok( b ) , "max facilities: the row is relaxed, yet the "
	 "abstract check fails" );
  check( ! phys_ok( b ) , "max facilities: the physical check reads the "
	 "relaxed row" );
  maxF->relax( false , eNoMod );
  }

 // a single facility open, able to serve both customers
 set_point( b , { 0 , 1 } , { 0 , 0 , 1 , 1 } );
 check( abs_ok( b ) && phys_ok( b ) ,
	"max facilities: a feasible point fails" );

 delete b;
 }

/*--------------------------------------------------------------------------*/

static void test_benders( void )
{
 // the "Benders friendly" formulation has no X: only y and v are set, and
 // the rows that need the hidden sub-problem (hence a Solver) are not built
 auto b = make_instance( 1 );
 SimpleConfiguration< int > ben( 2 );
 b->generate_abstract_variables( & ben );
 SimpleConfiguration< int > wc( 4 );  // the maximum facilities row only
 b->generate_abstract_constraints( & wc );

 check( b->get_x( 0 , 0 ) == nullptr , "benders: an X is there" );
 check( b->get_v() != nullptr , "benders: no v" );

 // facility 1 alone covers the demand of 8
 b->get_y( 0 )->set_value( 0 );
 b->get_y( 1 )->set_value( 1 );
 check( abs_ok( b ) && phys_ok( b ) , "benders: a feasible point fails" );

 // facility 0 alone does not: 6 < 8
 b->get_y( 0 )->set_value( 1 );
 b->get_y( 1 )->set_value( 0 );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"benders: open capacity short of the demand, the point passes" );

 SimpleConfiguration< std::pair< double , int > > abs3( { 3 , 0 } );
 check( abs_ok( b , & abs3 ) && phys_ok( b , & abs3 ) ,
	"benders: a shortage of 2 beyond an absolute 3" );

 // both open: the maximum facilities row is violated
 b->get_y( 1 )->set_value( 1 );
 check( ( ! abs_ok( b ) ) && ( ! phys_ok( b ) ) ,
	"benders: two open facilities pass" );

 auto maxF = b->get_static_constraint< FRowConstraint >( "maxF" );
 check( maxF , "benders: no maximum facilities row" );
 if( maxF ) {
  maxF->relax( true , eNoMod );
  check( abs_ok( b ) , "benders: the row is relaxed, yet the abstract "
	 "check fails" );
  check( ! phys_ok( b ) , "benders: the physical check reads the relaxed "
	 "row" );
  maxF->relax( false , eNoMod );
  }

 // the Solution holds y only
 b->get_y( 1 )->set_value( 0 );
 auto sol = b->get_Solution( nullptr , false );
 b->get_y( 0 )->set_value( 0 );
 sol->write( b );
 check( b->get_y( 0 )->get_value() == 1 , "benders: the Solution loses y" );
 delete sol;

 delete b;
 }

/*--------------------------------------------------------------------------*/
/// what each CapacitatedFacilityLocationBlockMod says it changes [see
/// Modification::changes()], against what the change does to the model,
/// written down here: the costs are the Objective, in either direction,
/// the capacities, the demands and the maximum number of facilities the
/// Constraint, closing a facility or fixing it open shrinks the feasible
/// region and re-opening it enlarges it, and the splittability is the
/// integrality of the Variable, unsplittable shrinking the region

static void test_changes( void )
{
 using M = Modification;
 using CM = CapacitatedFacilityLocationBlockMod;
 auto b = make_instance( 2 );
 auto fs = new FakeSolver();
 b->register_Solver( fs );

 // the call has to issue exactly one CapacitatedFacilityLocationBlockMod,
 // of type t, saying that it changes c
 const auto one = [ & ]( const std::function< void( void ) > & call , int t ,
			 M::ModConcern c , const std::string & what ) {
  fs->lock_Modification_list();
  fs->get_Modification_list().clear();
  fs->unlock_Modification_list();
  call();
  std::vector< std::shared_ptr< CM > > got;
  fs->lock_Modification_list();
  for( auto & mod : fs->get_Modification_list() )
   if( auto cm = std::dynamic_pointer_cast< CM >( mod ) )
    got.push_back( cm );
  fs->get_Modification_list().clear();
  fs->unlock_Modification_list();
  check( got.size() == 1 , "changes, " + what + ": " +
	 std::to_string( got.size() ) + " Modification instead of 1" );
  if( got.size() != 1 )
   return;
  check( got.front()->type() == t , "changes, " + what + ": type " +
	 std::to_string( got.front()->type() ) + " instead of " +
	 std::to_string( t ) );
  check( got.front()->changes() == c , "changes, " + what +
	 ": it says it changes " + std::to_string( got.front()->changes() ) +
	 " instead of " + std::to_string( c ) );
  };

 const auto obj = M::eModPhys | M::eModObj | M::eObjUp | M::eObjDown;
 const auto cns = M::eModPhys | M::eModCnsSide | M::eModCnsCoef |
                  M::eRegnShrink | M::eRegnGrow;
 one( [ & ]() { b->chg_facility_cost( 5 , 0 ); } , CM::eChgFCost , obj ,
      "a facility cost" );
 one( [ & ]() { b->chg_transportation_cost( 4 , 1 ); } , CM::eChgTCost ,
      obj , "a transportation cost" );
 one( [ & ]() { b->chg_facility_capacity( 7 , 0 ); } , CM::eChgCap , cns ,
      "a capacity" );
 one( [ & ]() { b->chg_customer_demand( 3 , 1 ); } , CM::eChgDem , cns ,
      "a demand" );
 one( [ & ]() { b->chg_max_facilities( 1 ); } , CM::eChgMxF , cns ,
      "the maximum number of facilities" );
 one( [ & ]() { b->close_facility( 0 ); } , CM::eCloseF ,
      M::eModPhys | M::eModVarData | M::eRegnShrink , "a closing" );
 one( [ & ]() { b->open_facility( 0 ); } , CM::eOpenF ,
      M::eModPhys | M::eModVarData | M::eRegnGrow , "a re-opening" );
 one( [ & ]() { b->fix_open_facility( 1 ); } , CM::eBuyF ,
      M::eModPhys | M::eModVarData | M::eRegnShrink , "a fixing open" );
 one( [ & ]() { b->chg_UnSplittable( true ); } , CM::eChgUnSplt ,
      M::eModPhys | M::eModVarType | M::eRegnShrink , "unsplittable" );
 one( [ & ]() { b->chg_UnSplittable( false ); } , CM::eChgSplt ,
      M::eModPhys | M::eModVarType | M::eRegnGrow , "splittable" );

 b->unregister_Solver( fs );
 delete fs;
 delete b;
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::vector< std::pair< std::string , void (*)( void ) > > tests = {
  { "natural" , test_natural } ,
  { "knapsack" , test_knapsack } ,
  { "max_facilities" , test_max_facilities } ,
  { "benders" , test_benders } ,
  { "changes" , test_changes } };

 for( auto & [ name , test ] : tests )
  try {
   test();
   }
  catch( std::exception & e ) {
   check( false , name + ": exception " + e.what() );
   }

 if( failures )
  std::cout << failures << " checks failed" << std::endl;
 else
  std::cout << "All tests passed!!" << std::endl;

 return( failures );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
