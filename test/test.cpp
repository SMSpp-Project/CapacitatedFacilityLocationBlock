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
 * the maximum number of open facilities, built as a row and not. The rows
 * are relaxed with eNoMod, the Block not handling the Modification of that.
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

#include "FRowConstraint.h"

#include <iostream>

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

int main( void )
{
 const std::vector< std::pair< std::string , void (*)( void ) > > tests = {
  { "natural" , test_natural } ,
  { "knapsack" , test_knapsack } ,
  { "max_facilities" , test_max_facilities } };

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
