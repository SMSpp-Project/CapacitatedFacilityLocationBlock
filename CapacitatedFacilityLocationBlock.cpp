/*--------------------------------------------------------------------------*/
/*------------- File CapacitatedFacilityLocationBlock.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CapacitatedFacilityLocationBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CapacitatedFacilityLocationBlock.h"

#include "MCFBlock.h"

#include "BinaryKnapsackBlock.h"

#include "AbstractBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * astract and the physical representations agree. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*
using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;
*/

using v_coeff_pair = LinearFunction::v_coeff_pair;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr unsigned char FormMsk = ~3;
// mask for removing the first four bits and only leaving the formulation
// (irrespective of if it is splittable or not)

static constexpr unsigned char StdForm = 0;
// the "standard" formulation is used

static constexpr unsigned char KskForm = 1;
// the "knapsack" formulation is used

static constexpr unsigned char FlwForm = 2;
// the "flow" formulation is used

static constexpr unsigned char UnSpltF = 4;
// fourth bit of AR == 1 if the problem is unsplittable (the X[] are integer)

static constexpr unsigned char HasVar = 8;
// 4th bit of AR == 1 if the Variable have been constructed

static constexpr unsigned char HasObj = 16;
// 5th bit of AR == 1 if the Objective has been constructed

static constexpr unsigned char HasSatCns = 32;
// 6th bit of AR == 1 if the customer satisfaction Constraints are constructed

static constexpr unsigned char HasCapCns = 64;
// 7th bit of AR == 1 if the capacity Constraints are constructed

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static BinaryKnapsackBlock * BKB( Block * b ) {
 return( static_cast< BinaryKnapsackBlock * >( b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static AbstractBlock * AB( Block * b ) {
 return( static_cast< AbstractBlock * >( b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static MCFBlock * MCFB( Block * b ) {
 return( static_cast< MCFBlock * >( b ) );
 }

/*----------------------------------------------------------------------------

// returns the number of elements where two vectors differ

template< typename T >
static Index countdiff( T beg , T end , T cmp )
{
 Index ndiff = 0;
 for( ; beg != end ; )
  if( *(beg++) != *(cmp++) )
   ndiff++;

 return( ndiff );
 }

------------------------------------------------------------------------------
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector<T> & vec , c_Subset & nms ,
		      typename std::vector<T>::const_iterator cmp ,
		      Index n_max )
{
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   return( false );
  }

 return( true );
 }

------------------------------------------------------------------------------
// returns the number of elements where two vectors differ, one of them
// being given as a base vector and a subset of indices

template< typename T >
static Index countdiff( std::vector<T> & vec , c_Subset & nms ,
			typename std::vector<T>::const_iterator cmp ,
			Index n_max )
{
 Index ndiff = 0;
 for( auto nm : nms ) {
  if( nm >= n_max )
   throw( std::invalid_argument( "invalid name in nms" ) );
  if( vec[ nm ] != *(cmp++) )
   ndiff++;
  }

 return( ndiff );
 }

------------------------------------------------------------------------------
// copys one vector to a given subset of another

template< typename T >
static void copyidx( std::vector<T> & vec , c_Subset & nms ,
		     typename std::vector<T>::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

----------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register CapacitatedFacilityLocationBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( CapacitatedFacilityLocationBlock );

/*--------------------------------------------------------------------------*/
/*--------------- METHODS OF CapacitatedFacilityLocationBlock --------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::load( Index n , Index m ,
					     DVector && Q , CVector && F ,
					     DVector && D , CMatrix && C )
{
 static const std::string _prfx = "CapacitatedFacilityLocationBlock::load: ";

 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( m == 0 )
  throw( std::invalid_argument( _prfx + "number of facilities too small" ) );

 if( n == 0 )
  throw( std::invalid_argument( _prfx + "number of customers too small" ) );

 if( Q.size() != m )
  throw( std::invalid_argument( _prfx + "capacity vector has wrong size" ) );

 if( std::any_of( Q.begin() , Q.begin() + m ,
		  []( auto qi ) { return( qi <= 0 ); } ) )
  throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

 if( F.size() != m )
  throw( std::invalid_argument( _prfx + "opening cost vector has wrong size"
				) );
 if( D.size() != n )
  throw( std::invalid_argument( _prfx + "demands vector has wrong size" ) );

 if( std::any_of( D.begin() , D.begin() + n ,
		  []( auto di ) { return( di <= 0 ); } ) )
  throw( std::invalid_argument( _prfx + "non-positive demand" ) );

 auto shp = C.shape();
 if( ( shp[ 0 ] != m ) || ( shp[ 1 ] != n ) )
  throw( std::invalid_argument( _prfx +
			   "transportation cost matrix has wrong shape"	) );

 // erase existing abstract representation, if any - - - - - - - - - - - - - -

 if( AR & ~7 )
  guts_of_destructor();
		   
 // move over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 f_n_facilities = m;
 f_n_customers = n;

 v_capacity = std::move( Q );
 v_fixed_cost = std::move( F );
 v_demand = std::move( D );
 v_transp_cost = std::move( C );

 f_cond_lower = f_cond_upper = dNaN;  // reset conditional bounds

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( CapacitatedFacilityLocationBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::load( std::istream & input )
{
 static const std::string _prfx = "CapacitatedFacilityLocationBlock::load: ";

 // erase existing abstract representation, if any - - - - - - - - - - - - - -

 if( AR & ~7 )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -

 input >> eatcomments >> f_n_facilities;
 if( input.fail() )
  goto( input_failure );

 if( f_n_facilities == 0 )
  throw( std::invalid_argument( _prfx + "number of facilities too small" ) );

 input >> eatcomments >> f_n_customers;
 if( input.fail() )
  goto( input_failure );

 if( f_n_customers == 0 )
  throw( std::invalid_argument( _prfx + "number of customers too small" ) );

 v_capacity.resize( f_n_facilities );
 v_fixed_cost.resize( f_n_facilities );

 for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
  input >> eatcomments >> v_capacity[ i ];
  if( input.fail() )
   goto( input_failure );
  if( v_capacity[ i ] <= 0 )
   throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

  input >> eatcomments >> v_fixed_cost[ i ];
  if( input.fail() )
   goto( input_failure );

  }  // end( for( each facility ) )

 v_demand.resize( f_n_customers );
 v_transp_cost.resize( boost::extents[ f_n_facilities ][ f_n_customers ] );

 for( Index j = 0 ; j < f_n_customers ; ++j ) {  // for( each customer )
  input >> eatcomments >> v_demand[ j ];
  if( input.fail() )
   goto( input_failure );
  if( v_demand[ j ] <= 0 )
   throw( std::invalid_argument( _prfx + "non-positive demand" ) );

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
   input >> eatcomments >> v_transp_cost[ i ][ j ];
   if( input.fail() )
    goto( input_failure );
  
   }  // end( for( each facility ) )
  }  // end( for( each customer ) )

 f_cond_lower = f_cond_upper = dNaN;  // reset conditional bounds

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 return;

 input_failure:

 throw( std::logic_error( _prfx + "error reading from stream" ) );

 }  // end( CapacitatedFacilityLocationBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::deserialize(
					      const netCDF::NcGroup & group )
{
 static const std::string _prfx =
                            "CapacitatedFacilityLocationBlock::deserialize: ";

 // erase existing abstract representation, if any - - - - - - - - - - - - - -

 if( AR & ~7 )
  guts_of_destructor();
		   
 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto nf = group.getDim( "NFacilities" );
 if( nf.isNull() )
  throw( std::logic_error( _prfx + "NFacilities dimension is required" ) );
 f_n_facilities = nf.getSize();
 if( f_n_facilities == 0 )
  throw( std::invalid_argument( _prfx + "number of facilities too small" ) );

 auto nc = group.getDim( "NCustomers" );
 if( nc.isNull() )
  throw( std::logic_error( _prfx + "NCustomers dimension is required" ) );
 f_n_customers = nc.getSize();
 if( f_n_customers == 0 )
  throw( std::invalid_argument( _prfx + "number of customers too small" ) );

 auto fq = group.getVar( "FacilityCapacity" );
 if( fq.isNull() )
  throw( std::logic_error( _prfx + "FacilityCapacity not found" ) );
 auto fqs = ::get_sizes_dimensions( fq );
 if( ( fqs.size() != 1 ) || ( fqs[ 0 ] != f_n_facilities ) )
  throw( std::logic_error( _prfx + "FacilityCapacity has wrong size" ) );
 
 v_capacity.resize( f_n_facilities );
 fq.getVar( v_capacity.data() );
 if( std::any_of( v_capacity.begin() , v_capacity.begin() + f_n_facilities ,
		  []( auto qi ) { return( qi <= 0 ); } ) )
  throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

 auto fc = group.getVar( "FacilityCost" );
 if( fc.isNull() )
  throw( std::logic_error( _prfx + "FacilityCost not found" ) );
 auto fcs = ::get_sizes_dimensions( fc );
 if( ( fcs.size() != 1 ) || ( fcs[ 0 ] != f_n_facilities ) )
  throw( std::logic_error( _prfx + "FacilityCost has wrong size" ) );
 
 v_fixed_cost.resize( f_n_facilities );
 fq.getVar( v_fixed_cost.data() );

 auto cd = group.getVar( "CustomerDemand" );
 if( cd.isNull() )
  throw( std::logic_error( _prfx + "CustomerDemand not found" ) );
 auto cds = ::get_sizes_dimensions( cd );
 if( ( cds.size() != 1 ) || ( cds[ 0 ] != f_n_customers ) )
  throw( std::logic_error( _prfx + "CustomerDemand has wrong size" ) );
 
 v_demand.resize( f_n_customers );
 cd.getVar( v_demand.data() );
 if( std::any_of( v_demand.begin() , v_demand.begin() + f_n_customers ,
		  []( auto dj ) { return( dj <= 0 ); } ) )
  throw( std::invalid_argument( _prfx + "non-positive demand" ) );

 auto tc = group.getVar( "TransportationCost" );
 if( tc.isNull() )
  throw( std::logic_error( _prfx + "TransportationCost not found" ) );
 auto tcs = ::get_sizes_dimensions( tc );
 if( ( tcs.size() != 2 ) ||
     ( tcs[ 0 ] != f_n_facilities ) || ( tcs[ 1 ] != f_n_customers ) )
  throw( std::logic_error( _prfx + "TransportationCost has wrong size" ) );

 v_transp_cost.resize( tcs );
 tc.getVar( std::vector< std::size_t >( 2 , 0 ) , tcs ,
	    v_transp_cost.data() );

 f_cond_lower = f_cond_upper = dNaN;  // reset conditional bounds

 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

 }  // end( CapacitatedFacilityLocationBlock::deserialize )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_abstract_variables(
						       Configuration * stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do

 AR |= HasVar;      // variables will be constructed now once and for all
 f_unsplittable = wf & 4;
 
 Index wf = 0;
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  wf = sci->f_value;
 
 if( ! ( wf & 3 ) ) {  // "natural formulation" (NF)- - - - - - - - - - - - -
                       // - - - - - - - - - - - - - - - - - - - - - - - - - -
  // AR |= StdForm;  does nothing
  v_y.resize( f_n_facilities );
  for( auto yi : v_y )
   yi.set_type( ColVariable::kBinary );
  add_static_variable( v_y , "y" );

  v_x.resize( { f_n_facilities , f_n_customers } );
  auto xt = ColVariable::kPosUnitary;
  if( f_unsplittable ) {
   AR |= UnSpltF;
   xt = ColVariable::kBinary;
   }

  for( auto xji : v_x )
   xji.set_type( xt );
  add_static_variable( v_x , "x" );

  return:
  }

 if( ( wf & 3 ) == 1 ) {  // "knapasck formulation" (KF)- - - - - - - - - - -
                          //- - - - - - - - - - - - - - - - - - - - - - - - -
  AR |= KskForm;
  // construct one knapsack problem for each facility
  v_Block.resize( f_n_facilities );

  BinaryKnapsackBlock::doubleVec W( f_n_customers + 1 );
  BinaryKnapsackBlock::doubleVec C( f_n_customers + 1 );
  BinaryKnapsackBlock::boolVec I;

  if( f_unsplittable ) {
   AR |= UnSpltF;
   I.resize( f_n_customers + 1 , true );
   }
  else {
   I.resize( f_n_customers + 1 , false );
   I[ f_n_customers ] = true;
   }

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   for( Index j = 0 ; j < f_n_customers ; ++j ) {
    W[ j ] = v_demand[ j ];
    C[ j ] = v_transp_cost[ i ][ j ];
    }
   W[ f_n_customers ] = - v_capacity[ i ];
   C[ f_n_customers ] = v_fixed_cost[ i ];

   auto ki = new BinaryKnapsackBlock( this );
   ki->load( f_n_customers + 1 , 0 , W , P , I );
   ki->generate_abstract_variables();
   v_Block[ i ] = ki;
   }

  return:
  }

 if( ( wf & 3 ) == 2 ) {  // "flow formulation" (FF)- - - - - - - - - - - - -
                          //- - - - - - - - - - - - - - - - - - - - - - - - -
  AR |= FlwForm;

  v_Block.resize( 2 );  // exactly two sub-Block

  // the first sub-Block is an AbstractBlock with the y[] variables
  auto ab = new AbstractBlock( this );
  v_Block[ 0 ] = ab;

  v_y.resize( f_n_facilities );
  for( auto yi : v_y )
   yi.set_type( ColVariable::kBinary );
  ab->add_static_variable( v_y , "y" );

  // the second Block is a MCFBlock as constructed by get_R3_Block
  static SimpleConfiguration< int > r3bc( 2 );
  auto mcfb = static_cast< MCFBlock * >( get_R3_Block( & r3bc ) );
  v_Block[ 1 ] = mcfb;

  // ... except the cost of the facility arcs are zerod
  MCFBlock::Vec_CNumber zero( f_n_facilities , 0 );
  mcfb->chg_costs( zero.begin() , Range( 0 , f_n_facilities ),
		   eNoMod , eNoMod );
  mcfb->generate_abstract_variables();
  mcfb->set_f_Block( this );

  return:
  }

 throw( std::invalid_argument(
	   "CapacitatedFacilityLocationBlock::generate_abstract_variables: "
	   "invalid formulation" ) );

 }  // end( CapacitatedFacilityLocationBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_abstract_constraints(
						       Configuration * stcc )
{
 const static std::string _prfx =
         "CapacitatedFacilityLocationBlock::generate_abstract_constraints: ";

 if( ! ( AR & HasVar ) )
  throw( std::logic_error( _prfx + "generate_abstract_variables not called"
			   ) );
 Index wc = 0;
 if( ( ! stcc ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  wc = sci->f_value;

 if( ( AR & FormMsk ) == StdForm ) {  // "natural formulation" (NF) - - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  if( ( ! ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) ) {
   // construct customer satisfaction constraints (not there already)
   v_sat.resize( f_n_customers );
   for( Index j = 0 ; j < f_n_customers ; ++j ) {
    v_coeff_pair coeffs( f_n_facilities );

    for( Index i = 0 ; i < f_n_facilities ; ++i )
     coeffs[ i ] = std::make_pair( & v_x[ i ][ j ] , double( 1 ) );

    v_sat[ j ].set_both( 1 );
    v_sat[ j ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
    }

   add_static_constraint( v_sat , "sat" );
   AR |= HasSatCns;
   }

  if( ( ! ( wc & 2 ) ) && ( ! ( AR & HasCapCns ) ) ) {
   // construct facility capacity constraints (not there already)
   v_cap.resize( f_n_facilities );
   for( Index i = 0 ; i < f_n_facilities ; ++i ) {
    v_coeff_pair coeffs( f_n_customers + 1 );

    for( Index j = 0 ; j < f_n_customers ; ++j )
     coeffs[ j ] = std::make_pair( & v_x[ i ][ j ] , v_demand[ j ] );

    coeffs[ f_n_customers ] = std::make_pair( & v_y[ i ] ,
					      - v_capacity[ i ] );
    v_cap[ i ].set_rhs( 0 );
    v_cap[ i ].set_lhs( 1 ) = - Inf< RHSValue >();
    v_cap[ i ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
    }

   add_static_constraint( v_cap , "cap" );
   AR |= HasCapCns;
   }

  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // "knapsack formulation" (KF)- - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  if( ( ! ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) ) {
   // construct customer satisfaction constraints (not there already)
   v_sat.resize( f_n_customers );
   for( Index j = 0 ; j < f_n_customers ; ++j ) {
    v_coeff_pair coeffs( f_n_facilities );

    for( Index i = 0 ; i < f_n_facilities ; ++i )
     coeffs[ i ] = std::make_pair( BKB( v_Block[ i ] )->get_Var( j ) ,
				   double( 1 ) );
    v_sat[ j ].set_both( 1 );
    v_sat[ j ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
    }

   add_static_constraint( v_sat , "sat" );
   AR |= HasSatCns;
   }

  if( ( ! ( wc & 2 ) ) && ( ! ( AR & HasCapCns ) ) ) {
   // construct facility capacity constraints (not there already)
   // these are just the constraint in the BinaryKnapsackBlock
   for( auto ki : v_Block )
    ki->generate_abstract_constraints();

   AR |= HasCapCns;
   }

  return;
  }

 // else it is the "flow formulation" (FF)- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) ) {
  // construct customer satisfaction constraints (not there already)
  // these are, just the MCFBlock ones
  v_Block[ 1 ]->generate_abstract_constraints();
  AR |= HasSatCns;
  }

 if( ( ! ( wc & 2 ) ) && ( ! ( AR & HasCapCns ) ) ) {
  // construct facility capacity constraints (not there already)
  // these are the ones linking the y[] variables to the
  // arc flow variables of facility arcs

  auto ab = AB( v_Block[ 0 ] );
  auto mcfb = MCFB( v_Block[ 1 ] );

  v_cap.resize( f_n_facilities );
  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   v_coeff_pair coeffs( 2 );

   coeffs[ 0 ] = std::make_pair( mcfb->i2p_x( i ) , double( 1 ) );
   coeffs[ 0 ] = std::make_pair( & v_y[ i ] , - v_capacity[ i ] );

   v_cap[ i ].set_rhs( 0 );
   v_cap[ i ].set_lhs( 1 ) = - Inf< RHSValue >();
   v_cap[ i ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
   }

  ab->add_static_constraint( v_sat , "sat" );

  AR |= HasCapCns;
  }

 /*!!
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif
 !!*/

 }  // end( CapacitatedFacilityLocationBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_objective(
						       Configuration * objc )
{
 const static std::string _prfx =
                    "CapacitatedFacilityLocationBlock::generate_objective: ";
 
 if( AR & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error( _prfx + "generate_abstract_variables not called"
			   ) );
 
 AR |= HasObj;      // Objective will be constructed now once and for all

 if( ( AR & FormMsk ) == StdForm ) {  // "natural formulation" (NF) - - - - -
                                      //- - - - - - - - - - - - - - - - - - -

  // construct a "dense" LinearFunction
  v_coeff_pair p( f_n_facilities * ( f_n_customers + 1 ) );
  auto pi = p.begin();

  // first the Y[ i ] variables
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   *(pi++) = std::make_pair( & v_y[ i ] , v_fixed_cost[ i ] );

  // then the X[ j ][ i ] ones
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    for( Index j = 0 ; j < f_n_customers ; ++j )
     *(pi++) = std::make_pair( & v_x[ i ][ j ] , v_transp_cost[ i ][ j ] );

  c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
  set_objective( & c , eNoMod );
  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // "knapsack formulation" (KF)- - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  for( auto ki : v_Block )    // the Objective is all in the sub-Block
   ki->generate_objective();

  return;
  }

 // else it is the "flow formulation" (FF)- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // construct a "dense" LinearFunction for the y[] variables only in the
 // first sub-Block

 v_coeff_pair p( f_n_facilities );
 for( Index i = 0 ; i < f_n_facilities ; ++i )
  p[ i ] = std::make_pair( & v_y[ i ] , v_fixed_cost[ i ] );

 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 AB( v_Block[ 0 ] )->set_objective( & c , eNoMod );

 // generate the objective in the MCFBlock
 MCFB( v_Block[ 1 ] )->generate_objective();

 /*!!
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif
 !!*/

 }  // end( CapacitatedFacilityLocationBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::is_feasible( bool useabstract ,
						    Configuration * fsbc )
{
 double eps = 1e-10;

 if( ( ! fsbc ) && f_BlockConfig )
  fsbc = f_BlockConfig->f_is_feasible_Configuration;

 if( auto tfsbc = dynamic_cast< SimpleConfiguration< double > * >( fsbc ) )
  eps = tfsbc->f_value;

 if( ! ( AR & HasVar ) )
  throw( std::logic_error( "CapacitatedFacilityLocationBlock::is_feasible"
			   "generate_abstract_variables not called" ) );

 // check variable feasibility
 for( Index i = 0 ; i < f_n_facilities ; ++i )
  if( ! get_y( i ).is_feasible( eps ) )
   return( false );

 for( Index i = 0 ; i < f_n_facilities ; ++i )
  for( Index j = 0 ; j < f_n_customers ; ++j )
   if( ! get_x( i , j ).is_feasible( eps ) )
    return( false );

 return( customer_feasible( eps , useabstract ) &&
	 facility_feasible( eps , useabstract ) );

 }  //  end( CapacitatedFacilityLocationBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::customer_feasible( double eps ,
							  bool useabstract )
{
 const static std::string _prfx =
                     "CapacitatedFacilityLocationBlock::customer_feasible: ";

 if( ! ( AR & HasSatCns ) )  // if customer constraints are not defined
  useabstract = false;       // you cannot use them to chek feasibility

 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == StdForm ) ||
      ( ( AR & FormMsk ) == KskForm ) ) {
   for( const auto & cnst : v_sat )
    if( cnst.rel_viol() > eps )
     return( false );

   return( true );
   }

  SimpleConfiguration< double > cfg( eps ); 
  return( MCFB( v_Block[ 1 ] )->is_feasible( true , & cfg );
  
  throw( std::logic_error( _prfx + "flow formulation not supported yet" ) );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   auto tot = 0;
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    tot += get_x( i , j ).get_value();

   if( std::abs( 1 - tot ) > eps )
    return( false );
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::customer_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::capacity_feasible( double eps ,
							  bool useabstract )
{
 const static std::string _prfx =
                     "CapacitatedFacilityLocationBlock::capacity_feasible: ";

 if( ! ( AR & HasCapCns ) )  // if capacity constraints are not defined
  useabstract = false;       // you cannot use them to chek feasibility

 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == StdForm ) ||
      ( ( AR & FormMsk ) == FlwForm ) ) {
   for( const auto & cnst : v_cap )
    if( cnst.rel_viol() > eps )
     return( false );

   return( true );
   }

  if( ( AR & FormMsk ) == KskForm ) {
   SimpleConfiguration< double > cfg( eps );

   for( Index i = 0 ; i < f_n_facilities ; ++i )
    if( ! v_Block[ i ]->is_feasible( true , & cfg ) )
     return( false );

   return( true );
   }

  throw( std::logic_error( _prfx + "flow formulation not supported yet" ) );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   auto tot = v_capacity[ i ] * get_y( i ).get_value();
   for( Index j = 0 ; j < f_n_customers ; ++j )
    tot -= v_demand[ j ] * get_x( i , j ).get_value();

   if( tot < - eps * v_capacity[ i ] )
    return( false );
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::capacity_feasible )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

Block * CapacitatedFacilityLocationBlock::get_R3_Block( Configuration * r3bc ,
					      Block * base , Block * father )
{
 const static std::string _prfx =
                           "CapacitatedFacilityLocationBlock::get_R3_Block: ";
 int wR3B = 0
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  CapacitatedFacilityLocationBlock * CFLB;
  if( base ) {
   CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( base );
   if( ! CFLB )
    throw( std::invalid_argument( _prfx +
			"base is not a CapacitatedFacilityLocationBlock" ) );
   }
  else
   CFLB = new CapacitatedFacilityLocationBlock( father );

  CFLB->load( f_n_facilities , f_n_customers , v_capacity , v_fixed_cost ,
	      v_demand , v_transp_cost );
 
  return( CFLB );
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFBlock * MCFB;
 if( base ) {
  MCFB = dynamic_cast< MCFBlock * >( base );
  if( ! MCFB )
   throw( std::invalid_argument( _prfx + "base is not a MCFBlock" ) );
   }
 else
  MCFB = new MCFBlock( father );

 Index NN = f_n_facilities + f_n_customers + 1;
 Index NA = f_n_facilities * ( f_n_customers + 1 )
          + wR3B > 1 ? f_n_customers : 0;

 Subset EN( NA );
 Subset SN( NA );
 MCFBlock::Vec_FNumber U( NA );
 MCFBlock::Vec_CNumber C( NA );
 MCFBlock::Vec_FNumber B( NN );

 // construct deficits vector
 Index i = 0;
 while( i <= f_n_facilities )
  B[ i++ ] = 0;

 for( Index j = 0 ; j < f_n_customers ; j++ ) {
  B[ 0 ] -= v_demand[ j ];
  B[ i++ ] = v_demand[ j ];
  }

 // construct arcs SN, EN, U, C: "common part" of the graph
 Index a = 0;

 // first the source -> facility arcs
 for( i = 0 ; i < f_n_facilities ; ++i ) {
  SN[ a ] = 1;
  EN[ a ] = i + 2;
  U[ a ] = v_capacity[ i ];
  C[ a++ ] = v_fixed_cost[ i ] / v_capacity[ i ];
  }

 // now the facility -> customers arcs
  for( i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j ) {
    SN[ a ] = i + 2;
    EN[ a ] = f_n_facilities + 2 + j;
    U[ a ] = Inf< MCFClass::FNumber >();
    C[ a++ ] = v_transp_cost[ i ][ j ] / v_demand[ j ];
    }

 if( wR3B > 1 ) {
  // now the artificial arcs to ensure feasibility

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   SN[ a ] = 1;
   EN[ a ] = f_n_facilities + 2 + j;
   U[ a ] = Inf< MCFClass::FNumber >();

   // compute an upper bound on the worst-case transportation cost
   MCFClass::CNumber maxc = 0;
   for( i = 0 ; i < f_n_facilities ; ++i , ++a )
    if( auto tci = C[ i ] + v_transp_cost[ i ][ j ] / v_demand[ j ] ;
	maxc < tci )
     maxc = tci;
   maxc += 1;    // ! +1
   maxc *= 100;  // ! *100 
   C[ a++ ] = maxc;
   }
  }

 MCFB->load( NN , NA , EN , SN , U , C , B );

 return( MCFB );

 }  // end( CapacitatedFacilityLocationBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::map_back_solution( Block * R3B ,
				 Configuration * r3bc , Configuration * solc )
{
 const static std::string _prfx =
                      "CapacitatedFacilityLocationBlock::map_back_solution: ";

 if( ! ( AR & HasVar ) )
  throw( std::invalid_argument(  _prfx + "variables not generated yet" ) );

 int ws = 3;
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  ws = tcfg->f_value;

 if( ! ( ws & 3 ) )  // actually nothing to map back
  return;            // silently (and cowardly) return

 int wR3B = 0
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( R3B );
  if( ! CFLB )
   throw( std::invalid_argument( _prfx +
			 "R3B is not a CapacitatedFacilityLocationBlock" ) );

  if( ( get_NFacilities() != CFLB->get_NFacilities() ) ||
      ( get_NCustomers() != CFLB->get_NCustomers() ) )
   throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

  if( ws & 1 )
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    get_y( i ).set_value( CFLB->get_y( i ).get_value() );

  if( ws & 2 )
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    for( Index j = 0 ; j < f_n_customers ; ++j )
     get_x( i , j ).set_value( CFLB->get_x( i , j ).get_value() );

  return;
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( _prfx + "R3B is not a MCFBlock" ) );

 if( ( MCFB->get_NNodes() != f_n_facilities + f_n_customers + 1 ) ||
     ( MCFB->get_NArcs() != f_n_facilities * ( f_n_customers + 1 ) +
                            wR3B == 2 ? f_n_customers : 0 ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

 Index l = ws & 1 ? f_n_facilities : 0;
 Index u = ws & 2 ? f_n_facilities * ( f_n_customers + 1 ) : f_n_facilities;
 
 MCFBlock::Vec_FNumber x( u - l );

 MCFB->get_x( x , Range( l , u ) );

 auto it = x.begin();
 if( ws & 1 )
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   get_y( i ).set_value( *(it++) / v_capacity[ i ] );

 if( ws & 2 )
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j )
    get_x( i , j ).set_value( *(it++) / v_demand[ j ] );

 }  // end( CapacitatedFacilityLocationBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::map_forward_solution( Block * R3B ,
			         Configuration * r3bc , Configuration * solc )
{
 const static std::string _prfx =
                   "CapacitatedFacilityLocationBlock::map_forward_solution: ";

 if( ! ( AR & HasVar ) )
  throw( std::invalid_argument(  _prfx + "variables not generated yet" ) );

 int ws = 3;
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  ws = tcfg->f_value;

 if( ! ( ws & 3 ) )  // actually nothing to map forward
  return;            // silently (and cowardly) return

 int wR3B = 0
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( R3B );
  if( ! CFLB )
   throw( std::invalid_argument( _prfx +
			 "R3B is not a CapacitatedFacilityLocationBlock" ) );

  // fantastically dirty trick: because the two objects are copies, mapping
  // forward a solution from this to R3B is the same as mapping back a
  // solution from R3B to this

  CFLB->map_back_solution( this , r3bc , solc );
  return;
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( _prfx + "R3B is not a MCFBlock" ) );

 if( ( MCFB->get_NNodes() != f_n_facilities + f_n_customers + 1 ) ||
     ( MCFB->get_NArcs() != f_n_facilities * ( f_n_customers + 1 ) +
                            wR3B == 2 ? f_n_customers : 0 ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

 Index l = ws & 1 ? f_n_facilities : 0;
 Index u = ws & 2 ? f_n_facilities * ( f_n_customers + 1 ) : f_n_facilities;
 
 MCFBlock::Vec_FNumber x( u - l );

 auto it = x.begin();
 if( ws & 1 )
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   *(it++) = get_y( i ).get_value() * v_capacity[ i ];

 if( ws & 2 )
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j )
    *(it++) = get_x( i , j ).get_value() * v_demand[ j ];

 MCFB->set_x( x , Range( l , u ) );

 }  // end( CapacitatedFacilityLocationBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::map_forward_Modification(
			   Block * R3B , c_p_Mod mod , Configuration * r3bc ,
			   ModParam issuePMod , ModParam issueAMod )
{
 if( mod->concerns_Block() )  // an abstract Modification
  return( false );            // none of my business

 const static std::string _prfx =
              "CapacitatedFacilityLocationBlock::map_forward_Modification: ";

 int wR3B = 0
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( R3B );
  if( ! CFLB )
   throw( std::invalid_argument( _prfx +
			 "R3B is not a CapacitatedFacilityLocationBlock" ) );

  if( ( get_NFacilities() != CFLB->get_NFacilities() ) ||
      ( get_NCustomers() != CFLB->get_NCustomers() ) )
   throw( std::invalid_argument( _prfx + "incompatible sizes in R3Block" ) );

  return( guts_of_map_f_Mod_copy( CFLB , mod , issuePMod , issueAMod ) );
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( _prfx + "R3B is not a MCFBlock" ) );

 if( ( MCFB->get_NNodes() != f_n_facilities + f_n_customers + 1 ) ||
     ( MCFB->get_NArcs() != f_n_facilities * ( f_n_customers + 1 ) +
                            wR3B == 2 ? f_n_customers : 0 ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3Block" ) );

 return( guts_of_map_f_Mod_MCF( MCFB , mod , issuePMod , issueAMod ) );


 }  // end( CapacitatedFacilityLocationBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::map_back_Modification( Block * R3B ,
			            c_p_Mod mod , Configuration *r3bc ,
				    ModParam issuePMod , ModParam issueAMod )
{
 const static std::string _prfx =
                 "CapacitatedFacilityLocationBlock::map_back_Modification: ";
 int wR3B = 0
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( R3B );
  if( ! CFLB )
   throw( std::invalid_argument( _prfx +
			"r3bc is not a CapacitatedFacilityLocationBlock" ) );

  // fantastically dirty trick: because the two objects are copies, mapping
  // back a Modification to this from R3B is the same as mapping forward a
  // Modification from R3B to this

  return( CFLB->map_forward_Modification( this , mod , r3bc , issuePMod ,
					  issueAMod ) );
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( _prfx + "R3B is not a MCFBlock" ) );

 if( ( MCFB->get_NNodes() != f_n_facilities + f_n_customers + 1 ) ||
     ( MCFB->get_NArcs() != f_n_facilities * ( f_n_customers + 1 ) +
                            wR3B == 2 ? f_n_customers : 0 ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

 // TODO:: implement
 // return( guts_of_map_b_Mod_MCF( MCFB , mod , issuePMod , issueAMod ) );

 return( false );  // currently, no modification is properly mapped back

 }  // end( CapacitatedFacilityLocationBlock::map_back_Modification )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * CapacitatedFacilityLocationBlock::get_Solution(
				         Configuration * solc , bool emptys )
{
 auto * sol = new CapacitatedFacilityLocationSolution();

 int wsol = 0;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 if( wsol != 2 )
  sol->v_y.resize( f_n_facilities );

 if( wsol != 1 )
  sol->v_x.resize( { f_n_customers , f_n_facilities } );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( CapacitatedFacilityLocationBlock::get_Solution )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::add_Modification( sp_Mod mod ,
							 ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*---- LOADING, PRINTING & SAVING THE CapacitatedFacilityLocationBlock -----*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::serialize( netCDF::NcGroup & group )
 const
{
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // now the CapacitatedFacilityLocationBlock data - - - - - - - - - - - - - -

 std::vector< netCDF::NcDim > dims( 2 );

 auto nf = group.addDim( "NFacilities" , f_n_facilities );
 auto nc = group.addDim( "NCustomers" , f_n_customers );

 ( group.addVar( "FacilityCapacity" , netCDF::NcUint64() , nf )
   ).putVar( v_capacity.data() );

 ( group.addVar( "FacilityCost" , netCDF::NcUint64() , nf )
   ).putVar( v_fixed_cost.data() );

 ( group.addVar( "CustomerDemand" , netCDF::NcUint64() , nc )
   ).putVar( v_demand.data() );

 ::serialize( group , "TransportationCost" , netCDF::NcDouble() ,
              { nc , nf } , v_transp_cost );

 }  // end( CapacitatedFacilityLocationBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_costs( c_Vec_CNumber_it NCost ,
						  Range rng ,
						  ModParam issueMod ,
						  ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( C.empty() ) {
  if( std::all_of( NCost , NCost + ( rng.second - rng.first ) ,
		   []( c_CNumber cst ) { return( cst == 0 ); } ) )
   return;

  C.assign( get_MaxNArcs() , 0 );
  }

 if( std::equal( NCost , NCost + ( rng.second - rng.first ) ,
		 C.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCost , NCost + ( rng.second - rng.first ) ,
	     C.begin() + rng.first );

  // note that modify_coefficients owns the vector, so a copy has to be made
  get_lfo()->modify_coefficients( Vec_CNumber( NCost , NCost +
					       ( rng.second - rng.first ) ) ,
				  rng , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCost , NCost + ( rng.second - rng.first ) ,
	      C.begin() + rng.first );

 f_cond_lower = NAN;  // reset conditional bounds
 
 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
					      CapacitatedFacilityLocationBlockMod::eChgCost , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_costs( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_costs( c_Vec_CNumber_it NCost , Subset && nms ,
			  const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( C.empty() ) {
  if( std::all_of( NCost , NCost + nms.size() ,
		   []( c_CNumber cst ) { return( cst == 0 ); } ) )
   return;

  C.assign( get_MaxNArcs() , 0 );
  }

 if( is_equal( C , nms , NCost , get_NArcs() ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( C , nms , NCost );

  // note that modify_coefficients owns both vectors, so two copies have
  // to be made
  get_lfo()->modify_coefficients( Vec_CNumber( NCost , NCost + nms.size() ) ,
				  Subset( nms ) , ordered , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( C , nms , NCost );

 f_cond_lower = NAN;  // reset conditional bounds

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockSbstMod>( this ,
                                  CapacitatedFacilityLocationBlockMod::eChgCost , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_costs( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_cost( c_CNumber NCost , c_Index arc , 
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( C.empty() && NCost )
  C.assign( get_MaxNArcs() , 0 );

 if( C[ arc ] == NCost )
  return;

 f_cond_lower = NAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  C[ arc ] = NCost;

  get_lfo()->modify_coefficient( arc , NCost , issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   C[ arc ] = NCost;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			   CapacitatedFacilityLocationBlockMod::eChgCost , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_cost )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_ucaps( c_Vec_FNumber_it NCap , Range rng ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( U.empty() ) {
  if( std::all_of( NCap , NCap + ( rng.second - rng.first ) ,
		   []( c_FNumber cap ) { return( cap >= Inf<FNumber>() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf<FNumber>() );
  }

 c_Index ndiff = countdiff( NCap , NCap + ( rng.second - rng.first ) ,
			    U.cbegin() + rng.first );
 if( ! ndiff )
  return;

 f_cond_lower = NAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ;  ++i , ++NCap )
   if( U[ i ] != *NCap ) {
    U[ i ] = *NCap;
    UB[ i ].set_rhs( *NCap , ampar );
    }

  // dynamic part
  for( auto dubi = std::next( dUB.begin() ,
			      rng.first >= get_NStaticArcs() ?
			      i - get_NStaticArcs() : 0 ) ;
       i < rng.second ; ++i , ++NCap , ++dubi )
   if( U[ i ] != *NCap ) {
    U[ i ] = *NCap;
    dubi->set_rhs( *NCap , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCap , NCap + ( rng.second - rng.first ) ,
	      U.begin() + rng.first );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
				              CapacitatedFacilityLocationBlockMod::eChgCaps , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_ucaps( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_ucaps( c_Vec_FNumber_it NCap , Subset && nms ,
			  const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( U.empty() ) {
  if( std::all_of( NCap , NCap + nms.size() ,
		   []( c_FNumber cap ) { return( cap >= Inf<FNumber>() ); }
		   ) )
   return;

  U.assign( get_MaxNArcs() , Inf<FNumber>() );
  }

 Index ndiff = countdiff( U , nms , NCap , get_NArcs() );
 if( ! ndiff )
  return;

 f_cond_lower = NAN;  // reset conditional bounds

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  if( HasDynamicX() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ;
	 ++NCap , ++nit ) {
     if( U[ *nit ] != *NCap ) {
      U[ *nit ] = *NCap;
      UB[ *nit ].set_rhs( *NCap , ampar );
      }
     }

    // dynamic part
    auto dubi = dUB.begin();
    for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dubi )
     if( *nit == i ) {
      if( U[ i ] != *NCap ) {
       U[ i ] = *NCap;
       dubi->set_rhs( *NCap , ampar );
       }
      nit++;
      NCap++;
      }
    }
   else {
    // make a vector of pairs < arc index , new capacity >
    typedef std::pair< Index , FNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i )
     pairs[ i ] = std::make_pair( nms[ i ] , *(NCap++) );

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end() ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( pit->first < get_NStaticArcs() ) ;
	 ++pit )
     if( U[ pit->first ] != pit->second ) {
      U[ pit->first ] = pit->second;
      UB[ pit->first ].set_rhs( *NCap , ampar );
      }

    // dynamic part
    auto dubi = dUB.begin();
    for( Index i = get_NStaticArcs() ; pit != pairs.end() ; ++i , ++dubi )
     if( pit->first == i ) {
      if( U[ i ] != pit->second ) {
       U[ i ] = pit->second;
       dubi->set_rhs( pit->second , ampar );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NCap , ++nit ) {
    if( U[ *nit ] != *NCap ) {
     U[ *nit ] = *NCap;
     UB[ *nit ].set_rhs( *NCap , ampar );
     }
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( U , nms , NCap );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  auto mod = std::make_shared<CapacitatedFacilityLocationBlockSbstMod>( this ,
                                  CapacitatedFacilityLocationBlockMod::eChgCaps , std::move( nms ) );

  Block::add_Modification( mod , Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_ucaps( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_ucap( c_FNumber NCap , c_Index arc ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( U.empty() && ( NCap < Inf<FNumber>() ) )
  U.assign( get_MaxNArcs() , Inf<FNumber>() );

 if( U[ arc ] == NCap )
  return;

 f_cond_lower = NAN;  // reset conditional bounds

 if( not_dry_run( issueMod ) )
  U[ arc ] = NCap;  // only change the physical representation - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  if( ! ( AR & HasBnd ) )
   throw( std::logic_error(
		"bound constraints not defined, cannot change capacity" ) );

  if( arc < get_NStaticArcs() )
   UB[ arc ].set_rhs( NCap , issueAMod );
  else
   std::next( dUB.begin() , arc - get_NStaticArcs() )->set_rhs( NCap ,
								issueAMod );
  }

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			   CapacitatedFacilityLocationBlockMod::eChgCaps , Range( arc , arc + 1 )  ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_ucap )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Range rng ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NNodes() );
 if( rng.second <= rng.first )  // nothing to change
  return;                 // cowardly (and silently) return

 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + ( rng.second - rng.first ) ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 c_Index ndiff = countdiff( NDfct , NDfct + ( rng.second - rng.first ) ,
			    B.cbegin() + rng.first );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticNodes() ) ;  ++i , ++NDfct )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    E[ i ].set_both( *NDfct , ampar );
    }

  // dynamic part
  for( auto dei = std::next( dE.begin() ,
			     rng.first >= get_NStaticNodes() ?
			     i - get_NStaticNodes() : 0 ) ;
       i < rng.second ; ++i , ++NDfct , ++dei )
   if( B[ i ] != *NDfct ) {
    B[ i ] = *NDfct;
    dei->set_both( *NDfct , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NDfct , NDfct + ( rng.second - rng.first ) ,
	      B.begin() + rng.first );

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
				              CapacitatedFacilityLocationBlockMod::eChgDfct , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_dfcts( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_dfcts( c_Vec_CNumber_it NDfct , Subset && nms ,
			  const bool ordered ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( B.empty() ) {
  if( std::all_of( NDfct , NDfct + nms.size() ,
		   []( c_FNumber dfct ) { return( dfct == 0 ); } ) )
   return;

  B.assign( get_MaxNNodes() , 0 );
  }

 Index ndiff = countdiff( B , nms , NDfct , get_NNodes() );
 if( ! ndiff )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  if( HasDynamicE() )
   if( ordered ) {
    // static part
    auto nit = nms.begin();
    for( ; ( nit != nms.end() ) && ( *nit < get_NStaticNodes() ) ;
	 ++NDfct , ++nit ) {
     if( B[ *nit ] != *NDfct ) {
      B[ *nit ] = *NDfct;
      E[ *nit ].set_both( *NDfct , ampar );
      }
     }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; nit != nms.end() ; ++i , ++dei )
     if( *nit == i ) {
      if( B[ i ] != *NDfct ) {
       B[ i ] = *NDfct;
       dei->set_both( *NDfct , ampar );
       }
      nit++;
      NDfct++;
      }
    }
   else {
    // make a vector of pairs < arc index , new capacity >
    typedef std::pair< Index , FNumber > index_pair;
    std::vector< index_pair > pairs( nms.size() );
    for( Index i = 0 ; i < nms.size() ; ++i )
     pairs[ i ] = std::make_pair( nms[ i ] , *(NDfct++) );

    // sort the vector for increasing index
    std::sort( pairs.begin() , pairs.end() ,
	       []( index_pair i , index_pair j )
	       { return( i.first < j.first ); } );

    // static part
    auto pit = pairs.begin();
    for( ; ( pit != pairs.end() ) && ( pit->first < get_NStaticArcs() ) ;
	 ++pit )
     if( B[ pit->first ] != pit->second ) {
      B[ pit->first ] = pit->second;
      E[ pit->first ].set_both( pit->second , ampar );
      }

    // dynamic part
    auto dei = dE.begin();
    for( Index i = get_NStaticNodes() ; pit != pairs.end() ; ++i , ++dei )
     if( pit->first == i ) {
      if( B[ i ] != pit->second ) {
       B[ i ] = pit->second;
       dei->set_both( pit->second , ampar );
       }
      pit++;
      }
    }
  else
   for( auto nit = nms.begin() ; nit != nms.end() ; ++NDfct , ++nit ) {
    if( B[ *nit ] != *NDfct ) {
     B[ *nit ] = *NDfct;
     E[ *nit ].set_both( *NDfct , ampar );
     }
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( B , nms , NDfct );

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockSbstMod>( this ,
                                 CapacitatedFacilityLocationBlockMod::eChgDfct , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_dfcts( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_dfct( c_CNumber NDfct , c_Index nde ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nde >= get_NNodes() )
  throw( std::invalid_argument( "invalid node name" ) );

 if( B.empty() && NDfct )
  B.assign( get_NNodes() , 0 );

 if( B[ nde ] == NDfct )
  return;

 if( not_dry_run( issueMod ) )
  B[ nde ] = NDfct;  // change the physical representation- - - - - - - - - -

 if( not_dry_run( issueAMod ) && ( AR & HasFlw ) ) {
  // change the abstract representation - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  if( nde < get_NStaticNodes() )
   E[ nde ].set_both( NDfct , issueAMod );
  else
   std::next( dE.begin() , nde - get_NStaticNodes() )->set_both( NDfct ,
								 issueAMod );
  }
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			   CapacitatedFacilityLocationBlockMod::eChgDfct , Range( nde , nde + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_dfct )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_arcs( Range rng ,
			   c_ModParam issueMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                 // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;
  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() )
    ndiff++;

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < rng.second ; )
   if( ! (dxi++)->is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( i = rng.first ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( ! x[ i ].is_fixed() ) {
    x[ i ].set_value( 0 );
    x[ i ].is_fixed( true , ampar );
    }

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < rng.second ; ++dxi )
   if( ! dxi->is_fixed() ) {
    dxi->set_value( 0 );
    dxi->is_fixed( true , ampar );
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
				             CapacitatedFacilityLocationBlockMod::eCloseArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_arcs( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_arcs( Subset && nms , const bool ordered  ,
			   c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() )
    ndiff++;

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( ! dxi->is_fixed() )
     ndiff++;
    ++nit;
    }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( nit = nms.begin() ; ( nit != nms.end() ) &&
	                   ( *nit < get_NStaticArcs() ) ; ++nit )
   if( ! x[ *nit ].is_fixed() ) {
    x[ *nit ].set_value( 0 );
    x[ *nit ].is_fixed( true , ampar );
    }

  // dynamic part
  dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( ! dxi->is_fixed() ) {
     dxi->set_value( 0 );
     dxi->is_fixed( true , ampar );
     }
    ++nit;
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockSbstMod>( this ,
                                 CapacitatedFacilityLocationBlockMod::eCloseArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_arc( c_Index arc ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( xa->is_fixed() )
   return;

  xa->set_value( 0 );

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( true , issueAMod );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			   CapacitatedFacilityLocationBlockMod::eCloseArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_arc )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_arcs( Range rng ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_NArcs() );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;
  Index i = rng.first;

  // static part
  for( ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( x[ i ].is_fixed() )
    ndiff++;

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < rng.second ; )
   if( (dxi++)->is_fixed() )
    ndiff++;

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( i = rng.first ; i < std::min( rng.second , get_NStaticArcs() ) ; ++i )
   if( x[ i ].is_fixed() )
    x[ i ].is_fixed( false , ampar );

  // dynamic part
  for( auto dxi = dx.begin() ; i++ < rng.second ; ++dxi )
   if( dxi->is_fixed() )
    dxi->is_fixed( false , ampar );

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 // TODO: if some changes are "fake", restrict the range

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
				             CapacitatedFacilityLocationBlockMod::eOpenArc , rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_arcs( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_arcs( Subset && nms , const bool ordered  ,
			  c_ModParam issueMod , c_ModParam issueAMod )
{
 if( nms.empty() )
  return;

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  Index ndiff = 0;

  // static part
  auto nit = nms.begin();
  for( ; ( nit != nms.end() ) && ( *nit < get_NStaticArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() )
    ndiff++;

  // dynamic part
  auto dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( dxi->is_fixed() )
     ndiff++;
    ++nit;
    }

  if( ! ndiff )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  c_ModParam ampar = make_amod_param( issueAMod , ndiff );

  // static part
  for( nit = nms.begin() ; ( nit != nms.end() ) &&
	                   ( *nit < get_NStaticArcs() ) ; ++nit )
   if( x[ *nit ].is_fixed() )
    x[ *nit ].is_fixed( false , ampar );

  // dynamic part
  dxi = dx.begin();
  for( Index i = get_NStaticArcs() ; nit != nms.end() ; ++i , ++dxi )
   if( *nit == i ) {
    if( dxi->is_fixed() )
     dxi->is_fixed( false , ampar );
    ++nit;
    }

  unmake_amod_param( issueAMod , ampar , ndiff );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 // TODO: eliminate from nms the "fake" changes

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockSbstMod>( this ,
                                  CapacitatedFacilityLocationBlockMod::eOpenArc , std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_arcs( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_arc( c_Index arc ,
			 c_ModParam issueMod , c_ModParam issueAMod )
{
 if( arc >= get_NArcs() )
  throw( std::invalid_argument( "invalid arc name" ) );

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway

 if( not_dry_run( issueAMod ) ) {
  auto xa = i2p_x( arc );

  if( ! xa->is_fixed() )
   return;

  // the physical and abstract representation are the same- - - - - - - - - -
  // change both (doh!), and if so instructed also issue abstract Modification

  xa->is_fixed( false , issueAMod );
  }

 f_cond_lower = NAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			   CapacitatedFacilityLocationBlockMod::eOpenArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_arc )

/*--------------------------------------------------------------------------*/

CapacitatedFacilityLocationBlock::Index CapacitatedFacilityLocationBlock::add_arc( c_Index sn , c_Index en ,
				   c_CNumber cst , c_FNumber cap ,
				   c_ModParam issueMod ,
				   c_ModParam issueAMod )
{
 if( ( sn < 1 ) || ( sn > get_NNodes() ) )
  throw( std::invalid_argument( "invalid starting node name" ) );
 
 if( ( en < 1 ) || ( en > get_NNodes() ) )
  throw( std::invalid_argument( "invalid ending node name" ) );

 Index arc = get_NStaticArcs();
 while( ( arc < get_NArcs() ) && ( ! is_deleted( arc ) ) )
  ++arc;

 if( arc >= get_MaxNArcs() )
  return( Inf<Index>() );

 // change the physical representation- - - - - - - - - - - - - - - - - - - -
 if( not_dry_run( issueMod ) ) {
  // set new arc cost
  if( C.empty() && cst )
   C.assign( get_MaxNArcs() , 0 );

  if( ! C.empty() )
   C[ arc ] = cst;

  // set new arc capacity
  if( U.empty() && ( cap < Inf<FNumber>() ) )
   U.assign( get_MaxNArcs() , Inf<FNumber>() );

  if( ! U.empty() )
   U[ arc ] = cap;

  // set contribution to flow constraint
  SN[ arc ] = sn;
  EN[ arc ] = en;
  }

 // change the abstract representation- - - - - - - - - - - - - - - - - - - -
 // in the meantime, if so instructed also issue abstract Modification(s)
 // note that this is *always* done, unless issueAMod says this is a dry
 // run, because at least the BlockModAdd corresponding to adding the
 // Variable, or unfixing it, is always issued since the Variable are
 // always present

 if( not_dry_run( issueAMod ) ) {

  c_ModParam ampar = make_amod_param( issueAMod ,
				      AR & ( HasFlw | HasObj ) ? 4 : 1 );
  ColVariable * nx;
  LB0Constraint * nUB;
  if( arc == get_NArcs() ) {
   // the new arc is physically constructed

   // create the new variable
   std::list< ColVariable > na;
   na.emplace_back( this , ColVariable::kNonNegative );
   nx = &(na.back());

   // now add it
   Block::add_dynamic_variables( dx , na , ampar );

   // add the new coefficient in the objective
   if( AR & HasObj )
    get_lfo()->add_variable( nx , cst , ampar );

   if( ( cap < Inf<FNumber>() ) && ( AR & HasFlw ) && ( ! ( AR & HasBnd ) ) )
    throw( std::logic_error( "cannot set finite capacity" ) );

   if( AR & HasBnd ) {
    // construct new arc capacity constraint
    std::list< LB0Constraint > nub;
    nub.emplace_back( this , nx );
    nUB = &(nub.back());

    // now add it
    Block::add_dynamic_constraints( dUB , nub , ampar );
    }
   }
  else {
   // the arc is just inserted in a previously deleted slot

   // recover pointer to the flow Variable
   nx = const_cast< ColVariable * >(
		   &( *std::next( dx.begin() , arc - get_NStaticArcs() ) ) );

   // un-fix the Variable
   nx->is_fixed( false , ampar );

   // recover pointer to the bound Constraint (if any)
   if( AR & HasBnd )
    nUB = const_cast< LB0Constraint * >(
		  &( *std::next( dUB.begin() , arc - get_NStaticArcs() ) ) );

   // change the cost coefficient in the objective
   if( AR & HasObj )
    get_lfo()->modify_coefficient( arc , cst , ampar );
   }

  // set new arc capacity: abstract part
  if( AR & HasBnd )
   nUB->set_rhs( cap , ampar );

  // set contribution to flow constraint: abstract part
  if( AR & HasFlw ) {
   get_lfc( i2p_e( sn - 1 ) )->add_variable( nx , -1 , ampar );
   get_lfc( i2p_e( en - 1 ) )->add_variable( nx ,  1 , ampar );
   }

  unmake_amod_param( issueAMod , ampar , AR & ( HasFlw | HasObj ) ? 4 : 1 );
  }

 if( arc == get_NArcs() )
  ++NArcs;  // increase arc count

 f_cond_lower = NAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
			    CapacitatedFacilityLocationBlockMod::eAddArc , Range( arc , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 return( arc );

 }  // end( CapacitatedFacilityLocationBlock::add_arc )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::remove_arc( c_Index arc , c_ModParam issueMod ,
		                         c_ModParam issueAMod )
{
 if( ( arc < get_NStaticArcs() ) || ( arc >= get_NArcs() ) )
  throw( std::invalid_argument( "invalid arc name" ) );

 if( is_deleted( arc ) )  // arc deleted already
  return;                 // nothing to do

 auto sn = SN[ arc ]; sn--;
 auto en = EN[ arc ]; en--;

 Index rmvdarcs = 1;  // how many arcs are removed in the end

 // change the physical representation- - - - - - - - - - - - - - - - - - - -
 if( not_dry_run( issueMod ) )
  SN[ arc ] = EN[ arc ] = Inf<Index>();

 // change the abstract representation- - - - - - - - - - - - - - - - - - - -
 // in the meantime, if so instructed also issue abstract Modification(s)
 // note that this is *always* done, unless issueAMod says this is a dry
 // run, because at least the BlockModAD corresponding to deleting the
 // Variable(s), or fixing them, is always issued since the Variable are
 // always present

 if( not_dry_run( issueAMod ) ) {
  c_ModParam ampar = make_amod_param( issueAMod ,
				      AR & ( HasFlw | HasObj ) ? 4 : 1 );
  if( arc == get_NArcs() - 1 ) {
   // removing the last arc (and possibly more)

   // reverse iterator into dx
   auto ritdx = dx.rbegin();

   // pointer to LinearFunction in the objective (if any)
   LinearFunction * lfo;
   if( AR & HasObj )
    lfo = get_lfo();

   // scan from the end backwards, eliminate all deleted arcs
   for( Index ai = arc ; ritdx != dx.rend() ; ++rmvdarcs , --ai ) {
    auto rxi = &(*(ritdx++));

    // delete contribution to objective (if any)
    if( AR & HasObj )
     lfo->remove_variable( ai , ampar );

    // delete contribution to flow constraint (if any)
    if( AR & HasFlw ) {
     auto snc = get_lfc( i2p_e( sn ) );
     auto sni = snc->is_active( rxi );
     if( sni >= snc->get_num_active_var() )
      throw( std::logic_error( "x variable not active in flow constraint" ) );
     snc->remove_variable( sni , ampar );
     auto enc = get_lfc( i2p_e( en ) );
     auto eni = enc->is_active( rxi );
     if( eni >= enc->get_num_active_var() )
      throw( std::logic_error( "x variable not active in flow constraint" ) );
     enc->remove_variable( eni , ampar );
     }

    if( rmvdarcs >= get_NArcs() - get_NStaticArcs() )
     break;

    if( ! is_deleted( get_NArcs() - rmvdarcs - 1 ) )
     break;
    }

   // define the range of removed stuff
   Range range( get_NArcs() - get_NStaticArcs() - rmvdarcs ,
		get_NArcs() - get_NStaticArcs() );
   
   // now actually remove and clear the UB Constraint(s) (if any)
   // do this before removing the flow Variable(s), so that if they are
   // processed in FIFO order it is seen before
   if( AR & HasBnd )
    Block::remove_dynamic_constraints( dUB , range , ampar );

   // now actually remove the flow Variable(s) (if any)
   Block::remove_dynamic_variables( dx , range , ampar );
   }
  else {
   // deleting one arc in the middle

   auto rx = const_cast< ColVariable * >(
		  &( *std::next( dx.begin() , arc - get_NStaticArcs() ) ) );

   rx->set_value( 0 );            // set the Variable to 0
   rx->is_fixed( true , ampar );  // fix it

   // delete contribution to flow constraint (if any)
   if( AR & HasFlw ) {
    auto snc = get_lfc( i2p_e( sn ) );
    auto sni = snc->is_active( rx );
    if( sni >= snc->get_num_active_var() )
     throw( std::logic_error( "x variable not active in flow constraint" ) );
    snc->remove_variable( sni , ampar );
    auto enc = get_lfc( i2p_e( en ) );
    auto eni = enc->is_active( rx );
    if( eni >= enc->get_num_active_var() )
     throw( std::logic_error( "x variable not active in flow constraint" ) );
    enc->remove_variable( eni , ampar );
    }
   }

  unmake_amod_param( issueAMod , ampar , AR & ( HasFlw | HasObj ) ? 4 : 1 );
  }
 else  // at the very least ensure the value is 0
  std::next( dx.begin() , arc - get_NStaticArcs() )->set_value( 0 );

 if( arc == get_NArcs() - 1 )
  NArcs -= rmvdarcs;  // decrease arc count

 f_cond_lower = NAN;  // reset conditional bounds

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<CapacitatedFacilityLocationBlockRngdMod>( this ,
				     CapacitatedFacilityLocationBlockMod::eRmvArc ,
				     Range( arc - rmvdarcs + 1 , arc + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::remove_arc )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::print( std::ostream &output ) const
{
 if( verbosity_lvl != Block::complete ) {  // non-complete version
  // only basic information 
  output << "CapacitatedFacilityLocationBlock with: " << NNodes << " nodes and " << SN.size()
	 << " arcs" << std::endl;

  if( verbosity_lvl == Block::high ) {     // print the graph
   for( Index i = 0 ; i < get_NNodes() ; ++i )
    if( B[ i ] != 0 )
     output << "B[ " << i + 1 << " ] = " << B[ i ] << std::endl;

   if( C.empty() )
    if( U.empty() )
     output << "all arcs have 0 cost and +Inf upper bound" << std::endl;
    else {
     for( Index i = 0 ; i < get_NArcs() ; ++i )
      if( ! is_deleted( i ) ) {
       output << "( " << SN[ i ] << " , " << EN[ i ] << " ): U = ";
       print_UB( output , U[ i ] );
       output << std::endl;
       }
     }
   else
    if( U.empty() )
     for( Index i = 0 ; i < get_NArcs() ; ++i ) {
      if( ! is_deleted( i ) )
       output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = " << C[ i ]
	      << std::endl;
      }
    else
     for( Index i = 0 ; i < get_NArcs() ; ++i )
      if( ! is_deleted( i ) ) {
       output << "( " << SN[ i ] << " , " << EN[ i ] << " ): C = " << C[ i ]
	      << ", U = ";
       print_UB( output , U[ i ] );
       output << std::endl;
       }
   }
  }
 else  {
  // print header in DIMACS standard format
  output << std::endl << "p min " << get_NNodes() << " ";
  if( HasDynamicX() ) {
   Index narcs = get_NStaticArcs();
   for( Index i = narcs ; i < get_NArcs() ; )
    if( ! is_deleted( i++ ) )
     ++narcs;
   
   output << narcs << std::endl;
   }
  else
   output << SN.size() << std::endl;

  // print node descriptors in DIMACS standard format
  for( Index i = 0 ; i < get_NNodes() ; ++i )
   if( B[ i ] != 0 )
    output << "n\t" << i + 1 << "\t" << - B[ i ] << std::endl;

  // print arc descriptors in DIMACS standard format
  if( C.empty() )
   if( U.empty() )
    for( Index i = 0 ; i < get_NArcs() ; ++i ) {
     if( ! is_deleted( i ) )
      output << "a\t" << SN[ i ] + 1 << "\t" << EN[ i ] + 1 << "\t0\t+Inf\t0"
	     << std::endl;
     }
   else {
    for( Index i = 0 ; i < get_NArcs()  ; ++i )
     if( ! is_deleted( i ) ) {
      output << "a\t" << SN[ i ] + 1 << "\t" << EN[ i ] + 1 << "\t0\t";
      print_UB( output , U[ i ] );
      output << "\t0" << std::endl;
      }
    }
  else
   if( U.empty() ) {
    for( Index i = 0 ; i < get_NArcs() ; ++i )
     if( ! is_deleted( i ) )
      output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t+Inf\t"
	     << C[ i ] << std::endl;
    }
   else
    for( Index i = 0 ; i < get_NArcs() ; ++i )
     if( ! is_deleted( i ) ) {
      output << "a\t" << SN[ i ] << "\t" << EN[ i ] << "\t0\t";
      print_UB( output , U[ i ] );
      output << "\t" << C[ i ] << std::endl;
      }
  }
 }  // end( CapacitatedFacilityLocationBlock::print )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_destructor( void )
{
 /* clear() all Constraint to ensure that they do not bother to un-register
  * themselves from Variable that are going to be deleted anyway. Then
  * deletes all the "abstract representation", if any.
  *
  * Note that this method is also called to reset the existing "abstract
  * representation" in case a nwe instance is loaded in the object; yet, even
  * in this case mo Modification pertaining to Variable and Constraint being
  * removed is necessary, because a NBModification is issued immediately
  * afterwards which means that any listening Observer already knows that
  * none of the previus Variable and Constraint are valid any longer. */

 for( auto & cnst : v_sat )  // clear the satisfaction constraints
  cnst.clear();
 for( auto & cnst : v_cap )  // clear the capacity constraints
  cnst.clear();
 for( auto & cnst : v_sfc )  // clear the strong forcin constraints
  cnst.clear();

 // then delete them all
 v_sfc.clear();
 v_cap.clear();
 v_sat.clear();

 // clear the objective function
 c.clear();

 // delete all Variable
 x.clear();
 y.clear();

 // delete all sub-Block
 for( auto bi : v_Block )
  delete bi;

 v_Block.clear();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer current) abstract representation 
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 AR = 0;  // no longer any abstract representation

 }  // end( CapacitatedFacilityLocationBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_Modification(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock is a "leaf" Block (has no sub-Block), this
  * method does not have to deal with GroupModification since these are
  * produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN CapacitatedFacilityLocationBlock WHEN THIS METHOD IS
  *   EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS ISSUED:
  *   NO COMPLCATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some of the logic here.*/

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast<LinearFunction * const>( tmod->function() );
  if( static_cast<LinearFunction * const>( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Range in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble

  if( tmod->range().second == tmod->range().first + 1 )
   // changing one cost only
   chg_cost( lfo->get_coefficient( tmod->range().first ) ,
	     tmod->range().first , make_par( eNoBlck , chnl ) , eDryRun );
  else {                            // changing many costs at once
   Vec_CNumber NC( tmod->range().second - tmod->range().first );
   auto NCit = NC.begin();
   for( Index i = tmod->range().first ; i < tmod->range().second ; )
    *(NCit++) = lfo->get_coefficient( i++ );

   chg_costs( NC.begin() , tmod->range() ,
	      make_par( eNoBlck , chnl ) , eDryRun );
   }

  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) ) {
  if( ! ( AR & HasObj ) )
   throw( std::invalid_argument( "Modification to non-constructed Objective"
				 ) );

  auto lfo = static_cast<LinearFunction * const>( tmod->function() );
  if( static_cast< LinearFunction * const >( c.get_function() ) != lfo )
   throw( std::invalid_argument( "Modification to non-Objective" ) );

  // note: in the following we can assume that the Subset in tmod is
  //       precisely the one we have to use since no Variable can have
  //       been added or deleted, which saves *a lot* of trouble
  // note: chg_costs() owns subset, so a copy has to be made

  Vec_CNumber NC( tmod->subset().size() );
  auto NCit = NC.begin();
  for( auto i : tmod->subset() )
   *(NCit++) = lfo->get_coefficient( i++ );

  chg_costs( NC.begin() , Subset( tmod->subset() ) , true ,
	     make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }


 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< RowConstraintMod * >( mod ) ) {
  if( ! ( AR & HasFlw ) )
   throw( std::invalid_argument(
			    "Modification to non-constructed Constraint" ) );

  if( tmod->type() == RowConstraintMod::eChgRHS ) {
   auto cp = dynamic_cast< LB0Constraint * const >( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_ucap( cp->get_rhs() , p2i_ub( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  if( tmod->type() == RowConstraintMod::eChgBTS ) {
   auto cp = static_cast<FRowConstraint * const>( tmod->constraint() );
   if( ! cp )
    throw( std::invalid_argument( "invalid Modification to Constraint" ) );

   chg_dfct( cp->get_rhs() , p2i_e( cp ) ,
	     make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  throw( std::invalid_argument( "illegal Modification to Constraint" ) );
  }

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  auto xi = dynamic_cast<ColVariable * const>( tmod->variable() );
  if( ! xi )
   throw( std::logic_error( "Modification to wrong type of Variable" ) );
  if( ( xi->get_type() != ColVariable::kNonNegative ) &&
      ( xi->get_type() != ColVariable::kNatural ) )
   throw( std::logic_error( "changing type of flow Variable not allowed" ) );
   
  auto i = p2i_x( xi );
  if( xi->is_fixed() )
   close_arc( i , make_par( eNoBlck , chnl ) , eDryRun );
  else
   open_arc( i , make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_copy(
		       CapacitatedFacilityLocationBlock * R3B , c_p_Mod mod ,
		       ModParam issuePMod , ModParam issueAMod )
{
 bool ok = true;  // final return value

 if( const auto tmod = dynamic_cast< GroupModification * const >( mod ) ) {
  // if the channels are the default ones, open new ones
  auto iPM = par2chnl( issuePMod ) ? issuePMod
           : make_par( par2mod( issuePMod ) , R3B->open_channel() );
  auto iPA = par2chnl( issueAMod ) ? issueAMod
           : make_par( par2mod( issueAMod ) , R3B->open_channel() );

  for( const auto & submod : tmod->sub_Modifications() )  // for each sub-Mod
   if( ! guts_of_guts_of_map_f_Mod_copy( R3B , submod.get() , iPM , iPA ) )
    ok = false;

  // now close the opened channels, if any
  if( ! par2chnl( issuePMod ) )
   R3B->close_channel( par2chnl( iPM ) );
  if( ! par2chnl( issueAMod ) )
   R3B->close_channel( par2chnl( iPA ) );
  }
 else  // any other Modification: just make the call
  ok = guts_of_guts_of_map_f_Mod_copy( R3B , mod , issuePMod , issueAMod );

 return( ok );
 
 }  // end( CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_copy )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_guts_of_map_f_Mod_copy(
		       CapacitatedFacilityLocationBlock * R3B , c_p_Mod mod ,
		       ModParam issuePMod , ModParam issueAMod )
{
 /* When a GroupModification is processed, if no channel is provided, then
    one is opened. This only happens "at root", after which in guts_of_mfM()
    whenever a GroupModification is processed, then the channel is nested.
    Indeed, if the "root" Modification is not a GroupModification, then there
    cannot be any GroupModification in it.

 ModParam iPM = issuePMod;
 ModParam iPA = make_par( std::min( ModParam( eNoBlck ) ,
				    par2mod( issueAMod ) ) ,
			  par2chnl( issueAMod ) );

 /* Use a Lambda to define a "guts" of the method that can be called
    recursively without having to pass "local globals". Note the trick of
    defining the std::function object and "passing" it to the lambda,
    which allows recursive calls. Note the need to explicitly capture
    "this" to use fields/methods of the class.

 std::function< bool( c_p_Mod )> guts_of_mfM;
 guts_of_mfM = [ this , & guts_of_mfM , & MCFB , & iPM , & iPA ]( c_p_Mod mod
								  ) {
  // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* This requires to patiently sift through the possible Modification types
     to find what this Modification exactly is, and call the appropriate
     method of either MCFB, for a "physical Modification", or of the "abstract
     representation" of MCFB for an "abstract Modification".

  //!! std::cout << *mod << std::endl;
  
  // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( const auto tmod = dynamic_cast< GroupModification * const >( mod ) ) {
   MCFB->nest_channel( par2chnl( iPM ) );  // nest the channel for PM
   MCFB->nest_channel( par2chnl( iPA ) );  // nest the channel for PA

   bool ok = true;
   for( const auto & submod : tmod->sub_Modifications() )
    if( ! guts_of_mfM( submod.get() ) )
     ok = false;

   MCFB->un_nest_channel( par2chnl( iPM ) );  // un-nest the channel for PM
   MCFB->un_nest_channel( par2chnl( iPA ) );  // un-nest the channel for PA

   return( ok );
   }

  // CapacitatedFacilityLocationBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note: in the following we can assume that C, B and U are nonempty. This
     is because they can be empty only if they are so when the object is
     loaded. But if a Modification has been issued they are no longer empty
     (a Modification changin nothing from the "empty" state is not issued).

  if( const auto tmod = dynamic_cast< CapacitatedFacilityLocationBlockRngdMod * const >( mod ) ) {
   switch( tmod->type() ) {
    case( CapacitatedFacilityLocationBlockMod::eChgCost ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_cost( C[ tmod->rng().first ] , tmod->rng().first ,
		      iPM , iPA );
     else
      MCFB->chg_costs( C.begin() + tmod->rng().first , tmod->rng() ,
		       iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eChgCaps ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_ucap( U[ tmod->rng().first ] , tmod->rng().first ,
		      iPM , iPA );
     else
      MCFB->chg_ucaps( U.begin() + tmod->rng().first , tmod->rng() ,
		       iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eChgDfct ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NNodes() ) ||
	  ( tmod->rng().second > MCFB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->chg_dfct( B[ tmod->rng().first ] , tmod->rng().first ,
		      iPM , iPA );
     else
      MCFB->chg_dfcts( B.begin() + tmod->rng().first , tmod->rng() ,
		       iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->open_arc( tmod->rng().first , iPM , iPA );
     else
      MCFB->open_arcs( tmod->rng() , iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->rng().second > get_NArcs() ) ||
	  ( tmod->rng().second > MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( tmod->rng().second == tmod->rng().first + 1 )
      MCFB->close_arc( tmod->rng().first , iPM , iPA );
     else
      MCFB->close_arcs( tmod->rng() , iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eAddArc ):
     #ifndef NDEBUG
      if( tmod->rng().first > get_NArcs() )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     if( MCFB->add_arc( get_SN( tmod->rng().first ) ,
			get_EN( tmod->rng().first ) ,
			get_C( tmod->rng().first ) ,
			get_U( tmod->rng().first ) , iPM , iPA )
	 != tmod->rng().first )
      throw( std::logic_error( "inconsistency between arc names" ) );       
     break;
    case( CapacitatedFacilityLocationBlockMod::eRmvArc ):
     #ifndef NDEBUG
      if( tmod->rng().first > MCFB->get_NArcs() )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     MCFB->remove_arc( tmod->rng().second - 1 , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown CapacitatedFacilityLocationBlockRngdMod type" ) );
    }
   return( true );
   }

  // CapacitatedFacilityLocationBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Note that tmod->f_nms need be copied, since the chg_*() methods
   * *in principle* "consume" the names vector. This is actually not true
   * if MCFB will *not* issue a physical modification, which one may
   * actually know beforehand, but it has to be done anyway because the
   * CapacitatedFacilityLocationBlockSbstMod only provides read-only access
 to the vector.

  if( const auto tmod = dynamic_cast< CapacitatedFacilityLocationBlockSbstMod * const >( mod ) ) {
   switch( tmod->type() ) {
    case( CapacitatedFacilityLocationBlockMod::eChgCost ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     Vec_CNumber NCost( tmod->nms().size() );
     for( Index i = 0 ; i < NCost.size() ; i++ )
      NCost[ i ] = C[ tmod->nms()[ i ] ];

     MCFB->chg_costs( NCost.begin() , Subset( tmod->nms() ) , iPM , iPA );
     break;
     }
    case( CapacitatedFacilityLocationBlockMod::eChgCaps ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     Vec_FNumber NCap( tmod->nms().size() );
     for( Index i = 0 ; i < NCap.size() ; i++ )
      NCap[ i ] = U[ tmod->nms()[ i ] ];

     MCFB->chg_ucaps( NCap.begin() , Subset( tmod->nms() ) , iPM , iPA );
     break;
     }
    case( CapacitatedFacilityLocationBlockMod::eChgDfct ): {
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NNodes() ) ||
	  ( tmod->nms().back() >= MCFB->get_NNodes() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     Vec_FNumber NDfct( tmod->nms().size() );
     for( Index i = 0 ; i < NDfct.size() ; i++ )
      NDfct[ i ] = B[ tmod->nms()[ i ] ];

     MCFB->chg_dfcts( NDfct.begin() , Subset( tmod->nms() ) , iPM , iPA );
     break;
     }
    case( CapacitatedFacilityLocationBlockMod::eOpenArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     MCFB->open_arcs( Subset( tmod->nms() ) , iPM , iPA );
     break;
    case( CapacitatedFacilityLocationBlockMod::eCloseArc ):
     #ifndef NDEBUG
      if( ( tmod->nms().back() >= get_NArcs() ) ||
	  ( tmod->nms().back() >= MCFB->get_NArcs() ) )
       throw( std::logic_error(
		     "map_forward_Modification:: incompatible CapacitatedFacilityLocationBlock" ) );
     #endif
     MCFB->close_arcs( Subset( tmod->nms() ) , iPM , iPA );
     break;
    default:
     throw( std::invalid_argument( "unknown CapacitatedFacilityLocationBlockSbstMod type" ) );
    }
   return( true );
   }

  // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is the "nuclear option": the CapacitatedFacilityLocationBlock has been re-loaded
  // one should check that the Block is this CapacitatedFacilityLocationBlock, but it cannot
  // be otherwise, can it?

  if( const auto tmod = dynamic_cast< NBModification * const >( mod ) ) {
   MCFB->load( get_NNodes() , get_NArcs() , EN , SN , U , C , B ,
	       get_NNodes() - get_NStaticNodes() ,
	       get_NArcs() - get_NStaticArcs() ,
	       get_MaxNNodes() - get_NStaticNodes() ,
	       get_MaxNArcs() - get_NStaticArcs() );
   return( true );
   }

  return( false );

  };  // end( guts_of_mfM )- - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // finally, call the "guts of"- - - - - - - - - - - - - - - - - - - - - - - -
 // this is done differently if mod is a GroupModification, since at the root
 // a channel has to be opened while further down it has to be nested

 */

 return( false );  // any other Modification is not mapped

 }  // end( CapacitatedFacilityLocationBlock::guts_of_guts_of_map_f_Mod_copy )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_MCF(
				    MCFBlock * R3B , c_p_Mod mod ,
			            ModParam issuePMod , ModParam issueAMod )
{
 bool ok = true;  // final return value

 if( const auto tmod = dynamic_cast< GroupModification * const >( mod ) ) {
  // if the channels are the default ones, open new ones
  auto iPM = par2chnl( issuePMod ) ? issuePMod
           : make_par( par2mod( issuePMod ) , R3B->open_channel() );
  auto iPA = par2chnl( issueAMod ) ? issueAMod
           : make_par( par2mod( issueAMod ) , R3B->open_channel() );

  for( const auto & submod : tmod->sub_Modifications() )  // for each sub-Mod
   if( ! guts_of_guts_of_map_f_Mod_MCF( R3B , submod.get() , iPM , iPA ) )
    ok = false;

  // now close the opened channels, if any
  if( ! par2chnl( issuePMod ) )
   R3B->close_channel( par2chnl( iPM ) );
  if( ! par2chnl( issueAMod ) )
   R3B->close_channel( par2chnl( iPA ) );
  }
 else  // any other Modification: just make the call
  ok = guts_of_guts_of_map_f_Mod_MCF( R3B , mod , issuePMod , issueAMod );

 return( ok );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_MCF )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_guts_of_map_f_Mod_MCF(
				    MCFBlock * R3B , c_p_Mod mod ,
			            ModParam issuePMod , ModParam issueAMod )
{


 return( false );  // any other Modification is not mapped
 
 }  // end( CapacitatedFacilityLocationBlock::guts_of_guts_of_map_f_Mod_MCF )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::compute_conditional_bounds( void )
{
 f_cond_lower = f_cond_upper = 0;

 for( Index i = 0 ; i < f_n_facilities ; ++i )
  if( v_fixed_cost[ i ] >= 0 )
   f_cond_upper += v_fixed_cost[ i ];
  else
   f_cond_lower += v_fixed_cost[ i ];

 for( Index j = 0 ; j < f_n_customers ; ++j ) {
  auto minj = Inf< TCost >();
  auto maxj = - Inf< TCost >();

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   if( minj > v_transp_cost[ j ][ i ] )
    minj = v_transp_cost[ j ][ i ];
   if( maxj < v_transp_cost[ j ][ i ] )
    maxj = v_transp_cost[ j ][ i ];
   }

  f_cond_lower += minj;
  f_cond_upper += maxj;
  }
 }  // end( CapacitatedFacilityLocationBlock::compute_conditional_bounds )

/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::get_y(
				 std::vector< T >::iterator Sol , Range rng )
{
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 #ifndef NDEBUG
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( "get_facility_solution: variables not generated"
			    ) );
  if( rng.second > f_n_facilities )
   rng.second = f_n_facilities;
 #endif

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( Index i = rng.first ; i < rng.second ; )
   *(Sol++) = BKB( v_Block[ i ] )->get_x( f_n_customers );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index i = rng.first ; i < rng.second ; )
  *(Sol++) = v_y[ i++ ].get_value();

 }  // end( get_y( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::get_y(
			      std::vector< T >::iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "get_facility_solution: ";
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std__any_of( nms.begin() , nms.end ,
		   [ & f_n_facilities ]( auto i ) {
		    return( i > f_n_facilities );
		    }
		   ) )
   throw( std::logic_error( _prfx + "invalid facility index in nms" ) );
 #endif

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( auto i : nms )
   *(Sol++) = BKB( v_Block[ i ] )->get_x( f_n_customers );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( auto i : nms )
  *(Sol++) = v_y[ i ].get_value();

 }  // end( get_y( Subset ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::get_x(
				 std::vector< T >::iterator Sol , Range rng )
{
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 #ifndef NDEBUG
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( "get_transportation_solution: "
			    "variables not generated" ) );
  if( rng.second > f_n_facilities * f_n_customers )
   rng.second = f_n_facilities * f_n_customers;
 #endif

 if( ( AR & FormMsk ) == StdForm ) {  // standard formulation- - - - - - - - -
  for( Index h = rng.first ; h < rng.second ; )
   *(Sol++) = (v_x.data())[ h++ ].get_value();
  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  Index i = rng.first / f_n_customers;
  Index frst = rng.first % f_n_customers;
  // simple case: it's always the same facility
  if( i == rng.second / f_n_customers ) {
   BKB( v_Block[ i ] )->get_x( Sol , Range( frst ,
					    rng.second % f_n_customers ) );
   return;
   }
  // it's at least two facilities
  for( Index h = rng.first ; h < rng.second ; ++h , ++i ) {
   Index scnd = std::min( rng.second - h , f_n_customers );
   BKB( v_Block[ i ] )->get_x( Sol , Range( frst , second ) );
   if( scnd < f_n_customers )
    break;
   scnd -= frst;
   frst = 0;
   h += scnd;
   Sol += scnd;
   }
  return;
  }

 // flow formulation: first get it from the MCFBlock, but it is scaled- - - -
 MCFB( v_Block[ 1 ] )->get_x( Sol , Range( rng.first + f_n_facilities ,
					   rng.second + f_n_facilities ) );
 // now de-scale it
 for( Index h = rng.first ; h < rng.second ; ++h )
  *(Sol++) /= v_demand[ h % f_n_customers ];

 }  // end( get_x( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::get_x(
			      std::vector< T >::iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "get_transportation_solution: ";
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std__any_of( nms.begin() , nms.end ,
		   [ & f_n_facilities ]( auto i ) {
		    return( i > f_n_facilities );
		    }
		   ) )
   throw( std::logic_error( _prfx + "invalid index in nms" ) );
 #endif


 if( ( AR & FormMsk ) == StdForm ) {  // standard formulation- - - - - - - - -
  for( auto i : nms )
   *(Sol++) = (v_x.data())[ i ].get_value();
  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  auto nbit = nms.begin();
  for( Index i = *nbit / f_n_customers ; nbit != nms.end() ; ++i ) {
   auto neit = ++nbit;
   while( ( neit != nms.end() ) && ( *neit / f_n_customers == i ) )
    ++neit;
   Subset nnms( nbit , neit );
   for( auto & nm : nnms )
    nm %= f_n_customers;
   BKB( v_Block[ i ] )->get_x( Sol , nnms );
   Sol += nnms.size();
   nbit = neit;
   }
  return;
  }

 // flow formulation: first get it from the MCFBlock, but it is scaled- - - -
 Subset nnms( nms );
 for( auto & nm : nnms )
  nm -= f_n_facilities;
 MCFB( v_Block[ 1 ] )->get_x( Sol , nnms );
 // now de-scale it
 for( Index h : nms )
  *(Sol++) /= v_demand[ h % f_n_customers ];

 }  // end( get_x( Subset ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::set_y(
				 std::vector< T >::iterator Sol , Range rng )
{
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 #ifndef NDEBUG
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( "set_facility_solution: variables not generated"
			    ) );
  if( rng.second > f_n_facilities )
   rng.second = f_n_facilities;
 #endif

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( Index i = rng.first ; i < rng.second ; )
   BKB( v_Block[ i ] )->set_x( f_n_customers , *(Sol++) );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index i = rng.first ; i < rng.second ; )
  v_y[ i++ ].set_value( *(Sol++) );

 }  // end( set_y( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::set_y(
			      std::vector< T >::iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "set_facility_solution: ";
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std__any_of( nms.begin() , nms.end ,
		   [ & f_n_facilities ]( auto i ) {
		    return( i > f_n_facilities );
		    }
		   ) )
   throw( std::logic_error( _prfx + "invalid facility index in nms" ) );
 #endif

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( auto i : nms )
   BKB( v_Block[ i ] )->set_x( f_n_customers , *(Sol++) );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( auto i : nms )
  v_y[ i ].set_value( *(Sol++) );

 }  // end( set_y( Subset ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::set_x(
				 std::vector< T >::iterator Sol , Range rng )
{
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 #ifndef NDEBUG
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( "set_transportation_solution: "
			    "variables not generated" ) );
  if( rng.second > f_n_facilities * f_n_customers )
   rng.second = f_n_facilities * f_n_customers;
 #endif

 if( ( AR & FormMsk ) == StdForm ) {  // standard formulation- - - - - - - - -
  for( Index h = rng.first ; h < rng.second ; )
   (v_x.data())[ h++ ].set_value( *(Sol++) );
  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  Index i = rng.first / f_n_customers;
  Index frst = rng.first % f_n_customers;
  // simple case: it's always the same facility
  if( i == rng.second / f_n_customers ) {
   BKB( v_Block[ i ] )->set_x( Sol , Range( frst ,
					    rng.second % f_n_customers ) );
   return;
   }
  // it's at least two facilities
  for( Index h = rng.first ; h < rng.second ; ++h , ++i ) {
   Index scnd = std::min( rng.second - h , f_n_customers );
   BKB( v_Block[ i ] )->set_x( Sol , Range( frst , second ) );
   if( scnd < f_n_customers )
    break;
   scnd -= frst;
   frst = 0;
   h += scnd;
   Sol += scnd;
   }
  return;
  }

 // flow formulation: first need to scale the solution- - - - - - - - - - - -
 MCFBlock::Vec_FNumber SSol( Sol , Sol + ( rng.second - rng.first ) );
 auto Sit = SSol.begin();
 for( Index h = rng.first ; h < rng.second ; ++h )
  *(Sit++) *= v_demand[ h % f_n_customers ];

 // now set it in the MCFBlock
 MCFB( v_Block[ 1 ] )->set_x( SSol.begin() ,
			      Range( rng.first + f_n_facilities ,
				     rng.second + f_n_facilities ) );
 }  // end( set_x( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::set_x(
			      std::vector< T >::iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "set_transportation_solution: ";
  if( ! ( AR &= HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std__any_of( nms.begin() , nms.end ,
		   [ & f_n_facilities ]( auto i ) {
		    return( i > f_n_facilities );
		    }
		   ) )
   throw( std::logic_error( _prfx + "invalid index in nms" ) );
 #endif

 if( ( AR & FormMsk ) == StdForm ) {  // standard formulation- - - - - - - - -
  for( auto i : nms )
   (v_x.data())[ i ].set_value( *(Sol++) );
  return;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  auto nbit = nms.begin();
  for( Index i = *nbit / f_n_customers ; nbit != nms.end() ; ++i ) {
   auto neit = ++nbit;
   while( ( neit != nms.end() ) && ( *neit / f_n_customers == i ) )
    ++neit;
   Subset nnms( nbit , neit );
   for( auto & nm : nnms )
    nm %= f_n_customers;
   BKB( v_Block[ i ] )->set_x( Sol , nnms );
   Sol += nnms.size();
   nbit = neit;
   }
  return;
  }

 // flow formulation: first need to scale the solution- - - - - - - - - - - -
 MCFBlock::Vec_FNumber SSol( Sol , Sol + nms.size() );
 auto Sit = SSol.begin();
 for( Index h : nms )
  *(Sit++) *= v_demand[ h % f_n_customers ];

 // now set it in the MCFBlock
 Subset nnms( nms );
 for( auto & nm : nnms )
  nm -= f_n_facilities;
 MCFB( v_Block[ 1 ] )->set_x( SSol , nnms );

 }  // end( set_x( Subset ) )
 
/*--------------------------------------------------------------------------*/

ModParam CapacitatedFacilityLocationBlock::make_amod_param(
					     ModParam issueAMod , Index num )
{
 if( issue_mod( issueAMod ) ) {
  ChnlName chnl = par2chnl( issueAMod );
  if( num > 1 ) {            // more than one Modification have to be issued
    if( chnl )               // and a channel is already provided
     nest_channel( chnl );   // nest the channel
    else                     // it was being sent to default channel
     chnl = open_channel();  /* open a new channel: note that the
			      * GroupModification will automatically be a
     * "physical Modification" (i.e., concerns_Block() == false) since such
     * are all the Modification there inside: in fact, all the inner
     * Modification will be issued with the return value, which i eNoBlck */
   }

  return( make_par( eNoBlck , chnl ) );
  }
 else
  return( eNoMod );
 }

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::unmake_amod_param( ModParam oldiAM ,
						ModParam newiAM , Index num )
{
 if( newiAM == eNoMod )
  return;

 ChnlName chnl = par2chnl( newiAM );
 if( num > 1 ) {               // a channel had been opened/nested
  if( par2chnl( oldiAM ) )     // that's "nested"
   un_nest_channel( chnl );    // un-nest it
  else                         // that's "opened"
   close_channel( chnl );      // close it
  }
 }

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void CapacitatedFacilityLocationBlock::CheckAbsVSPhys( void )
{
 // check that the (part that has actually been constructed of the) abstract
 // representation coincides with the physical representation


 }  // end( CapacitatedFacilityLocationBlock::CheckAbsVSPhys )

#endif

/*--------------------------------------------------------------------------*/
/*------------- METHODS OF CapacitatedFacilityLocationSolution -------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::deserialize(
					      const netCDF::NcGroup & group )
{
 const std::string _prfx = "CapacitatedFacilityLocationSolution: ";
 netCDF::NcDim nf = group.getDim( "NFacilities" );
 if( nf.isNull() )
  throw( std::invalid_argument( _pfrx + "NFacilities dimension required" ) );

 netCDF::NcDim nc = group.getDim( "NCustomers" );
 if( nc.isNull() )
  throw( std::invalid_argument( _pfrx + "NCustomers dimension required" ) );

 netCDF::NcVar fs = group.getVar( "FacilityCapacitySolution" );
 if( fs.isNull() )
  v_y.clear();
 else {
  v_y.resize( nf.getSize() );
  fs.getVar( v_y.data() );
  }

 netCDF::NcVar ts = group.getVar( "TransportationSolution" );
 if( ts.isNull() )
  v_x.clear();
 else {
  v_pi.resize( { nc.getSize() , nf.getSize() } );
  ts.getVar( v_x.data() );
  }
 }  // end( CapacitatedFacilityLocationSolution::deserialize )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::read( const Block * block )
{
 auto CFLB = dynamic_cast<const CapacitatedFacilityLocationBlock * >( block );
 if( ! CFLB )
  throw( std::invalid_argument(
		        "block is not a CapacitatedFacilityLocationBlock" ) );

 // read y- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_y.empty() ) {
  if( v_x.size() < CFLB->get_NFacilities() )
   v_x.resize( CFLB->get_NFacilities );

  CFLB->get_facility_solution( v_y.begin() );
  }

 // read x- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( ( (v_x.shape())[ 0 ] != CFLB->get_NFacilities() ) ||
      ( (v_x.shape())[ 1 ] != CFLB->get_NCustomers() ) )
   v_x.resize( { CFLB->get_NFacilities() , CFLB->get_NCustomers() } );

  CFLB->get_transportation_solution( v_x.data().begin() );
  }
 }  // end( CapacitatedFacilityLocationSolution::read )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::write( Block * const block ) 
{
 auto CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( block );
 if( ! CFLB )
  throw( std::invalid_argument(
		       "block is not a CapacitatedFacilityLocationBlock" ) );

 // write y - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_y.empty() ) {
  if( v_y.size() < CFLB->get_NFacilities() )
   throw( std::invalid_argument( "incompatible facility size" ) );

  CFLB->set_facility_solution( v_y.begin() );
  }

 // write x - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( ( (v_x.shape())[ 0 ] != CFLB->get_NFacilities() ) ||
      ( (v_x.shape())[ 1 ] != CFLB->get_NCustomers() ) )
   throw( std::invalid_argument( "incompatible transportation size" ) );

  CFLB->set_transportation_solution( v_x.data().begin() );
  }
 }  // end( CapacitatedFacilityLocationSolution::write )

/*--------------------------------------------------------------------------*/

void MCFSolution::serialize( netCDF::NcGroup & group )
{
 std::vector<size_t> startp = { 0 };

 if( ! v_x.empty() ) {
  netCDF::NcDim na = group.addDim( "NumArcs" , v_x.size() );

  std::vector<size_t> countpa = { v_x.size() };

  ( group.addVar( "FlowSolution" , netCDF::NcDouble() , na ) ).putVar(
					      startp , countpa , v_x.data() );
  }

 if( v_pi.empty() )
  return;

 netCDF::NcDim nn = group.addDim( "NumNodes" ,  v_pi.size() );
 std::vector<size_t> countpn = { v_pi.size() };
 ( group.addVar( "Potentials" , netCDF::NcDouble() , nn ) ).putVar(
					     startp , countpn , v_pi.data() );
 
 }  // end( MCFSolution::serialize )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::scale( double factor ) const
{
 auto * sol = MCFSolution::clone( true );

 if( ! v_x.empty() )
  for( CapacitatedFacilityLocationBlock::Index i = 0 ; i < v_x.size() ; ++i )
   sol->v_x[ i ] = v_x[ i ] * factor;

 if( ! v_pi.empty() )
  for( CapacitatedFacilityLocationBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   sol->v_pi[ i ] = v_pi[ i ] * factor;

 return( sol );

 }  // end( MCFSolution::scale )

/*--------------------------------------------------------------------------*/

void MCFSolution::sum( const Solution * solution , double multiplier )
{
 auto MCFS = dynamic_cast< const MCFSolution * >( solution );
 if( ! MCFS )
  throw( std::invalid_argument( "solution is not a MCFSolution" ) );

 if( ! v_x.empty() ) {
  if( v_x.size() != MCFS->v_x.size() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  for( CapacitatedFacilityLocationBlock::Index i = 0 ; i < v_x.size() ; ++i )
   v_x[ i ] = MCFS->v_x[ i ] * multiplier;
  }

 if( ! v_pi.empty() ) {
  if( v_pi.size() != MCFS->v_pi.size()  )
   throw( std::invalid_argument( "incompatible potential size" ) );

  for( CapacitatedFacilityLocationBlock::Index i = 0 ; i < v_pi.size() ; ++i )
   v_pi[ i ] = MCFS->v_pi[ i ] * multiplier;
  }
 }  // end( MCFSolution::sum )

/*--------------------------------------------------------------------------*/

MCFSolution * MCFSolution::clone( bool empty ) const
{
 auto *sol = new MCFSolution();

 if( empty ) {
  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );

  if( ! v_pi.empty() )
   sol->v_pi.resize( v_pi.size() );
  }
 else {
  sol->v_x = v_x;
  sol->v_pi = v_pi;
  }

 return( sol );

 }  // end( MCFSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------- End File CapacitatedFacilityLocationBlock.cpp --------------*/
/*--------------------------------------------------------------------------*/
