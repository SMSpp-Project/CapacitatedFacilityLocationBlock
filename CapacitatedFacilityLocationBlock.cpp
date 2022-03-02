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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static LinearFunction * LF( Function * f ) {
 return( static_cast< LinearFunction * >( f ) );
 }

/*--------------------------------------------------------------------------*/
// returns true if two vectors differ, one of them being given as a base
// vector and a subset of indices

template< typename T >
static bool is_equal( std::vector<T> & vec , c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
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

/*--------------------------------------------------------------------------*/
// copies one vector to a given subset of another

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  vec[ nm ] = *(cpy++);
 }

/*----------------------------------------------------------------------------
// re-order the parallel vec and nms in increasing order of nms

template< typename T >
static void copyidx( std::vector< T > & vec , c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
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
 v_f_cost = std::move( F );
 v_demand = std::move( D );
 v_t_cost = std::move( C );

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
 v_f_cost.resize( f_n_facilities );

 for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
  input >> eatcomments >> v_capacity[ i ];
  if( input.fail() )
   goto( input_failure );
  if( v_capacity[ i ] <= 0 )
   throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

  input >> eatcomments >> v_f_cost[ i ];
  if( input.fail() )
   goto( input_failure );

  }  // end( for( each facility ) )

 v_demand.resize( f_n_customers );
 v_t_cost.resize( boost::extents[ f_n_facilities ][ f_n_customers ] );

 for( Index j = 0 ; j < f_n_customers ; ++j ) {  // for( each customer )
  input >> eatcomments >> v_demand[ j ];
  if( input.fail() )
   goto( input_failure );
  if( v_demand[ j ] <= 0 )
   throw( std::invalid_argument( _prfx + "non-positive demand" ) );

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
   input >> eatcomments >> v_t_cost[ i ][ j ];
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
 
 v_f_cost.resize( f_n_facilities );
 fq.getVar( v_f_cost.data() );

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

 v_t_cost.resize( tcs );
 tc.getVar( std::vector< std::size_t >( 2 , 0 ) , tcs ,
	    v_t_cost.data() );

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

  // first construct the vector and sort it, so that the pointers are
  // increasing with the facility index i, which speeds up some operations
  for( auto & bi : v_Block )
   bi = new BinaryKnapsackBlock( this );

  std::sort( v_Block.begin() , v_Block.end() );

  // now load the appropriate data into each BinaryKnapsackBlock
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
    C[ j ] = v_t_cost[ i ][ j ];
    }
   W[ f_n_customers ] = - v_capacity[ i ];
   C[ f_n_customers ] = v_f_cost[ i ];

   v_Block[ i ]->load( f_n_customers + 1 , 0 , W , P , I );
   v_Block[ i ]->set_objective_sense( false , eNoMod , eNoMod );
   v_Block[ i ]->generate_abstract_variables();
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
   *(pi++) = std::make_pair( & v_y[ i ] , v_f_cost[ i ] );

  // then the X[ j ][ i ] ones
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    for( Index j = 0 ; j < f_n_customers ; ++j )
     *(pi++) = std::make_pair( & v_x[ i ][ j ] , v_t_cost[ i ][ j ] );

  f_obj.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
  set_objective( & f_obj , eNoMod );
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
  p[ i ] = std::make_pair( & v_y[ i ] , v_f_cost[ i ] );

 f_obj.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 AB( v_Block[ 0 ] )->set_objective( & f_obj , eNoMod );

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

  CFLB->load( f_n_facilities , f_n_customers , v_capacity , v_f_cost ,
	      v_demand , v_t_cost );
 
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
 while( i < f_n_facilities )  // facilities nodes
  B[ i++ ] = 0;

 FNumber todD = 0;
 for( Index j = 0 ; j < f_n_customers ; j++ ) {  // customers nodes
  todD -= v_demand[ j ];
  B[ i++ ] = v_demand[ j ];
  }

 B[ i ] = todD;  // super-source

 // construct arcs SN, EN, U, C: "common part" of the graph
 Index a = 0;
 const Index ss = f_n_facilities + f_n_customers + 1;

 // first the source -> facility arcs
 for( i = 0 ; i < f_n_facilities ; ++i ) {
  SN[ a ] = 22;
  EN[ a ] = i + 1;
  U[ a ] = v_capacity[ i ];
  C[ a++ ] = v_f_cost[ i ] / v_capacity[ i ];
  }

 // now the facility -> customers arcs
  for( i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j ) {
    SN[ a ] = i + 1;
    EN[ a ] = f_n_facilities + 1 + j;
    U[ a ] = Inf< MCFClass::FNumber >();
    C[ a++ ] = v_t_cost[ i ][ j ] / v_demand[ j ];
    }

 if( wR3B > 1 ) {
  // now the artificial arcs to ensure feasibility

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   SN[ a ] = ss;
   EN[ a ] = f_n_facilities + 1 + j;
   U[ a ] = Inf< MCFClass::FNumber >();

   // compute an upper bound on the worst-case transportation cost
   MCFClass::CNumber maxc = 0;
   for( i = 0 ; i < f_n_facilities ; ++i , ++a )
    if( auto tci = C[ i ] + v_t_cost[ i ][ j ] / v_demand[ j ] ;
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
  sol->v_x.resize( f_n_customers * f_n_facilities );

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
  switch( AR & FormMsk ) {
   case( StdForm ): guts_of_add_ModificationSF( mod.get() , chnl ); break;
   case( KskForm ): guts_of_add_ModificationKF( mod.get() , chnl ); break;
   default:         guts_of_add_ModificationFF( mod.get() , chnl );
   }
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
   ).putVar( v_f_cost.data() );

 ( group.addVar( "CustomerDemand" , netCDF::NcUint64() , nc )
   ).putVar( v_demand.data() );

 ::serialize( group , "TransportationCost" , netCDF::NcDouble() ,
              { nc , nf } , v_t_cost );

 }  // end( CapacitatedFacilityLocationBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_costs(
				     c_CV_it NCost , Range rng ,
				     ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_facilities );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 const Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NCost , NCost + num , v_f_cost.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCost , NCost + num , v_f_cost.begin() + rng.first );

  if( ( AR & FormMsk ) != KskForm )
   // since modify_coefficients owns the vector, a copy has to be made
   get_lfo()->modify_coefficients( CVector( NCost , NCost + num ) , rng ,
				   issueAMod );
  else
   for( Index i = rng.first ; i < rng.second ; ++i )
    KB( v_Block[ i ] )->chg_weight( *(NCost++) , f_n_customers , issueMod ,
				    issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCost , NCost + num , v_f_cost.begin() + rng.first );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgDCost ,
			    rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_facility_costs( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_costs(
		               c_CV_it NCost , Subset && nms , bool ordered ,
			       ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_f_cost , nms , NCost , f_n_facilities ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_f_cost , nms , NCost );

  if( ( AR & FormMsk ) != KskForm )
   // since modify_coefficients owns both vectors, two copies are made
   get_lfo()->modify_coefficients( CVector( NCost , NCost + nms.size() ) ,
				   Subset( nms ) , true , issueAMod );
  else
   for( Index i = rng.first ; i < rng.second ; ++i )
    KB( v_Block[ i ] )->chg_weight( *(NCost++) , f_n_customers , issueMod ,
				    issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_f_cost , nms , NCost );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
                            CapacitatedFacilityLocationBlockMod::eChgFCost ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_facility_costs( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_cost( Cost NCost ,
			   Index i , ModParam issueMod , ModParam issueAMod )
{
 if( i >= f_n_facilities )
  throw( std::invalid_argument( "invalid facility name" ) );

 if( v_f_cost[ i ] == NCost )
  return;

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  v_f_cost[ i ] = NCost;

  if( ( AR & FormMsk ) != KskForm )
   get_lfo()->modify_coefficient( i , NCost , issueAMod );
  else
   KB( v_Block[ i ] )->chg_weight( NCost , f_n_customers , issueMod ,
				   issueAMod );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   v_f_cost[ i ] = NCost;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgFCost ,
			    Range( i , i + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_facility_cost )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_costs(
				     c_CV_it NCost , Range rng ,
				     ModParam issueMod , ModParam issueAMod )
{
 const Index maxn = f_n_facilities * f_n_customers;
 rng.second = std::min( rng.second , maxn );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 const Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NCost , NCost + num , v_t_cost.data() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCost , NCost + num , v_t_cost.data() + rng.first );

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // since modify_coefficients owns the vector, a copy has to be made
    CVector NC( NCost , NCost + num );
    get_lfo()->modify_coefficients( std::move( NC ) ,
				    Range( rng.first + f_n_facilities ,
					   rng.second + f_n_facilities ) ,
				    issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    Index f = rng.first;
    Index i = f / f_n_customers;
    Index l = f % f_n_customers;
    if( ( ( rng.second - 1 ) / f_n_customers ) == i ) {
     // the range is all inside a single facility
     KB( v_Block[ i ] )->chg_weighs( NCost , Range( l , l + num ) ,
				     issueMod , issueAMod );
     break;
     }

    // the range of the first facility does not necessarily start from 0
    // but it surely ends at f_n_customers
    KB( v_Block[ i++ ] )->chg_weighs( NCost , Range( l , f_n_customers ) ,
				      issueMod , issueAMod );
    NCost += ( f_n_customers - l );
    f += ( f_n_customers - l );
 
    // the range of all other facilities starts from 0 but it does not
    // necessarily end at f_n_customers
    for( ; ; ++i , NCost += f_n_customers ) {
     Index nf = f + f_n_customers;
     if( nf >= rng.second ) {  // last facility
      KB( v_Block[ i ] )->chg_weighs( NCost , Range( 0 , rng.second - f ) ,
				      issueMod , issueAMod );
      break;
      }
     else {
      KB( v_Block[ i ] )->chg_profits( NCost , Range( 0 , f_n_customers ) ,
				       issueMod , issueAMod );
      f = nf;
      }
     }
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_costs( NCost ,
				     Range( rng.first + f_n_facilities ,
					    rng.second + f_n_facilities ) ,
				     issueMod , issueAMod );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCost , NCost + num , v_t_cost.data() + rng.first );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgTCost ,
			    rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end CapacitatedFacilityLocationBlock::chg_transportation_costs( range )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_costs(
			        c_CV_it NCost , Subset && nms , bool ordered
				ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 // TODO: eliminate from nms the "fake" changes
 const Index maxn = f_n_facilities * f_n_customers;
 if( is_equal( v_f_cost , nms , NCost , maxn ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_t_cost , nms , NCost );

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    Subset nnms( nms );     // copy and translate names
    for( auto & el : nnms )
     el += f_n_facilities;

    // since modify_coefficients owns both vectors, copies has to be made
    CVector NC( NCost , NCost + num );
    get_lfo()->modify_coefficients( std::move( NC ) , std::move( nnms ) ,
				    true , issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    Index f = nms.front();
    Index i = f / f_n_customers;
    Index l = f % f_n_customers;
    if( ( nms.back() / f_n_customers ) == i ) {
     // the subset is all inside a single facility
     Subset nnms( nms );     // copy and translate names
     for( auto & el : nnms )
      el %= f_n_customers;
     KB( v_Block[ i ] )->chg_weighs( NCost , std::move( nnms ) , true ,
				     issueMod , issueAMod );
     break;
     }

    for( auto bit = nms.begin() ; ; ++i ) {
     auto eit = bit;
     while( ( eit != nms.end() ) && ( *eit / f_n_customers ) == i )
      ++eit;

     Subset nnms( bit , eit );     // copy and translate names
     for( auto & el : nnms )
      el %= f_n_customers;

     KB( v_Block[ i ] )->chg_profits( NCost , std::move( nnms ) , true ,
				      issueMod , issueAMod );
     if( eit == nms.end() )
      break;

     NCost += std::distance( bit , eit );
     }
    break;
    }
   default: {  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    Subset nnms( nms );     // copy and translate names
    for( auto & el : nnms )
     el += f_n_facilities;

    MCFB( v_Block[ 1 ] )->chg_costs( NCost , nnms , true ,
				     issueMod , issueAMod );
    }
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_t_cost , nms , NCost );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSnstMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgTCost ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end CapacitatedFacilityLocationBlock::chg_transportation_costs( range )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_cost(
				     CNumber NCost , Index p ,
				     ModParam issueMod , ModParam issueAMod )
{
 const Index maxn = f_n_facilities * f_n_customers;
 if( p >= maxn )
  throw( std::invalid_argument( "invalid name of ( facility , customer ) pair"
				) );

 const Index i = p / f_n_customers;
 const Index j = p % f_n_customers;
 if( v_t_cost[ i ][ j ] == NCost )
  return;

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  v_t_cost[ i ][ j ] = NCost;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    get_lfo()->modify_coefficient( NCost , p + f_n_facilities , issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    KB( v_Block[ i ] )->chg_profit( NCost , j , issueMod , issueAMod );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_cost( NCost , p + f_n_facilities ,
				     issueMod , issueAMod );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   v_t_cost[ i ][ j ] = NCost;

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eChgTCost ,
			     Range( p , p + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_transportation_cost )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_capacities(
				     c_DV_it NCap , Range rng ,
				     ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_facilities );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 const Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NCap , NCap + num , v_capacity.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasCapCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCap , NCap + num , v_capacity.begin() + rng.first );

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod ,
			      ( AR & FormMsk ) != FlwForm ? num : 0 );
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( Index i = rng.first ; i < rng.second ; ++i )
     LF( v_cap[ i ].set_function() )->modify_coefficient( *(NCap++) ,
							  f_n_customers ,
							  iAM );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( Index i = rng.first ; i < rng.second ; ++i )
     KB( v_Block[ i ] )->chg_capacity( *(NCap++) , issueMod , iAM );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_ucaps( NCap , rng , issueMod , iAM );
   }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , num );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NCap , NCap + num , v_capacity.begin() + rng.first );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgCap ,
			    rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end CapacitatedFacilityLocationBlock::chg_facility_capacities( range )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_capacities(
				  c_DV_it NCap , Subset && nms , bool ordered
				  ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_capacity , nms , NCap , f_n_facilities ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasCapCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_capacity , nms , NCap );

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod ,
			      ( AR & FormMsk ) != FlwForm ? num : 0 );
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto i : nms )
     LF( v_cap[ i ].set_function() )->modify_coefficient( *(NCap++) ,
							  f_n_customers ,
							  iAM );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto i : nms )
     KB( v_Block[ i ] )->chg_capacity( *(NCap++) , issueMod , iAM );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_ucaps( NCap , nms , ordered , issueMod , iAM );
   }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , num );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_capacity , nms , NCap );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  // ensure the names are ordered even if they were not so originally
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );

  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgCap ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }
 }  // end CapacitatedFacilityLocationBlock::chg_facility_capacities( sbst )


/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_capacity( Demand NCap ,
			   Index i , ModParam issueMod , ModParam issueAMod )
{
 if( i >= f_n_facilities )
  throw( std::invalid_argument( "invalid facility name" ) );

 if( v_capacity[ i ] == NCap )
  return;

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;

 if( not_dry_run( issueAMod ) && ( AR & HasCapCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  v_capacity[ i ] = NCap;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    LF( v_cap[ i ].set_function() )->modify_coefficient( NCap , f_n_customers ,
							 issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    KB( v_Block[ i ] )->chg_capacity( NCap , issueMod , issueAMod );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_ucap( NCap , i , issueMod , issueAMod );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   v_capacity[ i ] = NCap;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgCap ,
			    Range( i , i + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif
 
 }  // end( CapacitatedFacilityLocationBlock::chg_facility_capacity )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_customers_demands( c_DV_it NDem ,
			 Range rng , ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_customers );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 const Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NDem , NDem + num , v_demand.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasSatCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NDem , NDem + num , v_demand.begin() + rng.first );

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod ,
			      ( AR & FormMsk ) != FlwForm ? num : 0 );
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // if appropriate, open a new channel to bunch up all abstract Modification
    auto iAM = make_amod_param( issueAMod , num );

    for( Index i = rng.first ; i < rng.second ; ++i )
     v_sat[ i ].set_both( *(NCap++) , iAM );

    // if a new channel had been opened, close it
    unmake_amod_param( issueAMod , iAM , num );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    // open a new channel to bunch up all abstract Modification
    auto iAM = make_amod_param( issueAMod , f_n_facilities );

    for( auto bi : v_Block )
     KB( vi )->chg_weights( NCap , rng , issueMod , iAM );

    // close the new channel
    unmake_amod_param( issueAMod , iAM , f_n_facilities );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_dfcts( NCap ,
				     Range( rng.first + f_n_facilities ,
					    rng.second + f_n_facilities ) ,
				     issueMod , issueAMod );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   std::copy( NDem , NDem + num , v_demand.begin() + rng.first );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgDem ,
			    rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_customers_demands( rng ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_customers_demands( c_DV_it NDem ,
			            Subset && nms , bool ordered ,
				    ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_demand , nms , NDem , f_n_customers ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasSatCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_demand , nms , NDem );

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // if appropriate, open a new channel to bunch up all abstract Modification
    auto iAM = make_amod_param( issueAMod , nms.size() );

    for( auto i : nms )
     v_sat[ i ].set_both( *(NDem++) , iAM );
 
    // if a new channel had been opened, close it
    unmake_amod_param( issueAMod , iAM , num );
    break;
   }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    // open a new channel to bunch up all abstract Modification
    auto iAM = make_amod_param( issueAMod , f_n_facilities );

    // note that chg_weights() acquire the vector, so copies must be made
    for( auto bi : v_Block )
     KB( vi )->chg_weights( NDem , Subset( nms ) , true , issueMod , iAM );

    // close the new channel
    unmake_amod_param( issueAMod , iAM , f_n_facilities );
    break;
    }
   default: {  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    Subset nnms( nms.begin() , nms.end() );  // copy and translate nms
    for( auto & el : nnms )
     el+= f_n_facilities;
    MCFB( v_Block[ 1 ] )->chg_dfcts( NDem , nnms , true ,
				     issueMod , issueAMod );
    }
   }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , num );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_demand , nms , NCap );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgDem ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_customers_demands( sbst ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_customers_demand( Demand NDem ,
			  Index j  , ModParam issueMod , ModParam issueAMod )
{
 if( j >= f_n_customers )
  throw( std::invalid_argument( "invalid customer name" ) );

 if( v_demand[ j ] == NDem )
  return;

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;

 if( not_dry_run( issueAMod ) && ( AR & HasSatCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
   v_demand[ j ] = NDem;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    v_sat[ j ].set_both( NDem , issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    // open a new channel to bunch up all abstract Modification
    auto iAM = make_amod_param( issueAMod , f_n_facilities );

    for( auto bi : v_Block )
     KB( vi )->chg_weights( NDem , j , issueMod , iAM );

    // close the new channel
    unmake_amod_param( issueAMod , iAM , f_n_facilities );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_dfct( NDem , j + f_n_facilities ,
				    issueMod , issueAMod );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   v_demand[ j ] = NDem;

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgDem ,
			    Range( i , j + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::chg_customers_demand )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_facilities( Range rng ,
			             ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_facilities );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "close_facilities unavailable if Variable not generated" ) );

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( ( ! get_y( i ).is_fixed() ) || ( get_y( i ).get_value() != 0 ) )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    v_y[ i ].set_value( 0 );
    v_y[ i ].is_fixed( true , iAM );
    }
  else
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 0 );
    KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod>( this ,
			     CapacitatedFacilityLocationBlockMod::eCloseF ,
								     rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_facilities( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_facilities( Subset && nms ,
		      bool ordered , ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "close_facilities unavailable if Variable not generated" ) );

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( ( ! get_y( i ).is_fixed() ) || ( get_y( i ).get_value() != 0 ) )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( auto i : nms ) {
    v_y[ i ].set_value( 0 );
    v_y[ i ].is_fixed( true , iAM );
    }
  else
   for( auto i : nms ) {
    KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 0 );
    KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod>( this ,
			     CapacitatedFacilityLocationBlockMod::eCloseF ,
			     std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_facilities( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::close_facility( Index i ,
			             ModParam issueMod , ModParam issueAMod )
{
 if( i >= f_n_facilities )
 throw( std::invalid_argument( "invalid facility name" ) );

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "close_facilities unavailable if Variable not generated" ) );

 if( get_y( i ).is_fixed() && ( get_y( i ).get_value() == 0 ) )  // fixed
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  if( ( AR & FormMsk ) != KskForm ) {
   v_y[ i ].set_value( 0 );
   v_y[ i ].is_fixed( true , issueAMod );
   }
  else {
   KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 0 );
   KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , issueAMod );
   }
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod>( this ,
			     CapacitatedFacilityLocationBlockMod::eCloseF ,
			     Range( i , i + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::close_facility )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_facilities( Range rng ,
			             ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_facilities );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "open_facilities unavailable if Variable not generated" ) );

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( get_y( i ).is_fixed() )
   ++cnt;

 if( ! cnt )  // all facilities are open already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( Index i = rng.first ; i < rng.second ; ++i )
    v_y[ i ].is_fixed( false , iAM );
  else
   for( Index i = rng.first ; i < rng.second ; ++i )
    KB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , iAM );

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod>( this ,
			     CapacitatedFacilityLocationBlockMod::eOpenF ,
								     rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_facilities( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_facilities( Subset && nms ,
		      bool ordered , ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	        "open_facilities unavailable if Variable not generated" ) );

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( get_y( i ).is_fixed() )
   ++cnt;

 if( ! cnt )  // all facilities are open already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( auto i : nms )
    v_y[ i ].is_fixed( true , iAM );
  else
   for( auto i : nms )
    KB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , iAM );

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod>( this ,
			     CapacitatedFacilityLocationBlockMod::eOpenF ,
			     std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_facilities( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::open_facility( Index i ,
			             ModParam issueMod , ModParam issueAMod )
{
 if( i >= f_n_facilities )
 throw( std::invalid_argument( "invalid facility name" ) );

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "open_facilities unavailable if Variable not generated" ) );

 if( ! get_y( i ).is_fixed() )  // fixed already
  return;                       // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  if( ( AR & FormMsk ) != KskForm )
   v_y[ i ].is_fixed( false , issueAMod );
  else
   KB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , issueAMod );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eCloseF ,
			     Range( i , i + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::open_facility )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::fix_open_facilities( Range rng ,
			             ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_facilities );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	     "fix_open_facilities unavailable if Variable not generated" ) );

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( ( ! get_y( i ).is_fixed() ) || ( get_y( i ).get_value() != 1 ) )
   ++cnt;

 if( ! cnt )  // all facilities are fixed open already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    v_y[ i ].set_value( 1 );
    v_y[ i ].is_fixed( true , iAM );
    }
  else
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 1 );
    KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eBuyF ,
								     rng ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::fix_open_facilities( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::fix_open_facilities( Subset && nms ,
		      bool ordered , ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	     "fix_open_facilities unavailable if Variable not generated" ) );

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( ( ! get_y( i ).is_fixed() ) || ( get_y( i ).get_value() != 0 ) )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  auto iAM = make_amod_param( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( auto i : nms ) {
    v_y[ i ].set_value( 1 );
    v_y[ i ].is_fixed( true , iAM );
    }
  else
   for( auto i : nms ) {
    KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 1 );
    KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }

  // if a new channel had been opened, close it
  unmake_amod_param( issueAMod , iAM , cnt );
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eBuyF ,
			     std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::fix_open_facilities( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::fix_open_facility( Index i ,
			             ModParam issueMod , ModParam issueAMod )
{
 if( i >= f_n_facilities )
 throw( std::invalid_argument( "invalid facility name" ) );

 if( ! ( AR & HasVar ) )
  throw( std::logic_error(
	       "close_facilities unavailable if Variable not generated" ) );

 if( get_y( i ).is_fixed() && ( get_y( i ).get_value() ==  ) )  // fixed open
  return;     // nothing to do

 // since the physical and abstract representation are the same, anything
 // that has to do with the abstract representation is skipped in the
 // "dry run" case; but the "phisical Modification" is issued anyway
 if( not_dry_run( issueAMod ) ) {

  if( ( AR & FormMsk ) != KskForm ) {
   v_y[ i ].set_value( 1 );
   v_y[ i ].is_fixed( true , issueAMod );
   }
  else {
   KB( v_Block[ i ] )->get_Var( f_n_customers ).set_value( 1 );
   KB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , issueAMod );
   }
  }

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eBuyF ,
			     Range( i , i + 1 ) ) ,
			   Observer::par2chnl( issueMod ) );
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::fix_open_facility )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::print( std::ostream & output ) const
{
 if( verbosity_lvl != Block::complete ) {  // non-complete version
  // only basic information 
  output << "CapacitatedFacilityLocationBlock with: " << f_n_facilities
	 << " facilities and " << f_n_cistomers() << " customers"
	 << std::endl;
  }
 else  {
  // complete version: file in standard ORLib format
  output << f_n_facilities << std::endl;
  output << f_n_customers << std::endl << std::endl;

  for( Index i = 0 ; i < f_n_facilities ; ++i )
   output << v_capacity[ i ] << "\t" << v_f_cost[ i ] << std::endl;
  
  output << std::endl;

  for( Index j = 0 ; i < f_n_customers ; ++i ) {
   output << v_demand[ j ] << std::endl;
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    output << v_t_cost[ i ][ j ] << "\t";
   output << std::endl;
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
 f_obj.clear();

 // delete all Variable
 x.clear();
 y.clear();

 // delete all sub-Block
 for( auto bi : v_Block )
  delete bi;

 v_Block.clear();  // then clear the vector

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

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationSF(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock in the "Standard
  * Formulation" is a "leaf" Block (has no sub-Block), this method does not
  * have to deal with GroupModification since these are produced by
  * Block::add_Modification(), but this method is called *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN CapacitatedFacilityLocationBlock
  *   WHEN THIS METHOD IS EXECUTED IS PRECISELY THE ONE IN WHICH THE
  *   Modification WAS ISSUED
  *
  * This assumption drastically simplifies some of the logic here. */

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinRngd * >( mod ) ) {
  Index f = tmod->range().first;
  const Index s = tmod->range().second;

  auto lfo = LF( tmod->function() );
  if( LF( f_obj.get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   if( f < f_n_facilities ) {  // facility costs are modified
    c_Index end = std::min( s , f_n_facilities );
    c_Index sz = end - f;
    if( sz == 1 )              // one facility
     chg_facility_cost( (lfo->get_v_var())[ f ].second , f ,
			make_par( eNoBlck , chnl ) , eDryRun );
    else {                     // many facilities
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < end ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_facility_costs( NC.begin() , Range( f , end ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
     }

    f = f_n_facilities;  // facilities costs accounted for
    }

   if( s >= f_n_facilities ) {  // transportation costs are modified
    c_Index sz = s - f;
    if( sz == 1 )               // one pair
     chg_transportation_cost( (lfo->get_v_var())[ f ].second ,
			      f - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_transportation_costs( NC.begin() , Range( f - f_n_facilities ,
						   s - f_n_facilities ) ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   return;
   }

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_observer() );
  if( ! cnst )
   throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

  if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
   throw( std::invalid_argument(
	    "        Modification to FRowConstraint not capacity one" ) );

  if( ( f != f_n_customers ) || ( s != f + 1 ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  // note that the coefficient of y is the opposite of the capacity
  chg_facility_capacity( - (lfo->get_v_var())[ f ].second ,
			 std::distance( & v_cap.front() , cap ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinSbst * >( mod ) ) {
  auto & nms = tmod->subset();

  auto lfo = LF( tmod->function() );
  if( LF( f_obj.get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   auto nmsit = nms.begin();

   if( nms.front() < f_n_facilities ) {  // facility costs are modified
    while( ( nmsit != nms.end() ) && ( *nmsit < f_n_facilities ) )
     ++nmsit;

    c_Index sz = std::distance( nms.begin() , nmsit );
    if( sz == 1 ) {            // one facility
     auto i = *std::prev( nmsit );
     chg_facility_cost( (lfo->get_v_var())[ i ].second , i ,
			make_par( eNoBlck , chnl ) , eDryRun );
     }
    else {                     // many facilities
     CVector NC( sz );
     Subset nnms( sz );
     auto NCit = NC.begin();
     auto nnmsit = nnms.begin();
     for( auto it = nms.begin() ; it != nmsit ; ) {
      auto i = *(it++);
      *(nnmsit++) = i;
      *(NCit++) = (lfo->get_v_var())[ i ].second;
      }
     chg_facility_costs( NC.begin() , std::move( nnms ) , true ,
			 make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   if( nms.back() >= f_n_facilities ) {  // transportation costs are modified
    c_Index sz = std::distance( nmsit , nms.end() );
    if( sz == 1 )               // one pair
     chg_transportation_cost( (lfo->get_v_var())[ *nmsit ].second ,
			      *nmsit - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( sz );
     Subset nnms( sz );
     auto NCit = NC.begin();
     auto nnmsit = nnms.begin();
     for( ; nmsit != nms.end() ; ++nmsit ) {
      auto h = *(it++);
      *(nnmsit++) = h - f_n_facilities;
      *(NCit++) = (lfo->get_v_var())[ h ].second;
      }
     chg_transportation_costs( NC.begin() , std::move( nnms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   return;
   }

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_observer() );
  if( ! cnst )
   throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

  if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
   throw( std::invalid_argument(
	            "Modification to FRowConstraint not capacity one" ) );

  if( ( nms.size() != 1 ) || ( nms.front() != f_n_customers ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  // note that the coefficient of y is the opposite of the capacity
  chg_facility_capacity( - (lfo->get_v_var())[ f_n_customers ].second ,
			 std::distance( & v_cap.front() , cap ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }


 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const RowConstraintMod * >( mod ) )
   throw( std::invalid_argument( "RowConstraintMod not allowed" ) );

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const VariableMod * >( mod ) )
  throw( std::invalid_argument( "VariableMod not allowed" ) );

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationSF )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationKF(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock in the "Knapsack
  * Formulation" is *not* a "leaf" Block, i.e., it has sub-Block, this
  * method must deal with GroupModification, since these can be produced by
  * Block::add_Modification() in the sub-Block. In fact, this method only
  * does that, and then relies on a further guts_of_guts to do the actual
  * processing. This introduces delay between the moment in which the
  * Modification is produced and the one in which it is processed, which
  * would in principle complicate the logic. However
  *
  *     CapacitatedFacilityLocationBlock IS A "STATIC" Block IN WHICH THE
  *     SIZE OF THE STUFF NEVER CHANGES (save if it is re-loaded whole)
  *
  * This means that the indices, sanges and subsets found in the Modification
  * are always still valid, which drastically simplifies some of the logic.
  */

 // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
  for( auto submod : tmod->sub_Modifications() )
   guts_of_add_ModificationFF( submod , chnl );

  return;
  }
 
 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the only possible changes are those of the Objective and of the
 // capacity (opposite of weight of last item) in some BinaryKnapsackBlock,
 // the only other Function in this formulation are in the customers
 // satisfaction constraints that cannot be changed

 if( auto tmod = dynamic_cast< const C05FunctionModLinRngd * >( mod ) ) {
  // must be coming from a BinaryKnapsackBlock
  auto bkb = dynamic_cast< BinaryKnapsackBlock * >( mod->get_Block() );
  if( ! bkb )
   throw( std::invalid_argument( "unsupportd C05FunctionModLinRngd" ) );

  // find which one
  auto it = std::lower_bound( v_Block.begin() , v_Block.end() , bkb );
  c_Index i = std::distance( v_Block.begin() , it ); 
  c_Index f = tmod->range().first;
  c_Index s = tmod->range().second;

  auto lfo = LF( tmod->function() );
  if( LF( static_cast< FRealObjective * >( bkb.get_objective()
					   )->get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   if( f < f_n_customers ) {  // transportation costs are modified
    c_Index end = std::min( s , f_n_customers )
    c_Index sz = s - f;
    c_Index offst = f_n_customers * i;
    if( sz == 1 )               // one pair
     chg_transportation_cost( (lfo->get_v_var())[ f ].second , f + offst ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < end ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_transportation_costs( NC.begin() ,
			       Range( f + offst , end + offst ) ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   if( s > f_n_customers )  // the facility cost is modified
    chg_facility_cost( (lfo->get_v_var())[ f_n_customers ].second , i ,
		       make_par( eNoBlck , chnl ) , eDryRun );

   return;
   }

  // the Modificatin can only come from the only FRowConstraint
  if( ( f != f_n_customers ) || ( s != f + 1 ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  // note that the coefficient of y is the opposite of the capacity
  chg_facility_capacity( (lfo->get_v_var())[ f ].second , i ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the only possible changes are those of the Objective and of the
 // capacity (opposite of weight of last item) in some BinaryKnapsackBlock,
 // the only other Function in this formulation are in the customers
 // satisfaction constraints that cannot be changed
 if( auto tmod = dynamic_cast< const C05FunctionModLinSbst * >( mod ) ) {
  // find which one
  auto it = std::lower_bound( v_Block.begin() , v_Block.end() , bkb );
  c_Index i = std::distance( v_Block.begin() , it ); 
  auto & nms = tmod->subset();

  auto lfo = LF( tmod->function() );
  if( LF( static_cast< FRealObjective * >( bkb.get_objective()
					   )->get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   if( nms.front() < f_n_customers ) {   // transportation costs are modified
    c_Index offst = f_n_customers * i;
    Subset nnms( nms.begin() , nms.back() == f_n_customers ?
		               std::prev( nms.end ) : nms.end() );
    if( nms.size() == 1 )               // only one pair
     chg_transportation_cost( (lfo->get_v_var())[ nms.front() ].second ,
			      nnms.front() + offst ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( nnms.size() );
     auto NCit = NC.begin();
     for( auto & j : nnms ) {
      *(NCit++) = (lfo->get_v_var())[ j ].second;
      j += offst;
      }
     chg_transportation_costs( NC.begin() , Subset( nnms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   if( nms.back() == f_n_customers )  // the facility cost is modified
    chg_facility_cost( (lfo->get_v_var())[ f_n_customers ].second , i ,
		       make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  // the Modificatin can only come from the only FRowConstraint
  if( ( nms.size() != 1 ) || ( nms.front() != f_n_customers ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  chg_facility_capacity( (lfo->get_v_var())[ f_n_customers ].second , i ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const RowConstraintMod * >( mod ) )
   throw( std::invalid_argument( "RowConstraintMod not allowed" ) );

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const VariableMod * >( mod ) )
  throw( std::invalid_argument( "VariableMod not allowed" ) );

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationKF )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationFF(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock in the "Flow
  * Formulation" is *not* a "leaf" Block, i.e., it has sub-Block, this
  * method must deal with GroupModification, since these can be produced by
  * Block::add_Modification() in the sub-Block. In fact, this method only
  * does that, and then relies on a further guts_of_guts to do the actual
  * processing. This introduces delay between the moment in which the
  * Modification is produced and the one in which it is processed, which
  * would in principle complicate the logic. However
  *
  *     CapacitatedFacilityLocationBlock IS A "STATIC" Block IN WHICH THE
  *     SIZE OF THE STUFF NEVER CHANGES (save if it is re-loaded whole)
  *
  * This means that the indices, sanges and subsets found in the Modification
  * are always still valid, which drastically simplifies some of the logic.
  */

 // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
  for( auto submod : tmod->sub_Modifications() )
   guts_of_add_ModificationFF( submod , chnl );

  return;
  }
 
 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinRngd * >( mod ) ) {
  Index f = tmod->range().first;
  const Index s = tmod->range().second;

  auto lfo = LF( tmod->function() );
  if( LF( f_obj.get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   if( f < f_n_facilities ) {  // facility costs are modified
    c_Index end = std::min( s , f_n_facilities );
    c_Index sz = end - f;
    if( sz == 1 )              // one facility
     chg_facility_cost( (lfo->get_v_var())[ f ].second , f ,
			make_par( eNoBlck , chnl ) , eDryRun );
    else {                     // many facilities
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < end ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_facility_costs( NC.begin() , Range( f , end ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
     }

    f = f_n_facilities;  // facilities costs accounted for
    }

   if( s >= f_n_facilities ) {  // transportation costs are modified
    c_Index sz = s - f;
    if( sz == 1 )               // one pair
     chg_transportation_cost( (lfo->get_v_var())[ f ].second ,
			      f - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_transportation_costs( NC.begin() , Range( f - f_n_facilities ,
						   s - f_n_facilities ) ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   return;
   }

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_observer() );
  if( ! cnst )
   throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

  if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
   throw( std::invalid_argument(
	    "        Modification to FRowConstraint not capacity one" ) );

  if( ( f != f_n_customers ) || ( s != f + 1 ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  // note that the coefficient of y is the opposite of the capacity
  chg_facility_capacity( - (lfo->get_v_var())[ f ].second ,
			 std::distance( & v_cap.front() , cap ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinSbst * >( mod ) ) {
  auto & nms = tmod->subset();

  auto lfo = LF( tmod->function() );
  if( LF( f_obj.get_function() ) == lfo ) {
   // Modification to the Objective - - - - - - - - - - - - - - - - - - - - -

   auto nmsit = nms.begin();

   if( nms.front() < f_n_facilities ) {  // facility costs are modified
    while( ( nmsit != nms.end() ) && ( *nmsit < f_n_facilities ) )
     ++nmsit;

    c_Index sz = std::distance( nms.begin() , nmsit );
    if( sz == 1 ) {            // one facility
     auto i = *std::prev( nmsit );
     chg_facility_cost( (lfo->get_v_var())[ i ].second , i ,
			make_par( eNoBlck , chnl ) , eDryRun );
     }
    else {                     // many facilities
     CVector NC( sz );
     Subset nnms( sz );
     auto NCit = NC.begin();
     auto nnmsit = nnms.begin();
     for( auto it = nms.begin() ; it != nmsit ; ) {
      auto i = *(it++);
      *(nnmsit++) = i;
      *(NCit++) = (lfo->get_v_var())[ i ].second;
      }
     chg_facility_costs( NC.begin() , std::move( nnms ) , true ,
			 make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   if( nms.back() >= f_n_facilities ) {  // transportation costs are modified
    c_Index sz = std::distance( nmsit , nms.end() );
    if( sz == 1 )               // one pair
     chg_transportation_cost( (lfo->get_v_var())[ *nmsit ].second ,
			      *nmsit - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( sz );
     Subset nnms( sz );
     auto NCit = NC.begin();
     auto nnmsit = nnms.begin();
     for( ; nmsit != nms.end() ; ++nmsit ) {
      auto h = *(it++);
      *(nnmsit++) = h - f_n_facilities;
      *(NCit++) = (lfo->get_v_var())[ h ].second;
      }
     chg_transportation_costs( NC.begin() , std::move( nnms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   return;
   }

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_observer() );
  if( ! cnst )
   throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

  if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
   throw( std::invalid_argument(
	            "Modification to FRowConstraint not capacity one" ) );

  if( ( nms.size() != 1 ) || ( nms.front() != f_n_customers ) )
   throw( std::invalid_argument(
	   "Modification to wrong coefficient in capacity constraint" ) );

  // note that the coefficient of y is the opposite of the capacity
  chg_facility_capacity( - (lfo->get_v_var())[ f_n_customers ].second ,
			 std::distance( & v_cap.front() , cap ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }


 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const RowConstraintMod * >( mod ) )
   throw( std::invalid_argument( "RowConstraintMod not allowed" ) );

 /*!!
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
  !!*/

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const VariableMod * >( mod ) )
  throw( std::invalid_argument( "VariableMod not allowed" ) );

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationSF )

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
  if( v_f_cost[ i ] >= 0 )
   f_cond_upper += v_f_cost[ i ];
  else
   f_cond_lower += v_f_cost[ i ];

 for( Index j = 0 ; j < f_n_customers ; ++j ) {
  auto minj = Inf< TCost >();
  auto maxj = - Inf< TCost >();

  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   if( minj > v_t_cost[ j ][ i ] )
    minj = v_t_cost[ j ][ i ];
   if( maxj < v_t_cost[ j ][ i ] )
    maxj = v_t_cost[ j ][ i ];
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
 if( ! num )
  return( issueAMod );

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
 if( ( newiAM == eNoMod ) || ( ! num ) )
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

 auto fs = group.getVar( "FacilitySolution" );
 if( fs.isNull() )
  v_y.clear();
 else {
  auto nf = group.getDim( "NFacilities" );
  if( nf.isNull() )
   throw( std::invalid_argument( _pfrx + "NFacilities dimension required"
				 ) );
  v_y.resize( nf.getSize() );
  fs.getVar( v_y.data() );
  }

 auto ts = group.getVar( "TransportationSolution" );
 if( ts.isNull() )
  v_x.clear();
 else {
  auto td = group.getDim( "TransportationDim" );
  if( td.isNull() )
   throw( std::invalid_argument( _pfrx +
				 "TransportationDim dimension required" ) );
  v_x.resize( td.getSize() );
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
  if( v_x.size() < CFLB->get_NFacilities() * CFLB->get_NCustomers() )
   v_x.resize( CFLB->get_NFacilities() * CFLB->get_NCustomers() );

  CFLB->get_transportation_solution( v_x.data().begin() );
  }
 }  // end( CapacitatedFacilityLocationSolution::read )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::write( Block * block ) 
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
  if( v_x.size() < CFLB->get_NFacilities() * CFLB->get_NCustomers() )
   throw( std::invalid_argument( "incompatible transportation size" ) );

  CFLB->set_transportation_solution( v_x.data().begin() );
  }
 }  // end( CapacitatedFacilityLocationSolution::write )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::serialize( netCDF::NcGroup & group )
{
 if( ! v_y.empty() ) {
  auto nf = group.addDim( "NFacilities" , v_y.size() );
  ::serialize( group , "FacilitySolution" , netCDF::NcDouble() , nf ,
	       v_y.data() );
  }

 if( ! v_x.empty() ) {
  auto td = group.addDim( "TransportationDim" , v_x.size() );
  ::serialize( group , "TransportationSolution" , netCDF::NcDouble() , td ,
	       v_x.data() );
  }
 }  // end( CapacitatedFacilityLocationSolution::serialize )

/*--------------------------------------------------------------------------*/

CapacitatedFacilityLocationSolution *
            CapacitatedFacilityLocationSolution::scale( double factor ) const
{
 auto * sol = CapacitatedFacilityLocationSolution::clone( true );

 if( ! v_y.empty() )
  for( Block::Index i = 0 ; i < v_y.size() ; ++i )
   sol->v_y[ i ] = v_y[ i ] * factor;

 if( ! v_x.empty() )
  for( Block::Index i = 0 ; i < v_x.size() ; ++i )
   sol->v_x[ i ] = v_x[ i ] * factor;

 return( sol );

 }  // end( CapacitatedFacilityLocationSolution::scale )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::sum( const Solution * solution ,
					       double multiplier )
{
 auto CFLS = dynamic_cast< const CapacitatedFacilityLocationSolution * >(
								  solution );
 if( ! CFLS )
  throw( std::invalid_argument(
		 "solution is not a CapacitatedFacilityLocationSolution" ) );

 if( ! v_y.empty() ) {
  if( v_y.size() != CFLS->v_y.size()  )
   throw( std::invalid_argument( "incompatible facility size" ) );

  auto yit = CFLS->v_y.begin();
  for( auto & yi : v_y )
   yi = *(yit++) * multiplier;
  }

 if( ! v_x.empty() ) {
  if( v_x.size() != CFLS->v_x.size() )
   throw( std::invalid_argument( "incompatible transportation size" ) );

  auto xit = CFLS->v_x.begin();
  for( auto & xi : v_x )
   xi = *(xit++) * multiplier;
  }
 }  // end( CapacitatedFacilityLocationSolution::sum )

/*--------------------------------------------------------------------------*/

CapacitatedFacilityLocationSolution *
               CapacitatedFacilityLocationSolution::clone( bool empty ) const
{
 auto *sol = new CapacitatedFacilityLocationSolution();

 if( empty ) {
  if( ! v_y.empty() )
   sol->v_y.resize( v_pi.size() );

  if( ! v_x.empty() )
   sol->v_x.resize( v_x.size() );
  }
 else {
  sol->v_y = v_y;
  sol->v_x = v_x;
  }

 return( sol );

 }  // end( CapacitatedFacilityLocationSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------- End File CapacitatedFacilityLocationBlock.cpp --------------*/
/*--------------------------------------------------------------------------*/
