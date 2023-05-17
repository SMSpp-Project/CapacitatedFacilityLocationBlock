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
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CapacitatedFacilityLocationBlock.h"

#include "AbstractBlock.h"

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
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

using v_coeff_pair = LinearFunction::v_coeff_pair;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr unsigned char FormMsk = 3;
// mask for removing all but the first two bits and only leaving the
// formulation (irrespective of if it is splittable or not)

static constexpr unsigned char StdForm = 0;
// the "standard" formulation is used

static constexpr unsigned char KskForm = 1;
// the "knapsack" formulation is used

static constexpr unsigned char FlwForm = 2;
// the "flow" formulation is used

static constexpr unsigned char UnSpltF = 4;
// 3rd bit of AR == 1 if the problem is unsplittable (the X[] are integer)

static constexpr unsigned char HasVar = 8;
// 4th bit of AR == 1 if the Variable have been constructed

static constexpr unsigned char HasObj = 16;
// 5th bit of AR == 1 if the Objective has been constructed

static constexpr unsigned char HasSatCns = 32;
// 6th bit of AR == 1 if the customer satisfaction Constraints are constructed

static constexpr unsigned char HasCapCns = 64;
// 7th bit of AR == 1 if the capacity Constraints are constructed

static constexpr unsigned char HasStrngCns = 128;
// 8th bit of AR == 1 if the strong Linking Constraints are separated

static constexpr unsigned char yFree = 0;  // facility is free

static constexpr unsigned char yFxd0 = 1;  // facility is closed

static constexpr unsigned char yFxd1 = 2;  // facility is fixed-open

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
// returns true if two vectors differ, one of them being given as a pointer
// to an array and a subset of indices

template< typename T >
static bool is_equal( T * vec , Block::c_Subset & nms ,
		      typename std::vector< T >::const_iterator cmp ,
		      Block::Index n_max )
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
static void copyidx( T * vec , Block::c_Subset & nms ,
		     typename std::vector< T >::const_iterator cpy )
{
 for( auto nm : nms )
  *( vec + nm ) = *(cpy++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register CapacitatedFacilityLocationBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( CapacitatedFacilityLocationBlock );

// register CapacitatedFacilityLocationSolution to the Solution factory

SMSpp_insert_in_factory_cpp_1( CapacitatedFacilityLocationSolution );

/*--------------------------------------------------------------------------*/
/*--------------- METHODS OF CapacitatedFacilityLocationBlock --------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::load( Index m , Index n ,
					     DVector && Q , CVector && F ,
					     DVector && D , CMatrix && C ,
					     bool unsplt )
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
 // this resize() should most definitely not be necessary, but it appears
 // that "=" between multi_arrays does not work as intended
 v_t_cost.resize( boost::extents[ m ][ n ] );
 v_t_cost = std::move( C );

 v_fxd.resize( m );
 std::fill( v_fxd.begin() , v_fxd.end() , yFree );

 f_cond_lower = f_cond_upper = dNaN;  // reset conditional bounds

 f_unsplittable = unsplt;

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 }  // end( CapacitatedFacilityLocationBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::load( std::istream & input ,
					     char frmt )
{
 static const std::string _prfx = "CapacitatedFacilityLocationBlock::load: ";

 // erase existing abstract representation, if any - - - - - - - - - - - - - -

 if( AR & ~7 )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -
 // the first part of the three formats at least is common

 input >> eatcomments >> f_n_facilities;
 if( input.fail() )
  goto input_failure;

 if( f_n_facilities == 0 )
  throw( std::invalid_argument( _prfx + "number of facilities too small" ) );

 input >> eatcomments >> f_n_customers;
 if( input.fail() )
  goto input_failure;

 if( f_n_customers == 0 )
  throw( std::invalid_argument( _prfx + "number of customers too small" ) );

 // size up internal data structures- - - - - - - - - - - - - - - - - - - - -
 v_capacity.resize( f_n_facilities );
 v_f_cost.resize( f_n_facilities );
 v_demand.resize( f_n_customers );
 v_t_cost.resize( boost::extents[ f_n_facilities ][ f_n_customers ] );
 v_fxd.resize( f_n_facilities );
 std::fill( v_fxd.begin() , v_fxd.end() , yFree );

 // now the format-specific parts
 switch( std::toupper( frmt ) ) {
  case( 'F' ):  // facility oriented, demands-first format- - - - - - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for( Index j = 0 ; j < f_n_customers ; ++j ) {  // read demands
    input >> eatcomments >> v_demand[ j ];
    if( input.fail() )
     goto input_failure;
    if( v_demand[ j ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive demand" ) );
    }

   for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // read capacity
    input >> eatcomments >> v_capacity[ i ];
    if( input.fail() )
     goto input_failure;
    if( v_capacity[ i ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive capacity" ) );
    }

   for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // read fixed cost
    input >> eatcomments >> v_f_cost[ i ];
    if( input.fail() )
     goto input_failure;
    }

   for( Index i = 0 ; i < f_n_facilities ; ++i )
    for( Index j = 0 ; j < f_n_customers ; ++j ) {
     input >> v_t_cost[ i ][ j ];
     if( input.fail() )
      goto input_failure;
     // these are unitary costs, make them total costs
     v_t_cost[ i ][ j ] *= v_demand[ j ];
     }

   break;

  case( 'L' ):  // facility oriented, demands-last format - - - - - - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
    input >> eatcomments >> v_capacity[ i ];
    if( input.fail() )
     goto input_failure;
    if( v_capacity[ i ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

    input >> eatcomments >> v_f_cost[ i ];
    if( input.fail() )
     goto input_failure;

    }  // end( for( each facility ) )

   for( Index j = 0 ; j < f_n_customers ; ++j ) {  // read demands
    input >> eatcomments >> v_demand[ j ];
    if( input.fail() )
     goto input_failure;
    if( v_demand[ j ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive demand" ) );
    }

   for( Index i = 0 ; i < f_n_facilities ; ++i )
    for( Index j = 0 ; j < f_n_customers ; ++j ) {
     input >> v_t_cost[ i ][ j ];
     if( input.fail() )
      goto input_failure;
     }

   break;

  default:  // the ORLib, customers-oriented format - - - - - - - - - - - - -
            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
    input >> eatcomments >> v_capacity[ i ];
    if( input.fail() )
     goto input_failure;
    if( v_capacity[ i ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive capacity" ) );

    input >> eatcomments >> v_f_cost[ i ];
    if( input.fail() )
     goto input_failure;

    }  // end( for( each facility ) )

   for( Index j = 0 ; j < f_n_customers ; ++j ) {  // for( each customer )
    input >> eatcomments >> v_demand[ j ];
    if( input.fail() )
     goto input_failure;
    if( v_demand[ j ] <= 0 )
     throw( std::invalid_argument( _prfx + "non-positive demand" ) );

    for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for( each facility )
     input >> eatcomments >> v_t_cost[ i ][ j ];
     if( input.fail() )
      goto input_failure;

     }  // end( for( each facility ) )
    }  // end( for( each customer ) )

  }  // end( switch( frmt ) ) - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 f_cond_lower = f_cond_upper = dNaN;  // reset conditional bounds
 f_unsplittable = false;

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

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
 fc.getVar( v_f_cost.data() );

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

 v_fxd.resize( f_n_facilities );
 auto ff = group.getVar( "FacilityFix" );
 if( ff.isNull() )
  std::fill( v_fxd.begin() , v_fxd.end() , yFree );
 else {
  auto ffs = ::get_sizes_dimensions( ff );
  if( ( ffs.size() != 1 ) || ( ffs[ 0 ] != f_n_facilities ) )
   throw( std::logic_error( _prfx + "FacilityFix has wrong size" ) );
  ff.getVar( v_fxd.data() );
  }

 f_unsplittable = false;

 auto unsplt = group.getDim( "UnSplittable" );
 if( ! unsplt.isNull() )
  f_unsplittable = ( unsplt.getSize() > 0 );

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

 Index wf = 0;
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  wf = sci->f_value;

 f_unsplittable = wf & UnSpltF;

 if( ! ( wf & FormMsk ) ) {  // "natural formulation" (NF)- - - - - - - - - -
                             // - - - - - - - - - - - - - - - - - - - - - - -
  // AR |= StdForm;  does nothing
  v_y.resize( f_n_facilities );
  auto fxdit = v_fxd.begin();
  for( auto & yi : v_y ) {
   yi.set_type( ColVariable::kBinary );
   if( auto fi = *(fxdit++) ) {
    yi.set_value( fi == yFxd0 ? 0 : 1 );
    yi.is_fixed( true , eNoMod );
    }
   }
  add_static_variable( v_y , "y" );

  v_x.resize( boost::extents[ f_n_facilities ][ f_n_customers ] );
  auto xt = ColVariable::kPosUnitary;
  if( f_unsplittable ) {
   AR |= UnSpltF;
   xt = ColVariable::kBinary;
   }

  Index cnt = f_n_facilities * f_n_customers;
  for( auto xij = v_x.data() ; xij != v_x.data() + cnt ; )
   (xij++)->set_type( xt );
  add_static_variable( v_x , "x" );

  return;
  }

 if( ( wf & FormMsk ) == KskForm ) {  // "knapsack formulation" (KF)- - - - -
                                      //- - - - - - - - - - - - - - - - - - -
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
  BinaryKnapsackBlock::doubleVec P( f_n_customers + 1 );
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
    P[ j ] = v_t_cost[ i ][ j ];
    }
   W[ f_n_customers ] = - v_capacity[ i ];
   P[ f_n_customers ] = v_f_cost[ i ];

   auto bi = BKB( v_Block[ i ] );
   bi->load( f_n_customers + 1 , 0 , W , P , I );
   if( v_fxd[ i ] != yFree )
    bi->fix_x(  v_fxd[ i ] == yFxd1 , i , eNoMod , eNoMod );    
   bi->set_objective_sense( false , eNoMod , eNoMod );
   bi->generate_abstract_variables();
   }

  return;
  }

 if( ( wf & FormMsk ) >= FlwForm ) {  // "flow formulation" (FF)- - - - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  if( f_unsplittable )
   throw( std::invalid_argument(
	   "unsplittable problem not supported with the Flow Formulation" ) );
  
  AR |= FlwForm;

  v_Block.resize( 2 );  // exactly two sub-Block

  // the first sub-Block is an AbstractBlock with the y[] variables
  auto ab = new AbstractBlock( this );
  v_Block[ 0 ] = ab;

  v_y.resize( f_n_facilities );
  auto fxdit = v_fxd.begin();
  for( auto & yi : v_y ) {
   yi.set_type( ColVariable::kBinary );
   if( auto fi = *(fxdit++) ) {
    yi.set_value( fi == yFxd0 ? 0 : 1 );
    yi.is_fixed( true , eNoMod );
    }
   }
  ab->add_static_variable( v_y , "y" );

  // the second Block is a MCFBlock as constructed by get_R3_Block
  SimpleConfiguration< int > r3bc( ( wf & FormMsk ) - 1 );
  auto mcfb = static_cast< MCFBlock * >( get_R3_Block( & r3bc ) );
  v_Block[ 1 ] = mcfb;

  // ... except the cost of the facility arcs are zeros
  MCFBlock::Vec_CNumber zero( f_n_facilities , 0 );
  mcfb->chg_costs( zero.begin() , Range( 0 , f_n_facilities ),
		   eNoMod , eNoMod );
  mcfb->generate_abstract_variables();
  mcfb->set_f_Block( this );

  return;
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
  throw( std::logic_error( _prfx + "generate_abstract_variables not called" ) );

 Index wc = 3;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  wc = sci->f_value;

 if( ( AR & FormMsk ) == StdForm ) {  // "natural formulation" (NF) - - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  if( ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) {
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

  if( ( wc & 2 ) && ( ! ( AR & HasCapCns ) ) ) {
   // construct facility capacity constraints (not there already)
   v_cap.resize( f_n_facilities );
   for( Index i = 0 ; i < f_n_facilities ; ++i ) {
    v_coeff_pair coeffs( f_n_customers + 1 );

    for( Index j = 0 ; j < f_n_customers ; ++j )
     coeffs[ j ] = std::make_pair( & v_x[ i ][ j ] , v_demand[ j ] );

    coeffs[ f_n_customers ] = std::make_pair( & v_y[ i ] ,
					      - v_capacity[ i ] );
    v_cap[ i ].set_rhs( 0 );
    v_cap[ i ].set_lhs( -Inf< RowConstraint::RHSValue >() );
    v_cap[ i ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
    }

   add_static_constraint( v_cap , "cap" );
   AR |= HasCapCns;
   }

  goto Strong_Linking;
  }

 if( ( AR & FormMsk ) == KskForm ) {  // "knapsack formulation" (KF)- - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  if( ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) {
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

  if( ( wc & 2 ) && ( ! ( AR & HasCapCns ) ) ) {
   // construct facility capacity constraints (not there already)
   // these are just the constraint in the BinaryKnapsackBlock
   for( auto ki : v_Block )
    ki->generate_abstract_constraints();

   AR |= HasCapCns;
   }

  goto Strong_Linking;
  }

 // else it is the "flow formulation" (FF)- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wc & 1 ) && ( ! ( AR & HasSatCns ) ) ) {
  // construct customer satisfaction constraints (not there already)
  // these are, just the MCFBlock ones
  v_Block[ 1 ]->generate_abstract_constraints();
  AR |= HasSatCns;
  }

 if( ( wc & 2 ) && ( ! ( AR & HasCapCns ) ) ) {
  // construct facility capacity constraints (not there already)
  // these are the ones linking the y[] variables to the
  // arc flow variables of facility arcs

  auto ab = AB( v_Block[ 0 ] );
  auto mcfb = MCFB( v_Block[ 1 ] );

  v_cap.resize( f_n_facilities );
  for( Index i = 0 ; i < f_n_facilities ; ++i ) {
   v_coeff_pair coeffs( 2 );

   coeffs[ 0 ] = std::make_pair( mcfb->i2p_x( i ) , double( 1 ) );
   coeffs[ 1 ] = std::make_pair( & v_y[ i ] , - v_capacity[ i ] );

   v_cap[ i ].set_rhs( 0 );
   v_cap[ i ].set_lhs( -Inf< RowConstraint::RHSValue >() );
   v_cap[ i ].set_function( new LinearFunction( std::move( coeffs ) , 0 ) );
   }

  ab->add_static_constraint( v_cap , "cap" );

  AR |= HasCapCns;
  }

 Strong_Linking:
 // construct "strong linking" Constraint, that are formulation-agnostic- - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( wc & 4 ) )
  return;

 v_sfc.resize( f_n_facilities );
 add_dynamic_constraint( v_sfc , "strong" );

 AR |= HasStrngCns;
 
 /*!!
 #if CHECK_DS
  CheckAbsVSPhys();
 #endif
 !!*/

 }  // end( CapacitatedFacilityLocationBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_dynamic_constraints(
						       Configuration * dycc )
{
 if( ! ( AR & HasStrngCns ) )
  return;

 // read separation parameters
 double eps = 1e-4;
 int max = -1;
 if( ( ! dycc ) && f_BlockConfig )
  dycc = f_BlockConfig->f_dynamic_constraints_Configuration;
 if( auto sc =
  dynamic_cast< SimpleConfiguration< std::pair< int , double > > * >( dycc ) ) {
  max = sc->f_value.first;
  if( max == 0 )  // weird
   return;        // what else?
  eps = sc->f_value.second;
  }

 // prepare data structures
 using Violij = std::tuple< double , Index , Index >;
 std::vector< Violij > found;

 // get fractional solution
 CntSolution y( f_n_facilities );
 get_facility_solution( y.begin() );

 CntSolution x( f_n_facilities * f_n_customers );
 get_transportation_solution( x.begin() );

 // perform separation loop
 auto yit = y.begin();
 auto xit = x.begin();
 for( Index i = 0 ; i < f_n_facilities ; ++i , ++yit )
  for( Index j = 0 ; j < f_n_customers ; ++j , ++xit )
   if( auto viol = *xit - *yit ; viol >= eps )
    found.push_back( std::make_tuple( viol , i , j ) );

 // found is naturally ordered by increasing i, but if it's too long it has
 // to be reduced, which implies sorting on viol first and then re-sorting
 // on i to restore the original property
 if( ( max > 0 ) && ( int( found.size() ) > max ) ) {
  // sort by decreasing violation
  std::sort( found.begin() , found.end() ,
	     []( auto & a , auto & b ) {
	      return( std::get< 0 >( a ) > std::get< 0 >( b ) );
	      } );

  found.resize( max );  // eliminate unwanted part

  // re-sort by increasing i
  std::sort( found.begin() , found.end() ,
	     []( auto & a , auto & b ) {
	      return( std::get< 1 >( a ) < std::get< 1 >( b ) );
	      } );
  }

 // finally the actual insertion of the Constraint
 auto fit = found.begin();
 for( Index i = 0 ; i < f_n_facilities ; ++i ) {  // for each i
  // construct the list of [FRow]Constraint corresponding to y_i
  std::list< FRowConstraint > lst;
  auto * yi = get_y( i );
  while( ( fit != found.end() ) && ( std::get< 1 >( *fit ) == i ) ) {
   std::list< FRowConstraint > li( 1 );
   li.back().set_rhs( 0 );
   li.back().set_lhs( -Inf< RowConstraint::RHSValue >() );

   Index j = std::get< 2 >( *(fit++) );
   v_coeff_pair p( 2 );

   // the constraint is x_{ij} - y_i <= 0, but in the Flow Formulation
   // x_{ij} is scaled by the demand of j, so we implement it as
   // x_{ij} - v_demand[ j ] y_i <= 0
   p[ 0 ] = std::make_pair( yi , ( AR & FormMsk ) == FlwForm
			         ? - v_demand[ j ] : -1 );
   p[ 1 ] = std::make_pair( get_x( i , j ) , 1 );

   /*!!
   std::cerr << std::endl << "i = " << i << ", j = " << j << ", y_i = "
	     << yi->get_value() << ", x_{ij} = "
	     << p[ 1 ].first->get_value();
	     !!*/
 
   li.back().set_function( new LinearFunction( std::move( p ) ) );
   lst.splice( lst.end() , li );
   }

  // if any constraint concerning y_i was separated, add them all in one blow
  // note the "eNoBlck", which makes sense because this is an "abstract
  // representation only matter" and therefore there is no point in having
  // the corresponding BlockModAdd< FRowConstraint > to be scanned inside
  // CapacitatedFacilityLocationBlock::add_Modification()
  if( ! lst.empty() )
   add_dynamic_constraints( v_sfc[ i ] , lst , eNoBlck );

  if( fit == found.end() )  // if we have depleted the list
   break;                   // nothing else to do
  }
 }  // end( CapacitatedFacilityLocationBlock::generate_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_objective(
						       Configuration * objc )
{
 const static std::string _prfx =
                    "CapacitatedFacilityLocationBlock::generate_objective: ";
 
 if( AR & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

 if( ! ( AR & HasVar ) )
  throw( std::logic_error( _prfx + "generate_abstract_variables not called" ) );
 
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

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

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
  throw( std::logic_error( "CapacitatedFacilityLocationBlock::is_feasible:"
			   " generate_abstract_variables not called" ) );

 // check variable feasibility
 for( Index i = 0 ; i < f_n_facilities ; ++i )
  if( ! get_y( i )->is_feasible( eps ) )
   return( false );

 for( Index i = 0 ; i < f_n_facilities ; ++i )
  for( Index j = 0 ; j < f_n_customers ; ++j )
   if( ! get_x( i , j )->is_feasible( eps ) )
    return( false );

 // now check constraints feasibility
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
   for( auto & cnst : v_sat ) {
    cnst.compute();
    if( cnst.rel_viol() > eps )
     return( false );
    }

   return( true );
   }

  SimpleConfiguration< double > cfg( eps ); 
  return( MCFB( v_Block[ 1 ] )->is_feasible( true , & cfg ) );
  
  throw( std::logic_error( _prfx + "flow formulation not supported yet" ) );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   double tot = 0;
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    tot += get_x( i , j )->get_value();

   if( std::abs( 1 - tot ) > eps )
    return( false );
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::customer_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::facility_feasible( double eps ,
							  bool useabstract )
{
 const static std::string _prfx =
                     "CapacitatedFacilityLocationBlock::capacity_feasible: ";

 if( ! ( AR & HasCapCns ) )  // if capacity constraints are not defined
  useabstract = false;       // you cannot use them to check feasibility

 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == StdForm ) ||
      ( ( AR & FormMsk ) == FlwForm ) ) {
   for( auto & cnst : v_cap ) {
    cnst.compute();
    if( cnst.rel_viol() > eps )
     return( false );
    }

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
   auto tot = v_capacity[ i ] * get_y( i )->get_value();
   for( Index j = 0 ; j < f_n_customers ; ++j )
    tot -= v_demand[ j ] * get_x( i , j )->get_value();

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
 int wR3B = 0;
 if( auto tcfg = dynamic_cast< SimpleConfiguration< int > * >( r3bc ) )
  wR3B = tcfg->f_value;

 if( ( wR3B < 0 ) || ( wR3B > 2 ) )
  throw( std::invalid_argument(  _prfx + "invalid R3B type" ) );

 if( ! wR3B ) {  // "copy" R3B- - - - - - - - - - - - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // create the CapacitatedFacilityLocationBlock or take it from base
  CapacitatedFacilityLocationBlock * CFLB;
  if( base ) {
   CFLB = dynamic_cast< CapacitatedFacilityLocationBlock * >( base );
   if( ! CFLB )
    throw( std::invalid_argument(
     _prfx + "base is not a CapacitatedFacilityLocationBlock" ) );
   }
  else
   CFLB = new CapacitatedFacilityLocationBlock( father );

  // load the data of the problem
  CFLB->load( f_n_facilities , f_n_customers , v_capacity , v_f_cost ,
	      v_demand , v_t_cost , f_unsplittable );

  // now copy over the fixed-to-0 status, if any
  if( auto nf = std::count_if( v_fxd.begin() , v_fxd.end() ,
			       []( auto el ) { return( el == yFxd0 ); } ) ) {
   Subset tfx( nf );
   auto tfxit = tfx.begin();
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    if( v_fxd[ i ] == yFxd0 )
     *(tfxit++) = i;
   CFLB->close_facilities( std::move( tfx ) , true , eNoMod , eNoMod );
   }

  // now copy over the fixed-to-1 status, if any
  if( auto nf = std::count_if( v_fxd.begin() , v_fxd.end() ,
			       []( auto el ) { return( el == yFxd1 ); } ) ) {
   Subset tfx( nf );
   auto tfxit = tfx.begin();
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    if( v_fxd[ i ] == yFxd1 )
     *(tfxit++) = i;
   CFLB->fix_open_facilities( std::move( tfx ) , true , eNoMod , eNoMod );
   }

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

 guts_of_get_R3B_MCF( MCFB , wR3B );

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

 if( ! ( ws & FormMsk ) )  // actually nothing to map back
  return;                  // silently (and cowardly) return

 int wR3B = 0;
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

  if( ws & 1 ) {
   CntSolution y( f_n_facilities );
   CFLB->get_facility_solution( y.begin() );
   set_facility_solution( y.begin() );
   }

  if( ws & 2 ) {
   CntSolution x( f_n_facilities * f_n_customers );
   CFLB->get_transportation_solution( x.begin() );
   set_transportation_solution( x.begin() );
   }

  return;
  }

 // "flow relaxation" R3B - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast< MCFBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( _prfx + "R3B is not a MCFBlock" ) );

 if( ( MCFB->get_NNodes() != f_n_facilities + f_n_customers + 1 ) ||
     ( MCFB->get_NArcs() != f_n_facilities * ( f_n_customers + 1 ) +
                            ( wR3B == 2 ? f_n_customers : 0  ) ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

 Index l = ws & 1 ? 0 : f_n_facilities;
 Index u = ws & 2 ? f_n_facilities * ( f_n_customers + 1 ) : f_n_facilities;
 
 MCFBlock::Vec_FNumber x( u - l );
 auto it = x.begin();

 MCFB->get_x( it , Range( l , u ) );

 if( ws & 1 ) {
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   *(it++) /= v_capacity[ i ];
  set_facility_solution( x.begin() );
  }

 if( ws & 2 ) {
  auto cit = it;
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j )
    *(it++) /= v_demand[ j ];
  set_transportation_solution( cit );
  }
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

 if( ! ( ws & FormMsk ) )  // actually nothing to map forward
  return;                  // silently (and cowardly) return

 int wR3B = 0;
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
                            ( wR3B == 2 ? f_n_customers : 0 ) ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3B" ) );

 Index l = ws & 1 ? f_n_facilities : 0;
 Index u = ws & 2 ? f_n_facilities * ( f_n_customers + 1 ) : f_n_facilities;
 
 MCFBlock::Vec_FNumber x( u - l );

 auto it = x.begin();
 if( ws & 1 ) {
  get_facility_solution( it );
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   *(it++) *= v_capacity[ i ];
  }

 if( ws & 2 ) {
  get_transportation_solution( it );
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   for( Index j = 0 ; j < f_n_customers ; ++j )
    *(it++) *= v_demand[ j ];
  }

 MCFB->set_x( x.begin() , Range( l , u ) );

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
 int wR3B = 0;
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
                            ( wR3B == 2 ? f_n_customers : 0 ) ) )
  throw( std::invalid_argument( _prfx + "incompatible sizes in R3Block" ) );

 // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is the "nuclear option": the CapacitatedFacilityLocationBlock has
 // been re-loaded, so the MCFBlock must be re-loaded. treat this here since
 // it is the only case where wR3B is needed

 if( auto tmod = dynamic_cast< const NBModification * >( mod ) ) {
  guts_of_get_R3B_MCF( MCFB , wR3B );
  return( true );
  }

 return( guts_of_map_f_Mod_MCF( MCFB , mod , issuePMod , issueAMod ) );

 }  // end( CapacitatedFacilityLocationBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::map_back_Modification( Block * R3B ,
			            c_p_Mod mod , Configuration *r3bc ,
				    ModParam issuePMod , ModParam issueAMod )
{
 const static std::string _prfx =
                 "CapacitatedFacilityLocationBlock::map_back_Modification: ";
 int wR3B = 0;
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
                            ( wR3B == 2 ? f_n_customers : 0 ) ) )
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

void CapacitatedFacilityLocationBlock::get_facility_solution( IS_it Sol ,
							      Range rng )
 const { get_y< bool >( Sol , rng ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::get_facility_solution( IS_it Sol ,
							      c_Subset & nms )
 const { get_y< bool >( Sol , nms ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::get_facility_solution( CS_it Sol ,
							      Range rng )
 const { get_y< double >( Sol , rng ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::get_facility_solution( CS_it Sol ,
							      c_Subset & nms )
 const { get_y< double >( Sol , nms ); }

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_transportation_solution(
					        CS_it Sol , Range rng ) const
{
 get_x< double >( Sol , rng );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::get_transportation_solution(
					    CS_it Sol , c_Subset & nms ) const
{
 get_x< double >( Sol , nms );
 }

/*--------------------------------------------------------------------------*/

RealObjective::OFValue
                CapacitatedFacilityLocationBlock::get_objective_value( void )
{
 if( ! ( AR & HasObj ) )  // the objective is not there
  return( Inf< RealObjective::OFValue >() );

 if( ( AR & FormMsk ) == StdForm ) {  // "natural formulation" (NF) - - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  f_obj.compute();
  return( f_obj.value() );
  }

 if( ( AR & FormMsk ) == KskForm ) {  // "knapsack formulation" (KF)- - - - -
                                      //- - - - - - - - - - - - - - - - - - -
  RealObjective::OFValue res = 0;
  for( auto ki : v_Block )
   res += BKB( ki )->get_objective_value();

  return( res );
  }

 // else it is the "flow formulation" (FF)- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the objective is split between the two sub-Block

 f_obj.compute();
 return( f_obj.value() + MCFB( v_Block[ 1 ] )->get_objective_value() );

 }  // end( CapacitatedFacilityLocationBlock::get_objective_value )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::set_facility_solution( c_IS_it Sol ,
							      Range rng )
{
 set_y< bool >( Sol , rng );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::set_facility_solution( c_IS_it Sol ,
							      c_Subset & nms )
{
 set_y< bool >( Sol , nms );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::set_facility_solution( c_CS_it Sol ,
							      Range rng )
{
 set_y< double >( Sol , rng );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::set_facility_solution( c_CS_it Sol ,
							      c_Subset & nms )
{
 set_y< double >( Sol , nms );
 }

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::set_transportation_solution(
					            c_CS_it Sol , Range rng )
{
 set_x< double >( Sol , rng );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void CapacitatedFacilityLocationBlock::set_transportation_solution(
					       c_CS_it Sol , c_Subset & nms )
{
 set_x< double >( Sol , nms );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::add_Modification( sp_Mod mod ,
							 ChnlName chnl )
{
 /* This add_Modification() is "nonstandard" in the sense that
  *
  *     IT ALSO USES "PHYSICAL Modification" TO UPDATE THE ABSTRACT
  *     REPRESENTATION
  *
  * The issue is that in the KF and FF part of the abstract representation of
  * the CapacitatedFacilityLocationBlock actually "lives" inside the inner
  * BinaryKnapsackBlock or MCFBlock. Indeed, said BinaryKnapsackBlock or
  * MCFBlock are themselves the "abstract representation" of the
  * CapacitatedFacilityLocationBlock rather than the "physical" one.
  *
  * Anyhow, any change in the abstract representation of the
  * BinaryKnapsackBlock or MCFBlock is "intercepted" by these :Block, which
  * update their physical representation
  *
  *     AND RESET THE concerns_Block() BIT
  *
  * before passing them to the father CapacitatedFacilityLocationBlock. This
  * means that the "standard trick" of relying upon the concerns_Block() bit
  * to identify the Modification that need (not) be processed by the
  * CapacitatedFacilityLocationBlock does not work. The solution to this
  * issue is that CapacitatedFacilityLocationBlock ensures that f_mod_skip
  * == true if and only if this is what is happening. That is, f_mod_skip is
  * false by default, it is set to true right before calling the chg_*
  * methods of the sub-Block (within which "physical" Modification are issued
  * and CapacitatedFacilityLocationBlock::add_Modification() is called) and
  * put back to false immediately after this happened. This would be dangerous
  * in case multiple changes would be happening at the same time, with some of
  * them "known already" by CapacitatedFacilityLocationBlock and others not,
  * but this is not supposed to happen since the mechanism is only activated
  * inside the chg_* methods of CapacitatedFacilityLocationBlock, which are
  * only supposed to be called when it is lock()-ed.
  *
  * Furthermore, this "nonstandard" mechanism implies that
  *
  *     CapacitatedFacilityLocationBlock CAN USE THE "PHYSICAL
  *     Modification" OF BinaryKnapsackBlock AND MCFBlock TO
  *     PERFORM THE UPDATE OF ITS PHYSICAL REPRESENTATION
  *
  * These are in fact less than the "abstract Modification" produced by
  * changing the abstract representation and "more informative", so the
  * code is actually quite a bit simpler because of this. */
 //!! std::cout << *mod << std::endl;

 if( f_mod_skip )           // this is a "physical-abstract" Modification
  goto do_the_usual_stuff;  // just do the usual stuff

 // otherwise process it for abstract representation -> physical
 // representation synchronization: this must be done on a case-by-case basis
 switch( AR & FormMsk ) {
  case( StdForm ):  // Standard Formulation
   if( mod->concerns_Block() ) {  // the usual drill
    mod->concerns_Block( false );
    guts_of_add_ModificationSF( mod.get() , chnl );
    }
   break;

  case( KskForm ):  // Knapsack Formulation
   if( mod->concerns_Block() )  // if it's abstract
    // it's an error: nothing can be changed in the "root" Block
    throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

   // it can't be a physical Modification coming from the "root" Block
   // either since these won't pass from here, hence it must be coming
   // from some BinaryKnapsackBlockMod
   if( auto bmod = dynamic_cast< BinaryKnapsackBlockMod * >( mod.get() ) )
    guts_of_add_ModificationKFP( bmod , chnl );
   else
    if( auto gmod = dynamic_cast< GroupModification * >( mod.get() ) )
     guts_of_add_ModificationKFG( gmod , chnl );

   break;

  default:          // Flow Formulation
   if( mod->concerns_Block() ) {  // if it's abstract, the usual drill
    mod->concerns_Block( false );
    guts_of_add_ModificationFFA( mod.get() , chnl );
    }
   else
    if( auto mmod = dynamic_cast< MCFBlockMod * >( mod.get() ) )
     guts_of_add_ModificationFFP( mmod , chnl );
    else
     if( auto gmod = dynamic_cast< GroupModification * >( mod.get() ) )
      guts_of_add_ModificationFFG( gmod , chnl );

  }  // end( switch )

 // finally, do the usual stuff (pass to the Solver and to the father)
 do_the_usual_stuff:

 Block::add_Modification( mod , chnl );

 }  // end( CapacitatedFacilityLocationBlock::add_Modification )

/*--------------------------------------------------------------------------*/
/*-------- PRINTING & SAVING THE CapacitatedFacilityLocationBlock ----------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::print( std::ostream & output ,
					      char vlvl ) const
{
 if( std::tolower( vlvl ) != 'c' ) {  // non-complete version
  // only basic information 
  output << "CapacitatedFacilityLocationBlock with: " << f_n_facilities
	 << " facilities and " << f_n_customers << " customers" << std::endl;
  }
 else  {
  // complete version: file in standard ORLib format
  output << f_n_facilities << std::endl;
  output << f_n_customers << std::endl << std::endl;

  output << std::setprecision( 16 );
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   output << v_capacity[ i ] << "\t" << v_f_cost[ i ] << std::endl;
  
  output << std::endl;

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   output << v_demand[ j ] << std::endl;
   for( Index i = 0 ; i < f_n_facilities ; ++i )
    output << v_t_cost[ i ][ j ] << "\t";
   output << std::endl;
   }
  }
 }  // end( CapacitatedFacilityLocationBlock::print )

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
              { nf , nc } , v_t_cost );

 if( std::any_of( v_fxd.begin() , v_fxd.end() ,
		  []( auto fi ) { return( fi != yFree ); } ) )
  ( group.addVar( "FacilityFix" , netCDF::NcChar() , nf )
    ).putVar( v_fxd.data() );

 if( f_unsplittable )
  group.addDim( "UnSplittable" , 1 );

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

 c_Index num = rng.second - rng.first;
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
				   un_ModBlock( issueAMod ) );
  else {
   f_mod_skip = true;
   for( Index i = rng.first ; i < rng.second ; ++i )
    BKB( v_Block[ i ] )->chg_profit( *(NCost++) , f_n_customers ,
				     issueMod , issueAMod );
   f_mod_skip = false;
   }
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
			    CapacitatedFacilityLocationBlockMod::eChgFCost ,
			    rng ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( CapacitatedFacilityLocationBlock::chg_facility_costs( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_facility_costs(
		               c_CV_it NCost , Subset && nms , bool ordered ,
			       ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( nms.back() >= f_n_facilities )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_f_cost.data() , nms , NCost , f_n_facilities ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_f_cost.data() , nms , NCost );

  if( ( AR & FormMsk ) != KskForm )
   // since modify_coefficients owns both vectors, two copies are made
   get_lfo()->modify_coefficients( CVector( NCost , NCost + nms.size() ) ,
				   Subset( nms ) , ordered ,
				   un_ModBlock( issueAMod ) );
  else {
   f_mod_skip = true;
   for( auto i : nms )
    BKB( v_Block[ i ] )->chg_profit( *(NCost++) , f_n_customers ,
				     issueMod , un_ModBlock( issueAMod ) );
   f_mod_skip = false;
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_f_cost.data() , nms , NCost );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;

 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
                            CapacitatedFacilityLocationBlockMod::eChgFCost ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }
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
   get_lfo()->modify_coefficient( i , NCost , un_ModBlock( issueAMod ) );
  else {
   f_mod_skip = true;
   BKB( v_Block[ i ] )->chg_profit( NCost , f_n_customers ,
				    issueMod , un_ModBlock( issueAMod ) );
   f_mod_skip = false;
   }
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

 }  // end( CapacitatedFacilityLocationBlock::chg_facility_cost )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_costs(
				     c_CV_it NCost , Range rng ,
				     ModParam issueMod , ModParam issueAMod )
{
 c_Index maxn = f_n_facilities * f_n_customers;
 rng.second = std::min( rng.second , maxn );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 c_Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NCost , NCost + num , v_t_cost.data() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCost , NCost + num , v_t_cost.data() + rng.first );

  f_mod_skip = true;
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // since modify_coefficients owns the vector, a copy has to be made
    CVector NC( NCost , NCost + num );
    get_lfo()->modify_coefficients( std::move( NC ) ,
				    Range( rng.first + f_n_facilities ,
					   rng.second + f_n_facilities ) ,
				    un_ModBlock( issueAMod ) );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    Index f = rng.first;
    Index i = f / f_n_customers;
    Index l = f % f_n_customers;
    if( ( ( rng.second - 1 ) / f_n_customers ) == i ) {
     // the range is all inside a single facility
     BKB( v_Block[ i ] )->chg_profits( NCost , Range( l , l + num ) ,
				       issueMod , issueAMod );
     break;
     }

    // open a new channel to bunch up all  Modification
    not_ModBlock( issueAMod );
    auto iAM = open_if_needed( issueAMod , 2 );

    // the range of the first facility does not necessarily start from 0,
    // but it surely ends at f_n_customers
    BKB( v_Block[ i++ ] )->chg_profits( NCost , Range( l , f_n_customers ) ,
					issueMod , iAM );
    NCost += ( f_n_customers - l );
    f += ( f_n_customers - l );
 
    // the range of all other facilities starts from 0, but it does not
    // necessarily end at f_n_customers
    for( ; ; ++i , NCost += f_n_customers ) {
     Index nf = f + f_n_customers;
     if( nf >= rng.second ) {  // last facility
      BKB( v_Block[ i ] )->chg_profits( NCost , Range( 0 , rng.second - f ) ,
					issueMod , iAM );
      break;
      }
     else {
      BKB( v_Block[ i ] )->chg_profits( NCost , Range( 0 , f_n_customers ) ,
					issueMod , iAM );
      f = nf;
      }
     }

    close_if_needed( iAM , 2 );  // close the new channel
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_tcost_MCF( MCFB( v_Block[ 1 ] ) , rng ,
			   issueMod , issueAMod );
   }

  f_mod_skip = false;
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

 }  // end CapacitatedFacilityLocationBlock::chg_transportation_costs( range )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_costs(
			       c_CV_it NCost , Subset && nms , bool ordered ,
			       ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 // ensure the names are ordered even if they were not so originally

 c_Index maxn = f_n_facilities * f_n_customers;
 if( nms.back() >= maxn )
  throw( std::invalid_argument( "invalid name of ( facility , customer ) pair"
				) );

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_t_cost.data() , nms , NCost , maxn ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_t_cost.data() , nms , NCost );

  f_mod_skip = true;
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    Subset nnms( nms );     // copy and translate names
    for( auto & el : nnms )
     el += f_n_facilities;

    // since modify_coefficients owns both vectors, copies has to be made
    CVector NC( NCost , NCost + nms.size() );
    get_lfo()->modify_coefficients( std::move( NC ) , std::move( nnms ) ,
				    ordered , un_ModBlock( issueAMod ) );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    // this operation is horribly complex if nms is not ordered, so order it;
    // but this means also re-ordering the values accordingly
    auto tNC = NCost;
    CVector oNC;
    if( ! ordered ) {
     using ICPair = std::pair< Index , Cost >;
     std::vector< ICPair > tmp( nms.size() );
     for( Index i = 0 ; i < nms.size() ; ++i )
      tmp[ i ] = std::make_pair( nms[ i ] , *(NCost++) );
     std::sort( tmp.begin() , tmp.end() ,
		[]( auto & a , auto & b ) { return( a.first < b.first ); } );
     oNC.resize( nms.size() );
     for( Index i = 0 ; i < nms.size() ; ++i ) {
      nms[ i ] = tmp[ i ].first;
      oNC[ i ] = tmp[ i ].second;
      }
     tNC = oNC.begin();
     ordered = true;
     }

    Index i = nms.front() / f_n_customers;
    if( ( nms.back() / f_n_customers ) == i ) {
     // the subset is all inside a single facility
     Subset nnms( nms );     // copy and translate names
     for( auto & el : nnms )
      el %= f_n_customers;
     BKB( v_Block[ i ] )->chg_profits( tNC , std::move( nnms ) , true ,
				       issueMod , issueAMod );
     break;
     }

    // open a new channel to bunch up all Modification
    not_ModBlock( issueAMod );
    auto iAM = open_if_needed( issueAMod , 2 );

    for( auto bit = nms.begin() , eit = bit ; ; ++i ) {
     while( ( eit != nms.end() ) && ( *eit / f_n_customers ) == i )
      ++eit;

     if( eit == bit )  // the subset is empty for this i
      continue;        // move to next i

     Subset nnms( bit , eit );  // copy and translate names
     for( auto & el : nnms )
      el %= f_n_customers;

     BKB( v_Block[ i ] )->chg_profits( tNC , std::move( nnms ) , true ,
				       issueMod , iAM );
     if( eit == nms.end() )
      break;

     tNC += std::distance( bit , eit );
     bit = eit;
     }

    close_if_needed( iAM , 2 );  // close the new channel
    break;
    }
   default:    // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_tcost_MCF( MCFB( v_Block[ 1 ] ) , nms , ordered ,
			   issueMod , issueAMod );
   }

  f_mod_skip = false;
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_t_cost.data() , nms , NCost );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgTCost ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }
 }  // end CapacitatedFacilityLocationBlock::chg_transportation_costs( range )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_transportation_cost( Cost NCost ,
			   Index p , ModParam issueMod , ModParam issueAMod )
{
 c_Index maxn = f_n_facilities * f_n_customers;
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

  f_mod_skip = true;
  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    get_lfo()->modify_coefficient( p + f_n_facilities , NCost ,
				   un_ModBlock( issueAMod ) );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    BKB( v_Block[ i ] )->chg_profit( NCost , j , issueMod , issueAMod );
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_tcost_MCF( MCFB( v_Block[ 1 ] ) , Range( p , p + 1 ) ,
			   issueMod , issueAMod );
   }

  f_mod_skip = false;
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

 Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NCap , NCap + num , v_capacity.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasCapCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NCap , NCap + num , v_capacity.begin() + rng.first );

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  if( ( AR & FormMsk ) == FlwForm )
   ++num;
  auto iAM = open_if_needed( issueAMod , num );
  f_mod_skip = true;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( Index i = rng.first ; i < rng.second ; ++i )
     LF( v_cap[ i ].get_function()
	 )->modify_coefficient( f_n_customers , - *(NCap++) , iAM );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( Index i = rng.first ; i < rng.second ; ++i )
     BKB( v_Block[ i ]
	  )->chg_weight( - *(NCap++) , f_n_customers , issueMod , iAM );
    break;
    }
   default: {  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_ucaps( NCap , rng , issueMod , iAM );
    for( Index i = rng.first ; i < rng.second ; ++i )
     LF( v_cap[ i ].get_function()
	 )->modify_coefficient( 1 , - *(NCap++) , iAM );
    }
   }

  f_mod_skip = false;

  // if a new channel had been opened, close it
  close_if_needed( iAM , num );
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
			        c_DV_it NCap , Subset && nms , bool ordered ,
				ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( nms.back() >= f_n_facilities )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_capacity.data() , nms , NCap , f_n_facilities ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasCapCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_capacity.data() , nms , NCap );

  // if appropriate, open a new channel to bunch up all abstract Modification
  Index num = nms.size() + ( ( AR & FormMsk ) == FlwForm ? 1 : 0 );
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , num );
  f_mod_skip = true;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto i : nms )
     LF( v_cap[ i ].get_function()
	 )->modify_coefficient( f_n_customers , - *(NCap++) , iAM );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto i : nms )
     BKB( v_Block[ i ]
	  )->chg_weight( - *(NCap++) , f_n_customers , issueMod , iAM );
    break;
    }
   default: {  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    MCFB( v_Block[ 1 ] )->chg_ucaps( NCap , Subset( nms ) , ordered ,
				     issueMod , iAM );
    for( auto i : nms )
     LF( v_cap[ i ].get_function()
	 )->modify_coefficient( 1 , - *(NCap++) , iAM ); 
    }
   }

  f_mod_skip = false;
  // if a new channel had been opened, close it
  close_if_needed( iAM , num );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_capacity.data() , nms , NCap );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
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
  f_mod_skip = true;
  not_ModBlock( issueAMod );

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    LF( v_cap[ i ].get_function()
	)->modify_coefficient( f_n_customers , - NCap , issueAMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    BKB( v_Block[ i ]
	 )->chg_weight( - NCap , f_n_customers , issueMod , issueAMod );
    break;
    }
   default: {  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    auto iAM = open_if_needed( issueAMod , 2 );
    MCFB( v_Block[ 1 ] )->chg_ucap( NCap , i , issueMod , iAM );
    LF( v_cap[ i ].get_function() )->modify_coefficient( 1 , - NCap , iAM );
    close_if_needed( iAM , 2 );
    }
   }

  f_mod_skip = false;
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

void CapacitatedFacilityLocationBlock::chg_customer_demands( c_DV_it NDem ,
			 Range rng , ModParam issueMod , ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_n_customers );
 if( rng.second <= rng.first )  // nothing to change
  return;                       // cowardly (and silently) return

 c_Index num = rng.second - rng.first;
 // TODO: if some changes are "fake", rather restrict the range
 if( std::equal( NDem , NDem + num , v_demand.begin() + rng.first ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasSatCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  std::copy( NDem , NDem + num , v_demand.begin() + rng.first );

  // if appropriate, open a new channel to bunch up all abstract Modification
  Index nc = ( ( AR & FormMsk ) == FlwForm ) ? 0 : f_n_facilities;
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , nc );
  f_mod_skip = true;

  switch( AR & FormMsk ) {
   case( StdForm ):    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto & capi : v_cap )
     LF( capi.get_function()
	 )->modify_coefficients( DVector( NDem , NDem + num ) , rng , iAM );
    break;

   case( KskForm ):    // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto bi : v_Block )
     BKB( bi )->chg_weights( NDem , rng , issueMod , iAM );

    break;

   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_dem_MCF( MCFB( v_Block[ 1 ] ) , rng , issueMod , issueAMod );
   }

  f_mod_skip = false;
  // if a new channel had been opened, close it
  close_if_needed( iAM , nc );
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

 }  // end( CapacitatedFacilityLocationBlock::chg_customer_demands( rng ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_customer_demands( c_DV_it NDem ,
			            Subset && nms , bool ordered ,
				    ModParam issueMod , ModParam issueAMod )
{
 if( nms.empty() )  // nothing to change
  return;           // cowardly (and silently) return

 if( nms.back() >= f_n_customers )
  throw( std::invalid_argument( "invalid customer name" ) );

 // TODO: eliminate from nms the "fake" changes
 if( is_equal( v_demand.data() , nms , NDem , f_n_customers ) )
  return;  // actually nothing changes, avoid issuing the Modification

 if( not_dry_run( issueAMod ) && ( AR & HasSatCns ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  copyidx( v_demand.data() , nms , NDem );

  // if appropriate, open a new channel to bunch up all abstract Modification
  Index nc = ( ( AR & FormMsk ) == FlwForm ) ? 0 : f_n_facilities;
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , nc );
  f_mod_skip = true;

  switch( AR & FormMsk ) {
   case( StdForm ):    // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto & capi : v_cap )
     LF( capi.get_function()
	 )->modify_coefficients( DVector( NDem , NDem + nms.size() ) ,
				 Subset( nms ) , ordered , iAM );
    break;
 
   case( KskForm ):    // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto bi : v_Block )
     BKB( bi )->chg_weights( NDem , Subset( nms ) , ordered ,
			     issueMod , iAM );

    break;

   default:    // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_dem_MCF( MCFB( v_Block[ 1 ] ) , nms , ordered ,
			 issueMod , issueAMod );
   }

  f_mod_skip = false;
  // if a new channel had been opened, close it
  close_if_needed( iAM , nc );
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   copyidx( v_demand.data() , nms , NDem );

 f_cond_lower = NAN;  // reset conditional bounds
 f_cond_upper = NAN;
 
 if( issue_pmod( issueMod ) ) {  // issue "physical Modification" - - - - - -
  if( ! ordered )
   std::sort( nms.begin() , nms.end() );
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
			    CapacitatedFacilityLocationBlockMod::eChgDem ,
			    std::move( nms ) ) ,
			   Observer::par2chnl( issueMod ) );
  }
 }  // end( CapacitatedFacilityLocationBlock::chg_customer_demands( sbst ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_customer_demand( Demand NDem ,
			   Index j , ModParam issueMod , ModParam issueAMod )
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

  // if appropriate, open a new channel to bunch up all abstract Modification
  Index nc = ( ( AR & FormMsk ) == FlwForm ) ? 0 : f_n_facilities;
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , nc );
  f_mod_skip = true;

  switch( AR & FormMsk ) {
   case( StdForm ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto & capi : v_cap )
     LF( capi.get_function() )->modify_coefficient( j , NDem , iAM );

    break;
 
   case( KskForm ):    // - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( auto bi : v_Block )
     BKB( bi )->chg_weight( NDem , j , issueMod , iAM );

    break;

   default:    // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - -
    guts_of_chg_dem_MCF( MCFB( v_Block[ 1 ] ) , Range( j , j + 1 ) ,
			 issueMod , issueAMod );
   }

  f_mod_skip = false;
  // if a new channel had been opened, close it
  close_if_needed( iAM , nc );
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
			    Range( j , j + 1 ) ) ,
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

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( v_fxd[ i ] != yFree )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm ) {
   for( Index i = rng.first ; i < rng.second ; ++i )
    if( v_fxd[ i ] == yFree ) {
     v_fxd[ i ] = yFxd0;
     v_y[ i ].set_value( 0 );
     v_y[ i ].is_fixed( true , iAM );
     }
   }
  else {
   f_mod_skip = true;
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd0;
    BKB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , iAM );
    }
   f_mod_skip = false;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueMod ) )
   for( Index i = rng.first ; i < rng.second ; ++i )
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd0;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( v_fxd[ i ] != yFree )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm ) {
   for( auto i : nms )
    if( v_fxd[ i ] == yFree ) {
     v_fxd[ i ] = yFxd0;
     v_y[ i ].set_value( 0 );
     v_y[ i ].is_fixed( true , iAM );
     }
   }
  else {
   f_mod_skip = true;
   for( auto i : nms ) {
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd0;
    BKB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , iAM );
    }
   f_mod_skip = false;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueMod ) )
   for( auto i : nms )
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd0;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 if( v_fxd[ i ] != yFree )  // fixed already
  return;                   // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {
  v_fxd[ i ] = yFxd0;

  if( ( AR & FormMsk ) != KskForm ) {
   v_y[ i ].set_value( 0 );
   v_y[ i ].is_fixed( true , un_ModBlock( issueAMod ) );
   }
  else {
   f_mod_skip = true;
   BKB( v_Block[ i ] )->fix_x( false , f_n_customers , issueMod , issueAMod );
   f_mod_skip = false;
   }
  }
 else
  if( not_dry_run( issueMod ) )
   v_fxd[ i ] = yFxd0;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( v_fxd[ i ] != yFree )
   ++cnt;

 if( ! cnt )  // all facilities are open already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm ) {
   for( Index i = rng.first ; i < rng.second ; ++i )
    if( v_fxd[ i ] != yFree ) {
     v_fxd[ i ] = yFree;
     v_y[ i ].is_fixed( false , iAM );
     }
   }
  else {
   f_mod_skip = true;
   for( Index i = rng.first ; i < rng.second ; ++i )
    if( v_fxd[ i ] != yFree ) {
     v_fxd[ i ] = yFree;
     BKB( v_Block[ i ] )->unfix_x( f_n_customers , issueMod , iAM );
     }
   f_mod_skip = false;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueMod ) )
   for( Index i = rng.first ; i < rng.second ; ++i )
    v_fxd[ i ] = yFree;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( v_fxd[ i ] != yFree )
   ++cnt;

 if( ! cnt )  // all facilities are open already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm ) {
   for( auto i : nms )
    if( v_fxd[ i ] != yFree ) {
     v_fxd[ i ] = yFree;
     v_y[ i ].is_fixed( false , iAM );
     }
   }
  else {
   f_mod_skip = true;
   for( auto i : nms )
    if( v_fxd[ i ] != yFree ) {
     v_fxd[ i ] = yFree;
     BKB( v_Block[ i ] )->unfix_x( f_n_customers , issueMod , iAM );
     }
   f_mod_skip = false;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueMod ) )
   for( auto i : nms )
    v_fxd[ i ] = yFree;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockSbstMod >( this ,
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

 if( v_fxd[ i ] == yFree )  // open already
  return;                   // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {
  v_fxd[ i ] = yFree;

  if( ( AR & FormMsk ) != KskForm )
   v_y[ i ].is_fixed( false , un_ModBlock( issueAMod ) );
  else {
   f_mod_skip = true;
   BKB( v_Block[ i ] )->unfix_x( f_n_customers , issueMod , issueAMod );
   f_mod_skip = false;
   }
  }
 else
  if( not_dry_run( issueMod ) )
   v_fxd[ i ] = yFree;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockRngdMod >( this ,
			     CapacitatedFacilityLocationBlockMod::eOpenF ,
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

 // TODO: if some changes are "fake", restrict the range
 Index cnt = 0;
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( v_fxd[ i ] == yFree )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm )
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    if( v_fxd[ i ] == yFree ) {
     v_fxd[ i ] = yFxd1;
     v_y[ i ].set_value( 1 );
     v_y[ i ].is_fixed( true , iAM );
     }
    }
  else {
   f_mod_skip = true;
   for( Index i = rng.first ; i < rng.second ; ++i ) {
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd1;
    BKB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }
   f_mod_skip = true;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueAMod ) )
   for( Index i = rng.first ; i < rng.second ; ++i )
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd1;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 // ensure the names are ordered even if they were not so originally
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );

 if( nms.back() >= get_NFacilities() )
  throw( std::invalid_argument( "invalid facility name" ) );

 // TODO: if some changes are "fake", restrict the subset
 Index cnt = 0;
 for( auto i : nms )
  if( v_fxd[ i ] == yFree )
   ++cnt;

 if( ! cnt )  // all facilities are fixed already
  return;     // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {

  // if appropriate, open a new channel to bunch up all abstract Modification
  not_ModBlock( issueAMod );
  auto iAM = open_if_needed( issueAMod , cnt );

  if( ( AR & FormMsk ) != KskForm ) {
   for( auto i : nms )
    if( v_fxd[ i ] == yFree ) {
     v_fxd[ i ] = yFxd1;
     v_y[ i ].set_value( 1 );
     v_y[ i ].is_fixed( true , iAM );
     }
   }
  else {
   f_mod_skip = true;
   for( auto i : nms ) {
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd1;
    BKB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , iAM );
    }
   f_mod_skip = false;
   }

  // if a new channel had been opened, close it
  close_if_needed( iAM , cnt );
  }
 else
  if( not_dry_run( issueMod ) )
   for( auto i : nms )
    if( v_fxd[ i ] == yFree )
     v_fxd[ i ] = yFxd1;
   

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

 if( v_fxd[ i ] != yFree )  // fixed already
  return;                   // nothing to do

 if( ( AR & HasVar ) && not_dry_run( issueAMod ) ) {
  v_fxd[ i ] = yFxd1;
 
  if( ( AR & FormMsk ) != KskForm ) {
   v_y[ i ].set_value( 1 );
   v_y[ i ].is_fixed( true , un_ModBlock( issueAMod ) );
   }
  else {
   f_mod_skip = true;
   BKB( v_Block[ i ] )->fix_x( true , f_n_customers , issueMod , issueAMod );
   f_mod_skip = false;
   }
  }
 else
  if( not_dry_run( issueMod ) )
   v_fxd[ i ] = yFxd1;

 // conditional bounds could be reset if they were computed looking at
 // facilities fixings, but they are not, and therefore they are not (reset)

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

void CapacitatedFacilityLocationBlock::chg_UnSplittable( bool unsplt ,
							 ModParam issueMod  ,
							 ModParam issueAMod )
{
 if( unsplt == f_unsplittable )  // changing to the same value
  return;                        // nothing to do

 // TODO: properly package the possibly very many individual Modification in
 //       some appropriate GroupModification
 
 if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
  // change abstract and physical representation together - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in the meantime, if so instructed also issue abstract Modification
  f_unsplittable = unsplt;

  switch( AR & FormMsk ) {
   case( StdForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    auto xit = v_x.data();
    for( const auto xend = xit + f_n_facilities * f_n_customers ; xit != xend ;
	 )
     (xit++)->is_integer( unsplt , issueMod );
    break;
    }
   case( KskForm ): {  // - - - - - - - - - - - - - - - - - - - - - - - - - -
    BinaryKnapsackBlock::boolVec vI( f_n_customers );
    std::fill( vI.begin() , vI.end() , unsplt );
    f_mod_skip = true;
    for( auto bi : v_Block )
     BKB( bi )->chg_integrality( vI.begin() , Range( 0 , f_n_customers ) ,
				 issueMod , issueAMod );
    f_mod_skip = false;
    break;
    }
   default:  // FlwForm - - - - - - - - - - - - - - - - - - - - - - - - - - -
    throw( std::invalid_argument(
	   "unsplittable problem not supported with the Flow Formulation" ) );
   }
  }
 else
  // only change the physical representation- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( not_dry_run( issueMod ) )
   f_unsplittable = unsplt;

 if( issue_pmod( issueMod ) )  // issue "physical Modification" - - - - - - -
  Block::add_Modification( std::make_shared<
			   CapacitatedFacilityLocationBlockMod >( this ,
								  unsplt
			   ? CapacitatedFacilityLocationBlockMod::eChgUnSplt
			   : CapacitatedFacilityLocationBlockMod::eChgSplt ) ,
			   Observer::par2chnl( issueMod ) );

 }  // end( CapacitatedFacilityLocationBlock::chg_UnSplittable )

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
  * afterward which means that any listening Observer already knows that
  * none of the previous Variable and Constraint are valid any longer. */

 for( auto & lst : v_sfc )  // clear the strong forcing constraints
  for( auto & cnst : lst )
   cnst.clear();
 for( auto & cnst : v_cap )  // clear the capacity constraints
  cnst.clear(); 
 for( auto & cnst : v_sat )  // clear the satisfaction constraints
  cnst.clear();
 f_obj.clear();              // clear the objective function

 if( ( AR & FormMsk ) == FlwForm ) {
  // AbstractBlock assumes it is the sole owner of its Constraint and
  // Variable, but this is not true here, so remove them before
  // deleting it to avoid double deletion
  AB( v_Block[ 0 ] )->reset_static_constraints();
  AB( v_Block[ 0 ] )->reset_static_variables();
  AB( v_Block[ 0 ] )->reset_objective();
  // clear the capacity constraints before the AbstractBlock is destroyed
  // to avoid the destruction looking at pointers to the deleted Block
  v_cap.clear();
  // detach the objective from the AbstractBlock, since the latter will
  // be deleted before the former
  f_obj.set_Block( nullptr );
  }

 // delete all sub-Block
 for( auto bi : v_Block )
  delete bi;

 v_Block.clear();  // then clear the vector

 // then delete them all
 v_sfc.clear();
 v_cap.clear();
 v_sat.clear();

 // delete all Variable
 v_x.resize( boost::extents[ 0 ][ 0 ] );
 v_y.clear();

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

void CapacitatedFacilityLocationBlock::guts_of_get_R3B_MCF( MCFBlock * mcfb ,
							    int wR3B )
{
 // structure of the graph (remember that node names start from 1):
 //
 // - nodes 1 .. f_n_facilities: facility nodes
 // - nodes f_n_facilities + 1 ... f_n_facilities + f_n_customers:
 //   customer nodes
 // - node f_n_facilities + f_n_customers + 1: super-source
 //
 // - arcs 0 ... f_n_facilities: from super-source to facility 
 // - arcs f_n_facilities ... f_n_facilities * ( f_n_customers + 1 ) - 1:
 //   from source to facility, arranged facility-wise (first f_n_customers
 //   arcs from 1st facility, then f_n_customers arcs from 2nd facility ...)
 //
 // optionally, arcs f_n_facilities * ( f_n_customers + 1 ) ...
 // f_n_facilities * ( f_n_customers + 1 ) + f_n_customers - 1:
 // from super-source directly to customers (huge cost, +INF capacity) to
 // ensure the MCF is never empty

 Index NN = f_n_facilities + f_n_customers + 1;
 Index NA = f_n_facilities * ( f_n_customers + 1 );
 if( wR3B > 1 )
  NA += f_n_customers;

 Subset EN( NA );
 Subset SN( NA );
 MCFBlock::Vec_FNumber U( NA );
 MCFBlock::Vec_CNumber C( NA );
 MCFBlock::Vec_FNumber B( NN );

 // construct deficits vector
 Index i = 0;
 while( i < f_n_facilities )  // facilities nodes
  B[ i++ ] = 0;

 MCFBlock::FNumber todD = 0;
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
  SN[ a ] = ss;
  EN[ a ] = i + 1;
  U[ a ] = v_capacity[ i ];
  // arcs corresponding to fixed-open facilities have 0 cost
  C[ a++ ] = v_fxd[ i ] != yFxd1 ? v_f_cost[ i ] / v_capacity[ i ] : 0;
  }

 // now the facility -> customers arcs
 for( i = 0 ; i < f_n_facilities ; ++i )
  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   SN[ a ] = i + 1;
   EN[ a ] = f_n_facilities + 1 + j;
   U[ a ] = Inf< MCFBlock::FNumber >();
   C[ a++ ] = v_t_cost[ i ][ j ] / v_demand[ j ];
   }

 if( wR3B > 1 ) {
  // now the artificial arcs to ensure feasibility

  for( Index j = 0 ; j < f_n_customers ; ++j ) {
   SN[ a ] = ss;
   EN[ a ] = f_n_facilities + 1 + j;
   U[ a ] = Inf< MCFBlock::FNumber >();

   // compute an upper bound on the worst-case transportation cost
   MCFBlock::CNumber maxc = 0;
   for( i = 0 ; i < f_n_facilities ; ++i )
    if( auto tci = C[ i ] + v_t_cost[ i ][ j ] / v_demand[ j ] ;
	maxc < tci )
     maxc = tci;
   maxc += 1;    // ! +1
   maxc *= 100;  // ! *100 
   C[ a++ ] = maxc;
   }
  }

 mcfb->load( NN , NA , EN , SN , U , C , B );

 if( auto nf = std::count_if( v_fxd.begin() ,  v_fxd.end() ,
			      [] ( auto el ) { return( el == yFxd0 ); } ) ) {
  Subset tfx( nf );
  auto tfxit = tfx.begin();
  for( Index i = 0 ; i < f_n_facilities ; ++i )
   if( v_fxd[ i ] == yFxd0 )
    *(tfxit++) = i;
  mcfb->close_arcs( std::move( tfx ) , true , eNoMod , eNoMod );
  }
 }  // end( CapacitatedFacilityLocationBlock::guts_of_get_R3B_MCF )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_chg_tcost_MCF(
				     MCFBlock * mcfb , Range rng ,
				     ModParam issueMod , ModParam issueAMod )
{
 c_Index f = rng.first;
 c_Index s = rng.second;

 if( s == f + 1 ) {
  mcfb->chg_cost( (v_t_cost.data())[ f ] / v_demand[ f % f_n_customers ] ,
		  f + f_n_facilities , issueMod , issueAMod );
  return;
  }

 CVector NSC( s - f );
 auto NSCit = NSC.begin();
 auto tcit =  v_t_cost.data() + f;
 for( Index h = f ; h < s ; ++h )
  *(NSCit++) = *(tcit++) / v_demand[ h % f_n_customers ];

 rng.first += f_n_facilities;
 rng.second += f_n_facilities;
 
 mcfb->chg_costs( NSC.begin() , rng , issueMod , issueAMod );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_chg_tcost_MCF( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_chg_tcost_MCF(
			    MCFBlock * mcfb , c_Subset & nms , bool ordered ,
			    ModParam issueMod , ModParam issueAMod )
{
 if( nms.size() == 1 ) {
  Index f = nms.front();
  mcfb->chg_cost( (v_t_cost.data())[ f ] / v_demand[ f % f_n_customers ] ,
		  f + f_n_facilities , issueMod , issueAMod );
  return;
  }

 CVector NSC( nms.size() );
 Subset nnms( nms.size() );
 auto NSCit = NSC.begin();
 auto nnmsit = nnms.begin();
 for( Index h : nms ) {
  *(NSCit++) = (v_t_cost.data())[ h ] / v_demand[ h % f_n_customers ];
  *(nnmsit++) = h + f_n_facilities;
  }

 mcfb->chg_costs( NSC.begin() , std::move( nnms ) , ordered ,
		  issueMod , issueAMod );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_chg_tcost_MCF( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_chg_dem_MCF( MCFBlock * mcfb ,
			 Range rng , ModParam issueMod , ModParam issueAMod )
{
 // this operation is complicated by the fact that the demands scale the
 // transportation costs, i.e., the (unitary flow) cost of the transportation
 // arc ( i , j ) is v_t_cost[ i ][ h ] / v_demand[ j ]

 // there will be three Modification, one for changing costs
 not_ModBlock( issueAMod );
 auto iAM = open_if_needed( issueAMod , 3 );
 f_mod_skip = true;

 // change the demands: these are the deficits of the corresponding demand
 // nodes plus the deficit of the super-source, to ensure that the sum of
 // all deficits always remains == 0 
 if( rng.second == rng.first + 1 ) 
  mcfb->chg_dfct( v_demand[ rng.first ] , f_n_facilities + rng.first ,
		  issueMod , iAM );
 else
  mcfb->chg_dfcts( v_demand.begin() + rng.first ,
		   Range( rng.first + f_n_facilities ,
			  rng.second + f_n_facilities ) , issueMod , iAM );

 // note the "Demand( 0 )": without it, std::accumulate() may decide to
 // accumulate on the integers, causing unfeasibility
 mcfb->chg_dfct( - std::accumulate( v_demand.begin() , v_demand.end() ,
				    Demand( 0 ) ) ,
		 f_n_facilities + f_n_customers , issueMod , iAM );

 // now update the costs of all arcs ( i , j ) s.t. the capacity of j changed
 Subset nms( ( rng.second - rng.first ) * f_n_facilities );
 auto nmsit = nms.begin();
 for( Index i = 0 ; i < f_n_facilities ; ++i )
  for( Index j = rng.first ; j < rng.second ; ++j )
   *(nmsit++) = i * f_n_customers + j;

 guts_of_chg_tcost_MCF( mcfb , nms , true , issueMod , iAM );

 f_mod_skip = false;
 close_if_needed( iAM , 3 );  // close the new channel
 }

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_chg_dem_MCF( MCFBlock * mcfb ,
		                     c_Subset & nms , bool ordered ,
				     ModParam issueMod , ModParam issueAMod )
{
 // this operation is complicated by the fact that the demands scale the
 // transportation costs, i.e., the (unitary flow) cost of the transportation
 // arc ( i , j ) is v_t_cost[ i ][ h ] / v_demand[ j ]

 // there will be two Modification, one being for changing costs
 not_ModBlock( issueAMod );
 auto iAM = open_if_needed( issueAMod , 2 );
 f_mod_skip = true;

 // change the demands: these are the deficits of the corresponding demand
 // nodes plus the deficit of the super-source, to ensure that the sum of
 // all deficits always remains == 0
 MCFBlock::Vec_FNumber ND( nms.size() + 1 );
 Subset nnms( nms.size() + 1 );
 std::copy( nms.begin() , nms.end() , nnms.begin() );
 auto NDit = ND.begin();
 for( auto & j : nnms ) {
  *(NDit++) = v_demand[ j ];
  j += f_n_facilities;
  }

 // note the "Demand( 0 )": without it, std::accumulate() may decide to
 // accumulate on the integers, causing unfeasibility
 ND.back() = - std::accumulate( v_demand.begin() , v_demand.end() ,
				Demand( 0 ) );
 nnms.back() = f_n_facilities + f_n_customers;

 mcfb->chg_dfcts( ND.begin() , std::move( nnms ) , ordered , issueMod , iAM );

 // now update the costs of all arcs ( i , j ) s.t. the capacity of j changed
 nnms.resize( nms.size() * f_n_facilities );
 auto nnmsit = nnms.begin();
 for( Index i = 0 ; i < f_n_facilities ; ++i )
  for( Index j : nms )
   *(nnmsit++) = i * f_n_customers + j;

 guts_of_chg_tcost_MCF( mcfb , nnms , ordered , issueMod , iAM );

 f_mod_skip = false;
 close_if_needed( iAM , 2 );  // close the new channel
 }

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationSF(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification for the Standard Formulation- - - - - - - -
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

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_Observer() );
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
			 std::distance( & v_cap.front() , cnst ) ,
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
     for( auto h : nms ) {
      *(nnmsit++) = h - f_n_facilities;
      *(NCit++) = (lfo->get_v_var())[ h ].second;
      }
     chg_transportation_costs( NC.begin() , std::move( nnms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }
    }

   return;
   }

  auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_Observer() );
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
			 std::distance( & v_cap.front() , cnst ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
  return;
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const RowConstraintMod * >( mod ) )
   throw( std::invalid_argument( "RowConstraintMod not allowed" ) );

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const VariableMod * >( mod ) ) {
  auto yi = static_cast< ColVariable * >( tmod->variable() );
  if( ( yi < & v_y.front() ) || ( yi > & v_y.back() ) )
   throw( std::invalid_argument( "VariableMod to not design variable" ) );

  c_Index i = std::distance( & v_y.front() , yi );

  auto new_state = tmod->new_state();  // get new state of the variable
  if( ( ! ColVariable::is_unitary( new_state ) ) ||
      ( ! ColVariable::is_positive( new_state ) ) )
   throw( std::invalid_argument( "invalid ColVariable Modification" ) );

  if( Variable::is_fixed( new_state ) )
   if( yi->get_value() == 1 )
    fix_open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
   else
    close_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
  else
   open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationSF )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationKFG(
			      const GroupModification * mod , ChnlName chnl )
{
 // process abstract Modification for the Knapsack Formulation- - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * This method does the part of the processing related to GroupModification
  * that come from some of the BinaryKnapsackBlock. Only "physical
  * Modification" or GroupModification further nested into mod need be
  * considered. */

 for( auto submod : mod->sub_Modifications() )
  if( auto bmod = dynamic_cast< BinaryKnapsackBlockMod * >( submod.get() ) )
   guts_of_add_ModificationKFP( bmod , chnl );
  else
   if( auto gmod = dynamic_cast< GroupModification * >( submod.get() ) )
    guts_of_add_ModificationKFG( gmod , chnl );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationKFG )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationKFP(
			 const BinaryKnapsackBlockMod * mod , ChnlName chnl )
{
 // process abstract Modification for the Knapsack Formulation- - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock in the "Knapsack
  * Formulation" is *not* a "leaf" Block, i.e., it has sub-Block, one must
  * deal with GroupModification, since these can be produced by
  * Block::add_Modification() in the sub-Block. While this is not directly
  * dealt with here (but in the *G version of the method), the consequence
  * is that GroupModification may introduce arbitrary delay between the
  * moment in which the Modification is produced and the one in which it is
  * processed, which would in principle complicate the logic. However
  *
  *     CapacitatedFacilityLocationBlock IS A "STATIC" Block IN WHICH THE
  *     SIZE OF THE STUFF NEVER CHANGES (save if it is re-loaded whole)
  *
  * This means that the indices, ranges and subsets found in the Modification
  * are always still valid, which drastically simplifies some of the logic. */

 // find the originating BinaryKnapsackBlock  - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto bkb = dynamic_cast< BinaryKnapsackBlock * >( mod->get_Block() );
 auto it = std::lower_bound( v_Block.begin() , v_Block.end() , bkb );
 if( *it != bkb )
  throw( std::invalid_argument(
			 "Modification from unknown BinaryKnapsackBlock" ) );

  c_Index i = std::distance( v_Block.begin() , it ); 

 // BinaryKnapsackBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the only supported changes are those of the profits and of the weight of
 // the last item  (opposite of capacity)
 if( auto tmod = dynamic_cast< const BinaryKnapsackBlockRngdMod * >( mod ) ) {
  c_Index f = tmod->rng().first;
  c_Index s = tmod->rng().second;

  switch( tmod->type() ) {
   case( BinaryKnapsackBlockMod::eChgWeight ):  //- - - - - - - - - - - - - -

    if( ( f == f_n_customers ) && ( s == f + 1 ) ) {
     // note that the coefficient of y is the opposite of the capacity
     chg_facility_capacity( - bkb->get_Weight( f ) , i ,
			    make_par( eNoBlck , chnl ) , eDryRun );
     return;
     }
    throw( std::invalid_argument(
		"unsupported weight Modification in BinaryKnapsackBlock" ) );

   case( BinaryKnapsackBlockMod::eChgProfit ):  //- - - - - - - - - - - - - -

    if( f < f_n_customers ) {  // transportation costs are modified
     c_Index end = std::min( s , f_n_customers );
     c_Index sz = s - f;
     c_Index offst = f_n_customers * i;
     if( sz == 1 )               // one pair
      chg_transportation_cost( bkb->get_Profit( f ) , f + offst ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     else {                      // many pairs
      CVector NC( bkb->get_Profits().begin() + f ,
		  bkb->get_Profits().begin() + end );
      chg_transportation_costs( NC.begin() ,
				Range( f + offst , end + offst ) ,
				make_par( eNoBlck , chnl ) , eDryRun );
      }
     }

    if( s > f_n_customers )  // the facility cost is modified
     chg_facility_cost( bkb->get_Profit( f_n_customers ) , i ,
			make_par( eNoBlck , chnl ) , eDryRun );
    return;

   case( BinaryKnapsackBlockMod::eFixX ):  // - - - - - - - - - - - - - - - -

    if( f < f_n_customers )
     throw( std::invalid_argument( "unsupported variable fixing" ) );

    if( bkb->get_Var( f_n_customers )->get_value() == 1 )
     fix_open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
    else
     close_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

    return;

   case( BinaryKnapsackBlockMod::eUnfixX ):  // - - - - - - - - - - - - - - -

    if( f < f_n_customers )
     throw( std::invalid_argument( "unsupported variable unfixing" ) );

    open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

    return;

   default:  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    throw( std::invalid_argument( "unsupported BinaryKnapsackBlockRngdMod" ) );

   }  // end( switch )
  }  // end( BinaryKnapsackBlockRngdMod )

 // BinaryKnapsackBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the only supported changes are those of the profits and of the weight of
 // the last item  (opposite of capacity)
 if( auto tmod = dynamic_cast< const BinaryKnapsackBlockSbstMod * >( mod ) ) {
  auto & nms = tmod->nms();

  switch( tmod->type() ) {
   case( BinaryKnapsackBlockMod::eChgWeight ):  //- - - - - - - - - - - - - -

    if( ( nms.back() == f_n_customers ) && ( nms.size() == 1 ) ) {
     // note that the coefficient of y is the opposite of the capacity
     chg_facility_capacity( - bkb->get_Weight( f_n_customers ) , i ,
			    make_par( eNoBlck , chnl ) , eDryRun );
     return;
     }
    throw( std::invalid_argument(
		"unsupported weight Modification in BinaryKnapsackBlock" ) );

   case( BinaryKnapsackBlockMod::eChgProfit ):  //- - - - - - - - - - - - - -

    if( nms.front() < f_n_customers ) {   // transportation costs changed
     c_Index offst = f_n_customers * i;
     Subset nnms( nms.begin() , nms.back() == f_n_customers ?
		                std::prev( nms.end() ) : nms.end() );
     if( nms.size() == 1 )               // only one pair
      chg_transportation_cost(  bkb->get_Profit( nms.front() ) ,
				nnms.front() + offst ,
				make_par( eNoBlck , chnl ) , eDryRun );
     else {                      // many pairs
      CVector NC( nnms.size() );
      auto NCit = NC.begin();
      for( auto & j : nnms ) {
       *(NCit++) = (bkb->get_Profits())[ j ];
       j += offst;
       }
      chg_transportation_costs( NC.begin() , std::move( nnms ) , true ,
				make_par( eNoBlck , chnl ) , eDryRun );
      }
     }

    if( nms.back() == f_n_customers )  // the facility cost is modified
     chg_facility_cost( bkb->get_Profit( f_n_customers ) , i ,
			make_par( eNoBlck , chnl ) , eDryRun );
    return;

   case( BinaryKnapsackBlockMod::eFixX ):  // - - - - - - - - - - - - - - - -

    if( ( nms.size() != 1 ) || ( nms.back() < f_n_customers ) )
     throw( std::invalid_argument( "unsupported variable fixing" ) );

    if( bkb->get_Var( f_n_customers )->get_value() == 1 )
     fix_open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
    else
     close_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

    return;

   case( BinaryKnapsackBlockMod::eUnfixX ):  // - - - - - - - - - - - - - - -

    if( ( nms.size() != 1 ) || ( nms.back() < f_n_customers ) )
     throw( std::invalid_argument( "unsupported variable unfixing" ) );

    open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

    return;

   default:  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    throw( std::invalid_argument( "unsupported BinaryKnapsackBlockSbstMod" ) );

   }  // end( switch )
  }  // end( BinaryKnapsackBlockSbstMod )
 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationKFP )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFA(
						c_p_Mod mod , ChnlName chnl )
{
 // process abstract Modification for the Flow Formulation- - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * This method only deals with the "abstract Modification" from either the
  * "root" CapacitatedFacilityLocationBlock or the "design AbstractBlock",
  * which arrive with concerns_Block() == true. This coming from the "root"
  * Block do so "directly", and therefore cannot be GroupModification, but
  * those from the AbstractBlock pass from Block::add_Modification() and
  * therefore can be GroupModification. This introduces delay between the
  * moment in which the Modification is produced and the one in which it is
  * processed, which would in principle complicate the logic. However
  *
  *     CapacitatedFacilityLocationBlock IS A "STATIC" Block IN WHICH THE
  *     SIZE OF THE STUFF NEVER CHANGES (save if it is re-loaded whole)
  *
  * This means that the indices, ranges and subsets found in the Modification
  * are always still valid, which drastically simplifies some of the logic. */

 // GroupModification - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that a GroupModification has concerns_Block() == true is any of its
 // sub-Modification has it, which means that not all the sub-Modification
 // have; it since this method only deals with Modification having
 // concerns_Block() == true, the others are skipped
 if( auto gmod = dynamic_cast< const GroupModification * >( mod ) ) {
  for( auto submod : gmod->sub_Modifications() )
   if( submod->concerns_Block() )
    guts_of_add_ModificationFFA( submod.get() , chnl );

  return;
  }
 
 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinRngd * >( mod ) ) {
  c_Index f = tmod->range().first;
  c_Index s = tmod->range().second;
  auto lfo = LF( tmod->function() );

  if( tmod->get_Block() == this ) {
   // Modification in the "root" CapacitatedFacilityLocationBlock- - - - - - -
   // it can only be in the facility capacity constraints
   // note: the facility capacity is actually found in *two* different places
   //       in the Flow Formulation, the linking constraints in the "root"
   //       Block (dealt with here) and the capacities of the facility arcs in
   //       the MCFBlock (dealt somewhere else). both Modification are
   //       dealt with independently, which is a bit wasteful but the second
   //       call to chg_facility_capacity() will do nothing as the capacity
   //       value is the same, so it's acceptable

   auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_Observer() );
   if( ! cnst )
    throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

   if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
    throw( std::invalid_argument(
	             "Modification to FRowConstraint not capacity one" ) );

   if( ( f != 1 ) || ( s != 2 ) )
    throw( std::invalid_argument(
	    "Modification to wrong coefficient in capacity constraint" ) );

   // note that the coefficient of y is the opposite of the capacity
   chg_facility_capacity( - (lfo->get_v_var())[ 1 ].second ,
			  std::distance( & v_cap.front() , cnst ) ,
			  make_par( eNoBlck , chnl ) , eDryRun );   
   return;
   }

  if( tmod->get_Block() == v_Block.front() ) {
   // Modification in the "design" AbstractBlock - - - - - - - - - - - - - - -
   if( LF( f_obj.get_function() ) == lfo ) {  // Modification to the Objective
    c_Index sz = s - f;
    if( sz == 1 )              // one facility
     chg_facility_cost( (lfo->get_v_var())[ f ].second , f ,
			make_par( eNoBlck , chnl ) , eDryRun );
    else {                     // many facilities
     CVector NC( sz );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; )
      *(NCit++) = (lfo->get_v_var())[ i++ ].second;
     chg_facility_costs( NC.begin() , Range( f , s ) ,
			 make_par( eNoBlck , chnl ) , eDryRun );
     }

    return;
    }

   throw( std::invalid_argument(
		       "unsupported Modification to design AbstractBlock" ) );
   }

  // then it must be a Modification in the MCFBlock- - - - - - - - - - - - - -
  // ... but no "abstract Modification" should ever reach here, ignore it
  return;
  }

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const C05FunctionModLinSbst * >( mod ) ) {
  auto & nms = tmod->subset();
  auto lfo = LF( tmod->function() );

  if( tmod->get_Block() == this ) {
   // Modification in the "root" CapacitatedFacilityLocationBlock- - - - - - -
   // it can only be in the facility capacity constraints

   auto cnst = dynamic_cast< FRowConstraint * >( lfo->get_Observer() );
   if( ! cnst )
    throw( std::invalid_argument( "Modification to not FRowConstraint" ) );

   if( ( cnst < & v_cap.front() ) || ( cnst > & v_cap.back() ) )
    throw( std::invalid_argument(
	             "Modification to FRowConstraint not capacity one" ) );

   if( ( nms.size() != 1 ) || ( nms.front() != 1 ) )
    throw( std::invalid_argument(
	    "Modification to wrong coefficient in capacity constraint" ) );

   // note that the coefficient of y is the opposite of the capacity
   chg_facility_capacity( - (lfo->get_v_var())[ 1 ].second ,
			  std::distance( & v_cap.front() , cnst ) ,
			  make_par( eNoBlck , chnl ) , eDryRun );
   return;
   }

  if( tmod->get_Block() == v_Block.front() ) {
   // Modification in the "design" AbstractBlock - - - - - - - - - - - - - - -
   if( LF( f_obj.get_function() ) == lfo ) {  // Modification to the Objective
    if( nms.size() == 1 )       // one pair
     chg_transportation_cost( (lfo->get_v_var())[ nms.front() ].second ,
			      nms.front() ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                      // many pairs
     CVector NC( nms.size() );
     auto NCit = NC.begin();
     for( auto h : nms )
      *(NCit++) = (lfo->get_v_var())[ h ].second;
     chg_transportation_costs( NC.begin() , Subset( nms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }

    return;
    }

   throw( std::invalid_argument(
		       "unsupported Modification to design AbstractBlock" ) );
   }

  // then it must be a Modification in the MCFBlock- - - - - - - - - - - - - -
  // ... but no "abstract Modification" should ever reach here, ignore it
  return;
  }

 // RowConstraintMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dynamic_cast< const RowConstraintMod * >( mod ) )
  throw( std::invalid_argument( "RowConstraintMod not allowed" ) );

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is the same code as for the Standard Formulation, except that the
 // design variables now live in the "design AbstractBlock"
 if( auto tmod = dynamic_cast< const VariableMod * >( mod ) ) {
  auto yi = static_cast< ColVariable * >( tmod->variable() );
  if( ( yi < & v_y.front() ) || ( yi > & v_y.back() ) )
   throw( std::invalid_argument( "VariableMod to not design variable" ) );

  c_Index i = std::distance( & v_y.front() , yi );

  auto new_state = tmod->new_state();  // get new state of the variable
  if( ( ! ColVariable::is_unitary( new_state ) ) ||
      ( ! ColVariable::is_positive( new_state ) ) )
   throw( std::invalid_argument( "invalid ColVariable Modification" ) );

  if( Variable::is_fixed( new_state ) )
   if( yi->get_value() == 1 )
    fix_open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
   else
    close_facility( i , make_par( eNoBlck , chnl ) , eDryRun );
  else
   open_facility( i , make_par( eNoBlck , chnl ) , eDryRun );

  return;
  }

 throw( std::invalid_argument(
	   "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFA )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFG(
			      const GroupModification * mod , ChnlName chnl )
{
 // process abstract Modification for the Flow Formulation- - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * This method does the part of the processing related to GroupModification
  * that come from some the MCFBlock. Only "physical Modification" or
  * GroupModification further nested into mod need be considered. */

 for( auto submod : mod->sub_Modifications() )
  if( auto mmod = dynamic_cast< MCFBlockMod * >( submod.get() ) )
   guts_of_add_ModificationFFP( mmod , chnl );
  else
   if( auto gmod = dynamic_cast< GroupModification * >( submod.get() ) )
    guts_of_add_ModificationFFG( gmod , chnl );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFG )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFP(
			           const MCFBlockMod * mod , ChnlName chnl )
{
 // process abstract Modification for the Flow Formulation- - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since CapacitatedFacilityLocationBlock in the "Flow
  * Formulation" is *not* a "leaf" Block, i.e., it has sub-Block, one must
  * deal with GroupModification, since these can be produced by
  * Block::add_Modification() in the sub-Block. While this is not directly
  * dealt with here (but in the *G version of the method), the consequence
  * is that GroupModification may introduce arbitrary delay between the
  * moment in which the Modification is produced and the one in which it is
  * processed, which would in principle complicate the logic. However
  *
  *     CapacitatedFacilityLocationBlock IS A "STATIC" Block IN WHICH THE
  *     SIZE OF THE STUFF NEVER CHANGES (save if it is re-loaded whole)
  *
  * This means that the indices, ranges and subsets found in the Modification
  * are always still valid, which drastically simplifies some of the logic. */

 // MCFBlockRngdMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the only supported changes are those of the costs (transportation ones,
 // since facility arcs are supposed to remain 0-cost), of the facility
 // capacities (capacity of facility arcs) and of the customers demands
 // (deficits of customers nodes). note that:
 // - the facility capacity is actually found in *two* different places in
 //   the Flow Formulation, the linking constraints in the "root" Block
 //   (dealt with somewhere else) and the capacities of the facility arcs in
 //   the MCFBlock (dealt )with here. both Modification are dealt with
 //   independently, which is a bit wasteful but the second call to
 //   chg_facility_capacity() will do nothing as the capacity value is the
 //   same, so it's acceptable
 // - we expect changes of node deficits to eventually maintain the
 //   zero-total property necessary for feasibility, but this may be true
 //   only at the end of the series of changes rather than at any point;
 //   hence, we disregard any change in the deficit of the super-source
 //   assuming that eventually it'll be what it necessarily needs be
 if( auto tmod = dynamic_cast< const MCFBlockRngdMod * >( mod ) ) {
  c_Index f = tmod->rng().first;
  Index s = tmod->rng().second;

  switch( tmod->type() ) {
   case( MCFBlockMod::eChgCost ): {  // - - - - - - - - - - - - - - - - - - -
    if( f < f_n_facilities )
     throw( std::invalid_argument(
			       "unsupported arc cost change in MCFBlock" ) );

    if( s - f == 1 )  // one pair
     chg_transportation_cost( MCFB( v_Block.back() )->get_C( f ) ,
			      f - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else              // many pairs
     // note that the cost vector is supposed to exist (costs not all 0)
     chg_transportation_costs( MCFB( v_Block.back() )->get_C().begin() + f ,
			       Range( f - f_n_facilities ,
				      s - f_n_facilities ) ,
			       make_par( eNoBlck , chnl ) , eDryRun );
    return;
    }

   case( MCFBlockMod::eChgCaps ): {  // - - - - - - - - - - - - - - - - - - -
    if( s > f_n_facilities )
     throw( std::invalid_argument(
			   "unsupported arc capacity change in MCFBlock" ) );

    if( s - f == 1 )  // one facility
     chg_facility_capacity( MCFB( v_Block.back() )->get_U( f ) , f ,
			    make_par( eNoBlck , chnl ) , eDryRun );
    else              // many facilities
     // note that the capacity vector is supposed to exist (not all +INF)
     chg_facility_capacities( MCFB( v_Block.back() )->get_U().begin() + f ,
			      tmod->rng() ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    return;
    }

   case( MCFBlockMod::eChgDfct ): {  // - - - - - - - - - - - - - - - - - - -
    if( f < f_n_facilities )
     throw( std::invalid_argument(
			   "unsupported node deficit change in MCFBlock" ) );

    if( s > f_n_facilities + f_n_customers )
     s = f_n_facilities + f_n_customers;  // ignore changes of the last node

    if( s - f == 1 )  // one customer
     chg_customer_demand( MCFB( v_Block.back() )->get_B( f ) ,
			  f - f_n_facilities ,
			  make_par( eNoBlck , chnl ) , eDryRun );
    else              // many customers
     // note that the deficits vector is supposed to exist (not all 0)
     chg_customer_demands( MCFB( v_Block.back() )->get_B().begin() + f ,
			   Range( f - f_n_facilities , s - f_n_facilities ) ,
			   make_par( eNoBlck , chnl ) , eDryRun );
    return;
    }

   default:  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    throw( std::invalid_argument( "unsupported MCFBlockRngdMod" ) );

   }  // end( switch )
  }  // end( MCFBlockRngdMod )

 // MCFBlockSbstMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const MCFBlockSbstMod * >( mod ) ) {
  auto & nms = tmod->nms();

  switch( tmod->type() ) {
   case( MCFBlockMod::eChgCost ): {  // - - - - - - - - - - - - - - - - - - -
    if( nms.front() < f_n_facilities )
     throw( std::invalid_argument(
			       "unsupported arc cost change in MCFBlock" ) );

    if( nms.size() == 1 )  // one pair
     chg_transportation_cost( MCFB( v_Block.back() )->get_C( nms.front() ) ,
			      nms.front() - f_n_facilities ,
			      make_par( eNoBlck , chnl ) , eDryRun );
    else {                 // many pairs
     // note that the cost vector is supposed to exist (costs not all 0)
     CVector NC( nms.size() );
     auto NCit = NC.begin();
     for( auto h : nms )
      *(NCit++) = ( MCFB( v_Block.back() )->get_C() )[ h ];
     Subset nnms( nms.begin() , nms.end() );
     for( auto & h : nnms )
      h -= f_n_facilities;
     chg_transportation_costs( NC.begin() , std::move( nnms ) , true ,
			       make_par( eNoBlck , chnl ) , eDryRun );
     }

    return;
    }

   case( MCFBlockMod::eChgCaps ): {  // - - - - - - - - - - - - - - - - - - -
    if( nms.back() >= f_n_facilities )
     throw( std::invalid_argument(
			   "unsupported arc capacity change in MCFBlock" ) );

    if( nms.size() == 1 )  // one facility
     chg_facility_capacity( MCFB( v_Block.back() )->get_U( nms.front() ) ,
			    nms.front() ,
			    make_par( eNoBlck , chnl ) , eDryRun );
    else {                 // many facilities
     // note that the capacity vector is supposed to exist (not all +INF)
     CVector NU( nms.size() );
     auto NUit = NU.begin();
     for( auto i : nms )
      *(NUit++) = ( MCFB( v_Block.back() )->get_U() )[ i ];
     chg_facility_capacities( NU.begin() , Subset( nms ) , true ,
			      make_par( eNoBlck , chnl ) , eDryRun );
      }

    return;
    }

   case( MCFBlockMod::eChgDfct ): {  // - - - - - - - - - - - - - - - - - - -
    if( nms.front() < f_n_facilities )
     throw( std::invalid_argument(
			   "unsupported node deficit change in MCFBlock" ) );

    if( nms.back() == f_n_facilities + f_n_customers ) {
     // changing also the last deficit, which must be ignored
     if( nms.size() == 1 )  // just the one
      return;

     if( nms.size() == 2 )  // one customer
      chg_customer_demand( MCFB( v_Block.back() )->get_B( nms.front() ) ,
			   nms.front() - f_n_facilities ,
			   make_par( eNoBlck , chnl ) , eDryRun );
     else {                 // many customers
      // note that the deficits vector is supposed to exist (not all 0)
      CVector NB( nms.size() - 1 );
      Subset nnms( nms.size() - 1 );
      auto NBit = NB.begin();
      auto nnmsit = nnms.begin();
      for( auto nmsit = nms.begin() ; nmsit != std::prev( nms.end() ) ; ) {
       *(NBit++) = ( MCFB( v_Block.back() )->get_B() )[ *nmsit ];
       *(nnmsit++) = *(nmsit++) - f_n_facilities;
       }
      chg_customer_demands( NB.begin() , std::move( nnms ) , true ,
			    make_par( eNoBlck , chnl ) , eDryRun );
      }
     }
    else  // not changing the last deficit
     if( nms.size() )  // one customer
      chg_customer_demand( MCFB( v_Block.back() )->get_B( nms.front() ) ,
			   nms.front() - f_n_facilities ,
			   make_par( eNoBlck , chnl ) , eDryRun );
     else {            // many customers
      // note that the deficits vector is supposed to exist (not all 0)
      CVector NB( nms.size() );
      auto NBit = NB.begin();
      for( auto i : nms )
       *(NBit++) = ( MCFB( v_Block.back() )->get_B() )[ i ];
      Subset nnms( nms.begin() , nms.end() );
      for( auto & h : nnms )
       h -= f_n_facilities;
      chg_customer_demands( NB.begin() , Subset( nms ) , true ,
			    make_par( eNoBlck , chnl ) , eDryRun );
      }

    return;
    }

   default:  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    throw( std::invalid_argument( "unsupported MCFBlockRngdMod" ) );

   }  // end( switch )
  }  // end( MCFBlockSbstMod )
 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_ModificationFFP )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_copy(
		       CapacitatedFacilityLocationBlock * R3B , c_p_Mod mod ,
		       ModParam issuePMod , ModParam issueAMod )
{
 /* The (Solver attached to the) R3Block may require some Modification to be
  * "neatly packaged" into appropriate GroupModification to work
  * (efficiently): hence, we make an effort to group the changes produced by
  * mod, if it's a GroupModification, by opening a new channel if the original
  * one is the default one or nesting it if it is already a GroupModification
  * one. Note that this is done independently for "physical Modification" and
  * "abstract Modification". */

 bool ok = true;  // final return value

 if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
  // open / nest two new the channels
  auto iPM = make_par( par2mod( issuePMod ) ,
		       R3B->open_channel( par2chnl( issuePMod ) ) );
  auto iPA = make_par( par2mod( issueAMod ) ,
		       R3B->open_channel( par2chnl( issueAMod ) ) );
  
  // run through each sub-Mod on these channels
  for( const auto & submod : tmod->sub_Modifications() )
   if( ! guts_of_map_f_Mod_copy( R3B , submod.get() , iPM , iPA ) )
    ok = false;

  // close / un-nest newly the opened channels
  R3B->close_channel( par2chnl( iPA ) );
  R3B->close_channel( par2chnl( iPM ) );
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
 // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this requires to patiently sift through the possible Modification types
 // to find what this Modification exactly is, and call the appropriate
 // method of R3B: note that we only consider "physical Modification" since
 // any change in the "abstract representation" is intercepted by the :Block
 // and a "physical Modification" is produced, so intercepting "abstract
 // Modification" is useless and wasteful (besides being more complicated)

 // CapacitatedFacilityLocationBlockRngdMod - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast<
                 const CapacitatedFacilityLocationBlockRngdMod * >( mod ) ) {
  c_Index f = tmod->rng().first;
  c_Index s = tmod->rng().second;

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgFCost ):  //- - - - - - - -
    if( s == f + 1 )
     R3B->chg_facility_cost( v_f_cost[ f ] , f , issuePMod , issueAMod );
    else
     R3B->chg_facility_costs( v_f_cost.begin() + f , tmod->rng() ,
			      issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgTCost ):  //- - - - - - - -
    if( s == f + 1 )
     R3B->chg_transportation_cost( (v_t_cost.data())[ f ] , f ,
				   issuePMod , issueAMod );
    else{
     CVector NC( s - f );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; ++i )
      *(NCit++) = (v_t_cost.data())[ i ]; 
     R3B->chg_transportation_costs( NC.begin() , tmod->rng() ,
				    issuePMod , issueAMod );
     }
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgCap ):  //- - - - - - - - -
    if( s == f + 1 )
     R3B->chg_facility_capacity( v_capacity[ f ] , f ,
				 issuePMod , issueAMod );
    else
     R3B->chg_facility_capacities( v_capacity.begin() + f , tmod->rng() ,
				   issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgDem ):  //- - - - - - - - -
    if( s == f + 1 )
     R3B->chg_customer_demand( v_demand[ f ] , f , issuePMod , issueAMod );
    else
     R3B->chg_customer_demands( v_demand.begin() + f , tmod->rng() ,
				issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eCloseF ):  //- - - - - - - - -
    if( s == f + 1 )
     R3B->close_facility( f , issuePMod , issueAMod );
    else
     R3B->close_facilities( tmod->rng() , issuePMod , issueAMod );
    break;
   case( CapacitatedFacilityLocationBlockMod::eOpenF ):  // - - - - - - - - -
    if( s == f + 1 )
     R3B->open_facility( f , issuePMod , issueAMod );
    else
     R3B->open_facilities( tmod->rng() , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eBuyF ):  //- - - - - - - - - -
    if( s == f + 1 )
     R3B->fix_open_facility( f , issuePMod , issueAMod );
    else
     R3B->fix_open_facilities( tmod->rng() , issuePMod , issueAMod );
    break;

    default:
     throw( std::invalid_argument(
		   "unknown CapacitatedFacilityLocationBlockRngdMod type" ) );

   }  // end( switch )

  return( true );
  }

 // CapacitatedFacilityLocationBlockSbstMod - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that tmod->nms() need be copied, since the chg_*() methods "consume"
 // the names vector
 if( auto tmod = dynamic_cast<
                 const CapacitatedFacilityLocationBlockSbstMod * >( mod ) ) {

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgFCost ):  //- - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->chg_facility_cost( v_f_cost[ tmod->nms().front() ] ,
			     tmod->nms().front() , issuePMod , issueAMod );
    else {
     CVector NC( tmod->nms().size() );
     auto NCit = NC.begin();
     for( auto i : tmod->nms() )
      *(NCit++) = v_f_cost[ i ]; 
     R3B->chg_facility_costs( NC.begin() , Subset( tmod->nms() ) , true ,
			      issuePMod , issueAMod );
     }
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgTCost ):  //- - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->chg_transportation_cost( (v_t_cost.data())[ tmod->nms().front() ] ,
				   tmod->nms().front() ,
				   issuePMod , issueAMod );
    else {
     CVector NC( tmod->nms().size() );
     auto NCit = NC.begin();
     for( auto i : tmod->nms() )
      *(NCit++) = (v_t_cost.data())[ i ]; 
     R3B->chg_transportation_costs( NC.begin() , Subset( tmod->nms() ) ,
				    true , issuePMod , issueAMod );
     }
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgCap ):  //- - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->chg_facility_capacity( v_capacity[ tmod->nms().front() ] ,
				 tmod->nms().front() ,
				 issuePMod , issueAMod );
    else {
     DVector NC( tmod->nms().size() );
     auto NCit = NC.begin();
     for( auto i : tmod->nms() )
      *(NCit++) = v_capacity[ i ]; 
     R3B->chg_facility_capacities( NC.begin() , Subset( tmod->nms() ) ,
				   true , issuePMod , issueAMod );
     }
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgDem ):  //- - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->chg_customer_demand( v_demand[ tmod->nms().front() ] ,
			       tmod->nms().front() ,
			       issuePMod , issueAMod );
    else {
     DVector ND( tmod->nms().size() );
     auto NDit = ND.begin();
     for( auto i : tmod->nms() )
      *(NDit++) = v_demand[ i ]; 
     R3B->chg_customer_demands( ND.begin() , Subset( tmod->nms() ) , true ,
				issuePMod , issueAMod );
     }
    break;

   case( CapacitatedFacilityLocationBlockMod::eCloseF ):  //- - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->close_facility( tmod->nms().front() , issuePMod , issueAMod );
    else
     R3B->close_facilities( Subset( tmod->nms() ) , true ,
			    issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eOpenF ):  // - - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->open_facility( tmod->nms().front() , issuePMod , issueAMod );
    else
     R3B->open_facilities( Subset( tmod->nms() ) , true ,
			   issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eBuyF ):  //- - - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->fix_open_facility( tmod->nms().front() , issuePMod , issueAMod );
    else
     R3B->fix_open_facilities( Subset( tmod->nms() ) , true ,
			       issuePMod , issueAMod );
    break;

   default:
     throw( std::invalid_argument(
		   "unknown CapacitatedFacilityLocationBlockRngdMod type" ) );

   }  // end( switch )

  return( true );
  }

 // CapacitatedFacilityLocationBlockMod - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // since the *Rngd and *Sbst versions derive from
 // CapacitatedFacilityLocationBlockMod this dynamic_cast<> would succeed on
 // these; but here we only want to catch the base class, so this has to be
 // done after the derived classes
 
 if( auto tmod = dynamic_cast<
                    const CapacitatedFacilityLocationBlockMod * >( mod ) ) {

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgUnSplt ):  //- - - - - - -
    R3B->chg_UnSplittable( true , issuePMod , issueAMod );
    break;
   case( CapacitatedFacilityLocationBlockMod::eChgSplt ):  //- - - - - - - -
    R3B->chg_UnSplittable( false , issuePMod , issueAMod );
    break;
   default:
    throw( std::invalid_argument(
		  "invalid type in CapacitatedFacilityLocationBlockMod" ) );
   }

  return( true );
  }

 // NBModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is the "nuclear option": the CapacitatedFacilityLocationBlock has
 // been re-loaded

 if( auto tmod = dynamic_cast< const NBModification * >( mod ) ) {
  R3B->load( f_n_facilities , f_n_customers , v_capacity , v_f_cost ,
	     v_demand , v_t_cost );
  return( true );
  }

 return( false );  // any other Modification is not mapped

 }  // end( CapacitatedFacilityLocationBlock::guts_of_guts_of_map_f_Mod_copy )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::guts_of_map_f_Mod_MCF(
				    MCFBlock * R3B , c_p_Mod mod ,
			            ModParam issuePMod , ModParam issueAMod )
{
 /* The (Solver attached to the) R3Block may require some Modification to be
  * "neatly packaged" into appropriate GroupModification to work
  * (efficiently): hence, we make an effort to group the changes produced by
  * mod, if it's a GroupModification, by opening a new channel if the original
  * one is the default one or nesting it if it is already a GroupModification
  * one. Note that this is done independently for "physical Modification" and
  * "abstract Modification". */

 bool ok = true;  // final return value

 if( auto tmod = dynamic_cast< const GroupModification * >( mod ) ) {
  // open / nest two new the channels
  auto iPM = make_par( par2mod( issuePMod ) ,
		       R3B->open_channel( par2chnl( issuePMod ) ) );
  auto iPA = make_par( par2mod( issueAMod ) ,
		       R3B->open_channel( par2chnl( issueAMod ) ) );
  
  // run through each sub-Mod on these channels
  for( const auto & submod : tmod->sub_Modifications() )
   if( ! guts_of_map_f_Mod_MCF( R3B , submod.get() , iPM , iPA ) )
    ok = false;

  // close / un-nest newly the opened channels
  R3B->close_channel( par2chnl( iPA ) );
  R3B->close_channel( par2chnl( iPM ) );
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
 // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this requires to patiently sift through the possible Modification types
 // to find what this Modification exactly is, and call the appropriate
 // method of R3B: note that we only consider "physical Modification" since
 // any change in the "abstract representation" is intercepted by the :Block
 // and a "physical Modification" is produced, so intercepting "abstract
 // Modification" is useless and wasteful (besides being more complicated)

 // CapacitatedFacilityLocationBlockRngdMod - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast<
                 const CapacitatedFacilityLocationBlockRngdMod * >( mod ) ) {
  c_Index f = tmod->rng().first;
  c_Index s = tmod->rng().second;

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgFCost ):  //- - - - - - - -
    if( s == f + 1 )
     R3B->chg_cost( v_f_cost[ f ] / v_capacity[ f ] , f ,
		    issuePMod , issueAMod );
    else {
     MCFBlock::Vec_CNumber NC( s - f );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; ++i )
      *(NCit++) = v_f_cost[ i ] / v_capacity[ i ];
     R3B->chg_costs( NC.begin() , tmod->rng() , issuePMod , issueAMod );
     }

    break;

   case( CapacitatedFacilityLocationBlockMod::eChgTCost ):  //- - - - - - - -
    guts_of_chg_tcost_MCF( R3B , tmod->rng() , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgCap ):  //- - - - - - - - -
    if( s == f + 1 ) {
     R3B->chg_ucap( v_capacity[ f ] , f , issuePMod , issueAMod );
     R3B->chg_cost( v_f_cost[ f ] / v_capacity[ f ] , f ,
		    issuePMod , issueAMod );
      }
    else {
     R3B->chg_ucaps( v_capacity.begin() + f , tmod->rng() ,
		     issuePMod , issueAMod );
     MCFBlock::Vec_CNumber NC( s - f );
     auto NCit = NC.begin();
     for( Index i = f ; i < s ; ++i )
      *(NCit++) = v_f_cost[ i ] / v_capacity[ i ];
     R3B->chg_costs( NC.begin() , tmod->rng() , issuePMod , issueAMod );
     }

    break;

   case( CapacitatedFacilityLocationBlockMod::eChgDem ):  //- - - - - - - - -
    guts_of_chg_dem_MCF( R3B , tmod->rng() , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eCloseF ):  //- - - - - - - - -
    if( s == f + 1 )
     R3B->close_arc( f , issuePMod , issueAMod );
    else
     R3B->close_arcs( tmod->rng() , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eOpenF ):  // - - - - - - - - -
    if( s == f + 1 )
     R3B->open_arc( f , issuePMod , issueAMod );
    else
     R3B->open_arcs( tmod->rng() , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eBuyF ):  //- - - - - - - - - -
    throw( std::invalid_argument(
     "mapping fix_open Modification to a MCFBlock R3Block not supported" ) );

    default:
     throw( std::invalid_argument(
		   "unknown CapacitatedFacilityLocationBlockRngdMod type" ) );

   }  // end( switch )

  return( true );
  }

 // CapacitatedFacilityLocationBlockSbstMod - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that tmod->nms() need be copied, since the chg_*() methods "consume"
 // the names vector
 if( auto tmod = dynamic_cast<
                 const CapacitatedFacilityLocationBlockSbstMod * >( mod ) ) {

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgFCost ):  //- - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->chg_cost( v_f_cost[ tmod->nms().front() ] /
		    v_capacity[ tmod->nms().front() ] ,
		    tmod->nms().front() , issuePMod , issueAMod );
    else {
     MCFBlock::Vec_CNumber NC( tmod->nms().size() );
     auto NCit = NC.begin();
     for( auto i : tmod->nms() )
      *(NCit++) = v_f_cost[ i ] / v_capacity[ i ]; 
     R3B->chg_costs( NC.begin() , Subset( tmod->nms() ) , true ,
		     issuePMod , issueAMod );
     }

    break;

   case( CapacitatedFacilityLocationBlockMod::eChgTCost ):  //- - - - - - - -
    guts_of_chg_tcost_MCF( R3B , tmod->nms() , true , issuePMod , issueAMod );
    break;

   case( CapacitatedFacilityLocationBlockMod::eChgCap ):  //- - - - - - - - -
    if( tmod->nms().size() == 1 ) {
     Index f = tmod->nms().front();
     R3B->chg_ucap( v_capacity[ f ] , f , issuePMod , issueAMod );
     R3B->chg_cost( v_f_cost[ f ] / v_capacity[ f ] , f ,
		    issuePMod , issueAMod );
     }
    else {
     MCFBlock::Vec_FNumber NCP( tmod->nms().size() );
     MCFBlock::Vec_CNumber NCS( tmod->nms().size() );
     auto NCPit = NCP.begin();
     auto NCSit = NCS.begin();
     for( auto i : tmod->nms() ) {
      *(NCPit++) = v_capacity[ i ];
      *(NCSit++) = v_f_cost[ i ] / v_capacity[ i ];
      }
     R3B->chg_ucaps( NCP.begin() , Subset( tmod->nms() ) , true ,
		     issuePMod , issueAMod );
     R3B->chg_costs( NCS.begin() , Subset( tmod->nms() ) , true ,
		     issuePMod , issueAMod );
     }

    break;

   case( CapacitatedFacilityLocationBlockMod::eChgDem ): {  //- - - - - - - -
    DVector ND( tmod->nms().size() );
    auto NDit = ND.begin();
    for( auto i : tmod->nms() )
     *(NDit++) = v_demand[ i ]; 

    guts_of_chg_dem_MCF( R3B , tmod->nms() , true , issuePMod , issueAMod );
    break;
    }

   case( CapacitatedFacilityLocationBlockMod::eCloseF ):  //- - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->close_arc( tmod->nms().front() , issuePMod , issueAMod );
    else
     R3B->close_arcs( Subset( tmod->nms() ) , true , issuePMod , issueAMod );

    break;

   case( CapacitatedFacilityLocationBlockMod::eOpenF ):  // - - - - - - - - -
    if( tmod->nms().size() == 1 )
     R3B->open_arc( tmod->nms().front() , issuePMod , issueAMod );
    else
     R3B->open_arcs( Subset( tmod->nms() ) , true , issuePMod , issueAMod );

    break;

   case( CapacitatedFacilityLocationBlockMod::eBuyF ):  //- - - - - - - - - -
    throw( std::invalid_argument(
     "mapping fix_open Modification to a MCFBlock R3Block not supported" ) );

   default:
     throw( std::invalid_argument(
		  "unknown CapacitatedFacilityLocationBlockRngdMod type" ) );

   }  // end( switch )

  return( true );
  }

 // CapacitatedFacilityLocationBlockMod - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // since the *Rngd and *Sbst versions derive from
 // CapacitatedFacilityLocationBlockMod this dynamic_cast<> would succeed on
 // these; but here we only want to catch the base class, so this has to be
 // done after the derived classes
 // the Modification corresponding to changing the type of the problem are
 // ignored here since the MCF R3Block represents a continuous relaxation of
 // the original CFL, and therefore it is identical in the splittable and
 // unsplittable case

 if( auto tmod = dynamic_cast<
                    const CapacitatedFacilityLocationBlockMod * >( mod ) ) {

  switch( tmod->type() ) {
   case( CapacitatedFacilityLocationBlockMod::eChgUnSplt ):  //- - - - - - -
   case( CapacitatedFacilityLocationBlockMod::eChgSplt ):  //- - - - - - - -
    break;    // nothing to do, although it's a weird case
   default:
    throw( std::invalid_argument(
		  "invalid type in CapacitatedFacilityLocationBlockMod" ) );
   }
  }

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
  auto minj = Inf< Cost >();
  auto maxj = -Inf< Cost >();

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
		  typename std::vector< T >::iterator Sol , Range rng ) const
{
 #ifndef NDEBUG
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( "get_facility_solution: variables not generated" ) );
 #endif

 if( rng.second > f_n_facilities )
  rng.second = f_n_facilities;
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( Index i = rng.first ; i < rng.second ; )
   *(Sol++) = BKB( v_Block[ i++ ] )->get_x( f_n_customers );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index i = rng.first ; i < rng.second ; )
  *(Sol++) = v_y[ i++ ].get_value();

 }  // end( get_y( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::get_y(
	       typename std::vector< T >::iterator Sol , c_Subset nms ) const
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "get_facility_solution: ";
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std::any_of( nms.begin() , nms.end() ,
		   [ & ]( auto i ) { return( i > f_n_facilities ); } ) )
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
		  typename std::vector< T >::iterator Sol , Range rng ) const
{
 #ifndef NDEBUG
  if( ! ( AR & HasVar ) )
   throw( std::logic_error(
		  "get_transportation_solution: variables not generated" ) );
 #endif

 if( rng.second > f_n_facilities * f_n_customers )
  rng.second = f_n_facilities * f_n_customers;
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

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
  // it's at least two facilities: deal with the first
  BKB( v_Block[ i++ ] )->get_x( Sol , Range( frst , f_n_customers ) );
  Index k = f_n_customers - frst;
  Index h = rng.first + k;
  Sol += k;
  // now from the second on
  while( h < rng.second ) {
   k = std::min( rng.second - h , f_n_customers );
   BKB( v_Block[ i++ ] )->get_x( Sol , Range( 0 , k ) );
   if( k < f_n_customers )
    break;
   h += f_n_customers;
   Sol += f_n_customers;
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
	       typename std::vector< T >::iterator Sol , c_Subset nms ) const
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "get_transportation_solution: ";
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std::any_of( nms.begin() , nms.end() ,
		   [ & ]( auto i ) { return( i > f_n_facilities ); } ) )
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
		  typename std::vector< T >::const_iterator Sol , Range rng )
{
 #ifndef NDEBUG
  if( ! ( AR & HasVar ) )
   throw( std::logic_error(
		         "set_facility_solution: variables not generated" ) );
 #endif

 if( rng.second > f_n_facilities )
  rng.second = f_n_facilities;
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

 if( ( AR & FormMsk ) == KskForm ) {  // knapsack formulation- - - - - - - - -
  for( Index i = rng.first ; i < rng.second ; )
   BKB( v_Block[ i++ ] )->set_x( f_n_customers , *(Sol++) );
  return;
  }

 // all other formulations- - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index i = rng.first ; i < rng.second ; )
  v_y[ i++ ].set_value( *(Sol++) );

 }  // end( set_y( Range ) )
 
/*--------------------------------------------------------------------------*/

template< typename T >
void CapacitatedFacilityLocationBlock::set_y(
	       typename std::vector< T >::const_iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  const std::string _prfx = "set_facility_solution: ";
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std::any_of( nms.begin() , nms.end() ,
		  [ & ]( auto i ) { return( i > f_n_facilities ); } ) )
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
		  typename std::vector< T >::const_iterator Sol , Range rng )
{
 #ifndef NDEBUG
  if( ! ( AR & HasVar ) )
   throw( std::logic_error(
		  "set_transportation_solution: variables not generated" ) );
 #endif

 if( rng.second > f_n_facilities * f_n_customers )
  rng.second = f_n_facilities * f_n_customers;
 if( rng.second <= rng.first )  // Range is empty
  return;                       // nothing to do

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
  // it's at least two facilities: deal with the first
  BKB( v_Block[ i++ ] )->set_x( Sol , Range( frst , f_n_customers ) );
  Index k = f_n_customers - frst;
  Index h = rng.first + k;
  Sol += k;
  // now from the second on
  while( h < rng.second ) {
   k = std::min( rng.second - h , f_n_customers );
   BKB( v_Block[ i++ ] )->set_x( Sol , Range( 0 , k ) );
   if( k < f_n_customers )
    break;
   h += f_n_customers;
   Sol += f_n_customers;
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
	       typename std::vector< T >::const_iterator Sol , c_Subset nms )
{
 if( nms.empty() )  // Subset is empty
  return;           // nothing to do

 #ifndef NDEBUG
  auto nmx = f_n_facilities * f_n_customers;
  const std::string _prfx = "set_transportation_solution: ";
  if( ! ( AR & HasVar ) )
   throw( std::logic_error( _prfx + "variables not generated" ) );
  if( std::any_of( nms.begin() , nms.end() ,
		   [ & ]( auto i ) { return( i > nmx ); } ) )
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
 MCFB( v_Block[ 1 ] )->set_x( SSol.begin() , nnms );

 }  // end( set_x( Subset ) )
 
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
   throw( std::invalid_argument( _prfx + "NFacilities dimension required"
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
   throw( std::invalid_argument( _prfx +
				 "TransportationDim dimension required" ) );
  v_x.resize( td.getSize() );
  ts.getVar( v_x.data() );
  }
 }  // end( CapacitatedFacilityLocationSolution::deserialize )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::read( const Block * block )
{
 auto CFLB = dynamic_cast< const CapacitatedFacilityLocationBlock * >( block );
 if( ! CFLB )
  throw( std::invalid_argument(
		        "block is not a CapacitatedFacilityLocationBlock" ) );

 // read y- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_y.empty() ) {
  if( v_x.size() < CFLB->get_NFacilities() )
   v_x.resize( CFLB->get_NFacilities() );

  CFLB->get_facility_solution( v_y.begin() );
  }

 // read x- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < CFLB->get_NFacilities() * CFLB->get_NCustomers() )
   v_x.resize( CFLB->get_NFacilities() * CFLB->get_NCustomers() );

  CFLB->get_transportation_solution( v_x.begin() );
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

  CFLB->set_transportation_solution( v_x.begin() );
  }
 }  // end( CapacitatedFacilityLocationSolution::write )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationSolution::serialize( netCDF::NcGroup & group
						     ) const
{
 if( ! v_y.empty() ) {
  auto nf = group.addDim( "NFacilities" , v_y.size() );
  ::serialize( group , "FacilitySolution" , netCDF::NcDouble() , nf , v_y );
  }

 if( ! v_x.empty() ) {
  auto td = group.addDim( "TransportationDim" , v_x.size() );
  ::serialize( group , "TransportationSolution" , netCDF::NcDouble() , td ,
	       v_x );
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
   sol->v_y.resize( v_y.size() );

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
