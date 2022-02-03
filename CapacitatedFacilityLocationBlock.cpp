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

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using FNumber = CapacitatedFacilityLocationBlock::FNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
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

void CapacitatedFacilityLocationBlock::load( Index n , Index m , c_Subset & pEn , c_Subset & pSn ,
		     c_Vec_FNumber & pU , c_Vec_CNumber & pC ,
		     c_Vec_FNumber & pB , Index dn , Index dm ,
		     Index mdn , Index mdm )
{
 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( pSn.size() < m )
  throw( std::invalid_argument( "pSn too small" ) );

 if( pEn.size() < m )
  throw( std::invalid_argument( "pEn too small" ) );

 if( ( pC.size() > 0 ) && ( pC.size() < m ) )
  throw( std::invalid_argument( "pC nonempty but too small" ) );

 if( ( pU.size() > 0 ) && ( pU.size() < m ) )
  throw( std::invalid_argument( "pU nonempty but too small" ) );

 if( ( pB.size() > 0 ) && ( pB.size() < n ) )
  throw( std::invalid_argument( "pB nonempty but too small" ) );

 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // copy over problem data - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNodes = n;
 NArcs = m;
 MaxNNodes = NNodes + ( mdn > dn ? mdn - dn : 0 );
 c_Index MaxNArcs = NArcs + ( mdm > dm ? mdm - dm : 0 );
 NStaticNodes = dn > n ? 0 : n - dn;
 NStaticArcs = dm > NArcs ? 0 : NArcs - dm;

 SN.resize( MaxNArcs );
 if( std::any_of( pSn.begin() , pSn.begin() + m ,
		  [ n ]( c_Index sn ) { return( ( sn < 1 ) || ( sn > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong starting node" ) );
 std::copy( pSn.begin() , pSn.begin() + m , SN.begin() );

 EN.resize( MaxNArcs );
 if( std::any_of( pEn.begin() , pEn.begin() + m ,
		  [ n ]( c_Index en ) { return( ( en < 1 ) || ( en > n ) ); }
		  ) )
  throw( std::invalid_argument( "wrong ending node" ) );
 std::copy( pEn.begin() , pEn.begin() + m , EN.begin() );

 if( std::any_of( pC.begin() , pC.begin() + m ,
		  []( c_CNumber ci ) { return( ci != 0 ); } ) ) {
  C.resize( MaxNArcs );
  std::copy( pC.begin() , pC.begin() + m , C.begin() );
  }
 else
  C.clear();

 if( std::any_of( pU.begin() , pU.begin() + m ,
		  []( c_FNumber ui ) { return( ui < Inf<FNumber>() ); } ) ) {
  U.resize( MaxNArcs );
  std::copy( pU.begin() , pU.begin() + m , U.begin() );
  }
 else
  U.clear();

 if( std::any_of( pB.begin() , pB.begin() + n ,
		  []( c_FNumber bi ) { return( bi != 0 ); } ) ) {
  B.resize( MaxNNodes );
  std::copy( pB.begin() , pB.begin() + n , B.begin() );
  }
 else
  B.clear();

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // throw Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 // the arc whose cost is infinite has to be closed,
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0; j < C.size() ; j++ )
  if( C[j] >= Inf<double>() ) {
   close_arc(j,eNoMod);
   C[j] = 0;
   }

 }  // end( CapacitatedFacilityLocationBlock::load( memory ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::load( std::istream &input )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();

 // read first non-comment line - - - - - - - - - - - - - - - - - - - - - - -

 char c;
 if( ! ( input >> eatDMXcomments >> c ) )
  throw( std::invalid_argument( "error reading the input stream" ) );

 if( c != 'p' )
  throw( std::invalid_argument( "format error in the input stream" ) );

 input >> eatDMXcomments;
 input.ignore( 3 , ' ' );  // skip "min"

 if( ! ( input >> eatDMXcomments >> NNodes ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of nodes" ) );

 Index tm;
 if( ! ( input >> eatDMXcomments >> NArcs ) )
  throw( std::invalid_argument( "LoadDMX: error reading number of arcs" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SN.resize( NArcs );
 EN.resize( NArcs );
 C.assign( NArcs , 0 );
 U.assign( NArcs , Inf<FNumber>() );
 B.assign( NNodes , 0 );

 NStaticNodes = MaxNNodes = NNodes;
 NStaticArcs = NArcs;

 // read problem data - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i = 0;  // arc counter
 for(;;) {
  if( ! ( input >> eatDMXcomments >> c ) )  // read next descriptor
   break;                                   // if none, end

  switch( c ) {
   case( 'n' ):  // description of a node
    Index j;
    if( ! ( input >> j ) )
     throw( std::invalid_argument( "error reading node name" ) );

    if( ( j < 1 ) || ( j > NNodes ) )
     throw( std::invalid_argument( "invalid node name" ) );

    FNumber Dfctj;
    if( ! ( input >> Dfctj ) )
     throw( std::invalid_argument( "error reading deficit" ) );

    B[ j - 1 ] -= Dfctj;
    break;

   case( 'a' ):  // description of an arc
    if( i == NArcs )
     throw( std::invalid_argument( "too many arc descriptors" ) );

    if( ! ( input >> SN[ i ] ) )
     throw( std::invalid_argument( "error reading start node" ) );

    if( ( SN[ i ] < 1 ) || ( SN[ i ] > NNodes ) )
     throw( std::invalid_argument( "invalid start node" ) );

    if( ! ( input >> EN[ i ] ) )
     throw( std::invalid_argument( "error reading end node" ) );

    if( ( EN[ i ] < 1 ) || ( EN[ i ] > NNodes ) )
     throw( std::invalid_argument( "LoadDMX: invalid end node" ) );

    if( SN[ i ] == EN[ i ] )
     throw( std::invalid_argument( "self-loops not permitted" ) );

    FNumber LB;
    if( ! ( input >> LB ) )
     throw( std::invalid_argument( "error reading lower bound" ) );

    U[ i ] = read_UB( input );

    if( ! ( input >> C[ i ] ) )
     throw( std::invalid_argument( "error reading arc cost" ) );

    if( U[ i ] < LB )
     throw( std::invalid_argument( "lower bound > upper bound" ) );

    if( LB > 0 ) {
     if( U[ i ] < Inf<CapacitatedFacilityLocationBlock::FNumber>() )
      U[ i ] -= LB;
     B[ SN[ i ] - 1 ] += LB;
     B[ EN[ i ] - 1 ] -= LB;
     }
    i++;
    break; 

   default:  // invalid code- - - - - - - - - - - - - - - - - - - - - - - - -
    throw( std::invalid_argument( "invalid DMX code" ) );

   }  // end( switch( c ) )
  }  // end( for( ever ) )

 if( i < NArcs )
  throw( std::invalid_argument( "too few arc descriptors" ) );

 f_cond_lower = NAN;  // reset conditional bounds

 // simplify out the deta structures- - - - - - - - - - - - - - - - - - - - -

 if( std::all_of( C.begin() , C.end() ,
		  []( c_CNumber ci ) { return( ci == 0 ); } ) )
  C.clear();

 if( std::all_of( B.begin() , B.end() ,
		  []( c_FNumber bi ) { return( bi == 0 ); } ) )
  B.clear();

 if( std::all_of( U.begin() , U.end() ,
		  []( c_FNumber ui ) { return( ui == Inf<FNumber>() ); } ) )
  U.clear();
 
 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // issue Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this is a NBModification, the "nuclear option"

 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );

 // the arc whose cost is infinite has to be closed,
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0; j < C.size() ; j++ )
  if( C[j] >= Inf<double>() ) {
   close_arc(j,eNoMod);
   C[j] = 0;
   }

 }  // end( CapacitatedFacilityLocationBlock::load( istream ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any- - - - - - - - - - - - - - - - - - - - - -

 if( MaxNNodes || get_MaxNArcs() )
  guts_of_destructor();
		   
 // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.getDim( "NNodes" );
 if( nn.isNull() )
  throw( std::logic_error( "NNodes dimension is required" ) );
 NNodes = nn.getSize();

 netCDF::NcDim na = group.getDim( "NArcs" );
 if( na.isNull() )
  throw( std::logic_error( "NArcs dimension is required" ) );
 NArcs = na.getSize();

 Index DynNNodes = 0;
 netCDF::NcDim dn = group.getDim( "DynNNodes" );
 if( dn.isNull() )
  NStaticNodes = NNodes;
 else {
  DynNNodes = dn.getSize();
  NStaticNodes = DynNNodes > NNodes ? 0 : NNodes - DynNNodes;
  }

 Index DynNArcs = 0;
 netCDF::NcDim dm = group.getDim( "DynNArcs" );
 if( dm.isNull() )
  NStaticArcs = NArcs;
 else {
  DynNArcs = dm.getSize();
  NStaticArcs = DynNArcs > NArcs ? 0 : NArcs - DynNArcs;
  }

 MaxNNodes = NNodes;
 netCDF::NcDim mdn = group.getDim( "MaxDynNNodes" );
 if( ( ! mdn.isNull() ) && ( mdn.getSize() > DynNNodes ) )
  MaxNNodes += mdn.getSize() - DynNNodes;

 Index MaxNArcs = NArcs;
 netCDF::NcDim mdm = group.getDim( "MaxDynNArcs" );
 if( ( ! mdm.isNull() ) && ( mdm.getSize() > DynNArcs ) )
  MaxNArcs += mdm.getSize() - DynNArcs;
 
 netCDF::NcVar sn = group.getVar( "SN" );
 if( sn.isNull() )
  throw( std::logic_error( "Starting Nodes not found" ) );

 SN.resize( MaxNArcs );

 sn.getVar( SN.data() );

 netCDF::NcVar en = group.getVar( "EN" );
 if( en.isNull() )
  throw( std::logic_error( "Ending Nodes not found" ) );

 EN.resize( MaxNArcs );
 en.getVar( EN.data() );

 netCDF::NcVar cst = group.getVar( "C" );
 if( ! cst.isNull() ) {
  C.resize( MaxNArcs );
  cst.getVar( C.data() );
  if( std::all_of( C.begin() , C.begin() + NArcs ,
		   []( c_CNumber ci ) { return( ci == 0 ); } ) )
   C.clear();
  }

 netCDF::NcVar cap = group.getVar( "U" );
 if( ! cap.isNull() ) {
  U.resize( MaxNArcs );
  cap.getVar( U.data() );
  if( std::all_of( U.begin() , U.begin() + NArcs ,
		   []( c_FNumber ui ) { return( ui == Inf<FNumber>() ); } ) )
   U.clear();
  }

 netCDF::NcVar dfc = group.getVar( "B" );
 if( ! dfc.isNull() ) {
  B.resize( MaxNNodes );
  std::vector<size_t> countn = { NNodes };
  dfc.getVar( B.data() );
  if( std::all_of( B.begin() , B.begin() + NNodes ,
		   []( c_FNumber bi ) { return( bi == 0 ); } ) )
   B.clear();
  }

 f_cond_lower = NAN;  // reset conditional bounds

 // allocate flow variables - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_abstract_variables();

 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -
 // inside this the NBModification, the "nuclear option",  is issued

 Block::deserialize( group );

 // the arc whose cost is infinite has to be closed,
 // in addition the cost has to be set to 0 - - - - - - - - - - - - - - - - -

 for( Index j = 0; j < C.size() ; j++ )
  if( C[j] >= Inf<double>() ) {
   close_arc(j,eNoMod);
   C[j] = 0;
   }

 }  // end( CapacitatedFacilityLocationBlock::deserialize )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_abstract_variables( Configuration *stvv )
{
 if( AR & HasVar )  // the variables are there already
  return;           // nothing to do

 if( HasStaticX() ) {
  x.resize( get_NStaticArcs() );
  for( auto & var : x )
   var.is_positive( true , eNoBlck );

  add_static_variable( x );
  }

 if( MayHaveDynX() ) {
  dx.resize( get_NArcs() - get_NStaticArcs() );
  for( auto & var : dx )
   var.is_positive( true , eNoBlck );

  add_dynamic_variable( dx );
  } 

 AR |= HasVar;

 }  // end( CapacitatedFacilityLocationBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_abstract_constraints( Configuration *stcc )
{
 if( ! ( AR & HasFlw ) ) {
  // count number of nonzeroes in each constraint, i.e., #FS( i ) + #BS( i )
  Subset count( get_NNodes() );
 
  for( Index i = 0 ; i < get_NStaticArcs() ; ++i ) {
   count[ SN[ i ] - 1 ]++;
   count[ EN[ i ] - 1 ]++;
   }

  for( Index i = NStaticArcs ; i < get_NArcs() ; ++i )
   if( ! is_deleted( i ) ) {
    count[ SN[ i ] - 1 ]++;
    count[ EN[ i ] - 1 ]++;
    }

  // initialize the vectors of coefficients, and reset count[]
  std::vector< LinearFunction::v_coeff_pair > coeffs( get_NNodes() );

  for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   coeffs[ i ].resize( count[ i ] );
   count[ i ] = 0;
   }

  // construct the vector of coefficients, static phase
  Index i = 0;
  for( ; i < get_NStaticArcs() ; ++i ) {
   coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( -1 ) );
   coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &x[ i ] , double( 1 ) );
   }

  // construct the vector of coefficients, dynamic phase
  if( MayHaveDynX() )
   for( auto dxi = dx.begin() ; i < get_NArcs() ; ++i , ++dxi )
    if( ! is_deleted( i ) ) {
     coeffs[ SN[ i ] - 1 ][ count[ SN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( -1 ) );
     coeffs[ EN[ i ] - 1 ][ count[ EN[ i ] - 1 ]++ ] =
                                    std::make_pair( &(*dxi) , double( 1 ) );
     }

  // generate the node-arc incidence matrix - - - - - - - - - - - - - - - - -
  // each constraint is an equality, i.e., LHS = RHS = B[ i ]

  // static part
  if( HasStaticE() ) {
   E.resize( get_NStaticNodes() );

   for( Index i = 0 ; i < get_NStaticNodes() ; ++i ) {
    E[ i ].set_both( B.empty() ? 0 : B[ i ] );
    E[ i ].set_function( new LinearFunction( std::move( coeffs[ i ] ) , 0 ) );
    }

   add_static_constraint( E );
   }

  // dynamic part
  if( MayHaveDynE() ) {
   dE.resize( get_NNodes() - get_NStaticNodes() );

   Index i = get_NStaticNodes();
   for( auto & cnst : dE ) {
    cnst.set_both( B.empty() ? 0 : B[ i ] );
    cnst.set_function( new LinearFunction( std::move( coeffs[ i++ ] ) , 0 ) );
    }

   add_dynamic_constraint( dE );
   }

  AR |= HasFlw;
  }

 // generate the bound constraints- - - - - - - - - - - - - - - - - - - - - -

 if( AR & HasBnd )  // bound constraints there already
  return;           // nothing to do

 if( U.empty() ) {
  // if upper bounds are not there and the Configuration says so, the
  // LB0Constraintare not constructed

  auto tstcc = dynamic_cast<SimpleConfiguration<int> *>( stcc );

  if( ( ! tstcc ) && f_BlockConfig &&
      f_BlockConfig->f_static_constraints_Configuration )
   tstcc = dynamic_cast<SimpleConfiguration<int> *>(
                         f_BlockConfig->f_static_constraints_Configuration );
  if( tstcc && ( tstcc->f_value != 0 ) )
   return;
  }

 // static part
 if( HasStaticX() ) {
  UB.resize( get_NStaticArcs() );
  for( Index i = 0 ; i < get_NStaticArcs() ; ++i ) {
   UB[ i ].set_variable( & x[ i ] , eNoBlck );
   UB[ i ].set_rhs( U[ i ] , eNoBlck );
   }

  add_static_constraint( UB );
  }

 // dynamic part
 if( MayHaveDynX() ) {
  dUB.resize( get_NArcs() - get_NStaticArcs() );

  auto dxi = dx.begin();
  auto ui = U.begin() + get_NStaticArcs();
  for( auto & cnst : dUB ) {
   cnst.set_variable( &(*(dxi++)) , eNoBlck );
   cnst.set_rhs( *(ui++) , eNoBlck );
   }

  add_dynamic_constraint( dUB );
  }

 AR |= HasBnd;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::generate_objective( Configuration *objc )
{
 if( AR & HasObj )  // the objective is there already
  return;           // cowardly (and silently) return

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 LinearFunction::v_coeff_pair p( get_NArcs() );

 // construct a "dense" LinearFunction- - - - - - - - - - - - - - - - - - - -

 Index i = 0;

 // static part
 if( HasStaticX() ) {
  if( C.empty() )
   for( ; i < get_NStaticArcs() ; ++i ) {
    p[ i ].first = &x[ i ];
    p[ i ].second = 0;
    }
  else
   for( ; i < get_NStaticArcs() ; ++i ) {
    p[ i ].first = &x[ i ];
    p[ i ].second = C[ i ];
    }
  }

 // dynamic part
 if( HasDynamicX() ) {
  auto dxi = dx.begin();
  if( C.empty() )
   for( ; i < get_NArcs() ; ++i ) {
    p[ i ].first = &(*(dxi++));
    p[ i ].second = 0;
    }
  else
   for( ; i < get_NArcs() ; ++i ) {
    p[ i ].first = &(*(dxi++));
    p[ i ].second = C[ i ];
    }
  }

 // ensure no Modification is issued: this may happen in case a CapacitatedFacilityLocationBlock
 // is re-loaded, so that set_objective( c ) had already been called
 c.set_function( new LinearFunction( std::move( p ) , 0 ) , eNoMod );
 c.set_Block( this );

 set_objective( & c , eNoMod );

 AR |= HasObj;

 #if CHECK_DS
  CheckAbsVSPhys();
 #endif

 }  // end( CapacitatedFacilityLocationBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::flow_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ! ( AR & HasFlw ) )
   throw( std::logic_error( "Constraint required for flow_feasible( , true )"
			    ) );
  // static part
  for( const auto & cnst : E )
   if( cnst.rel_viol() > feps )
    return( false );

  // dynamic part
  for( const auto & cnst : dE )
   if( cnst.rel_viol() > feps )
    return( false );
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Index i = 0;
  Vec_FNumber tB = B;

  // static part
  for( ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber xi = x[ i ].get_value();
   tB[ SN[ i ] - 1 ] += xi;
   tB[ EN[ i ] - 1 ] -= xi;
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for( ; i < get_NArcs() ; ++i , ++dxi )
    if( ! is_deleted( i ) ) {
     c_FNumber xi = dxi->get_value();
     tB[ SN[ i ] - 1 ] += xi;
     tB[ EN[ i ] - 1 ] -= xi;
     }
   }

 for( Index i = 0 ; i < get_NNodes() ; ++i ) {
   c_FNumber slck = B[ i ] == 0 ? std::abs( tB[ i ] ) :
                                  std::abs( tB[ i ] / B[ i ] );
   if( slck > feps )
    return( false );
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::flow_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::bound_feasible( c_FNumber feps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ! ( AR & ( HasFlw | HasVar ) ) ) )
   throw( std::logic_error(
	 "abstract representation not there in bound_feasible( , true )" ) );

  // static part
  if( HasStaticX() ) {
   if( UB.empty() ) {
    for( const auto & var : x )
     if( ! var.is_feasible( feps ) )
      return( false );
    }
   else
    for( const auto & cnst : UB )
     if( cnst.rel_viol() > feps )
      return( false );
   }

  // dynamic part
  if( HasDynamicX() ) {
   if( dUB.empty() ) {
    for( const auto & var : dx )
     if( ! var.is_feasible( feps ) )
      return( false );
    }
   else
    for( const auto & cnst : dUB )
     if( cnst.rel_viol() > feps )
      return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for( ; i < get_NStaticArcs() ; ++i ) {
   c_FNumber Ui = get_U( i );
   c_FNumber xi = x[ i ].get_value();
   if( Ui >= Inf<FNumber>() ) {
    if( xi < - feps )
     return( false );
    }
   else {
    c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                               std::max( - xi , xi - Ui ) / std::abs( Ui );
    if( slck > feps )
     return( false );
    }
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for(  ; i < get_NArcs() ; ++i ) {
    c_FNumber Ui = get_U( i );
    c_FNumber xi = (*(dxi++)).get_value();
    if( Ui >= Inf<FNumber>() ) {
     if( xi < - feps )
      return( false );
     }
    else {
     c_FNumber slck = Ui == 0 ? std::abs( xi ) :
                                std::max( - xi , xi - Ui ) / std::abs( Ui );
     if( slck > feps )
      return( false );
     }
    }
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::bound_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::dual_feasible( c_CNumber ceps , bool useabstract )
{
 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ! ( AR & HasFlw ) ) || ( ! ( AR & HasObj ) ) )
   throw( std::logic_error(
	 "abstract representation not there in dual_feasible( , true )" ) );

  auto obj = static_cast<FRealObjective *>( get_objective() );
  assert( obj );
  auto lfo = get_lfo();

  for( auto & pi : (*lfo).get_v_var() ) {
   auto xi = pi.first;
   auto RCi = pi.second;
   for( Index j = 0 ; j < xi->get_num_active() ; ++j ) {
    ThinVarDepInterface * ci = xi->get_active( j );
    if( auto rci = dynamic_cast< FRowConstraint * >( ci ) ) {
     auto lrci = get_lfc( rci );
     RCi -= rci->get_dual() * lrci->get_coefficient( lrci->is_active( xi ) );
     }
    else {
     auto bci = dynamic_cast< BoxConstraint * >( ci );
     assert( bci );
     RCi -= bci->get_dual();
     }
    }

   RCi = std::abs( RCi );
   if( pi.second != 0 )
    RCi /= pi.second;

   if( RCi > ceps )
    return( false );
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -

  Vec_CNumber RC;
  get_rc( RC );
  Vec_CNumber Pi;
  get_pi( Pi );

  for( Index i = 0 ; i < get_NArcs() ; ++i ) {
   if( is_deleted( i ) )
    continue;

   c_CNumber Ci = get_C( i );
   c_CNumber RCi = Ci + Pi[ SN[ i ] - 1 ] - Pi[ EN[ i ] - 1 ];
   c_CNumber df = std::abs( RCi - RC[ i ] );
   c_CNumber mx = std::max( std::abs( Ci ) , df );
   if( mx == 0 ) {
    if( df > ceps )
     return( false );
    }
   else
    if( df > ceps * mx )
     return( false );
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::dual_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::complementary_slackness( c_CNumber ceps , c_FNumber feps ,
					bool useabstract )
{
 Vec_CNumber RC;
 get_rc( RC );

 if( useabstract ) {
  // do it using the abstract representation- - - - - - - - - - - - - - - - -

  if( ( ! ( AR & HasVar ) ) || ( ! ( AR & HasObj ) ) )
   throw( std::logic_error(
    "abstract representation not there in complementary_slackness(( , true )"
			   ) );

  auto obj = static_cast<FRealObjective *>( get_objective() );
  assert( obj );
  auto lfo = get_lfo();
  Index i = 0;

  // static part
  if( HasStaticX() ) {
   if( UB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i ) {
     c_CNumber Ci = lfo->get_coefficient( i );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     if( ( x[ i ].get_value() > feps ) && ( RCi < - ceps ) )
      return( false );
     }
    }
   else
    for( ; i < get_NStaticArcs() ; ++i ) {
     c_CNumber Ci = lfo->get_coefficient( i );
     CNumber RCi = RC[ i ];
     if( Ci )
      RCi /= Ci;
     c_FNumber xiv = x[ i ].get_value();
     c_FNumber UBi = UB[ i ].get_rhs();
     if( UBi >= Inf<RowConstraint::RHSValue>() ) {
      if( ( xiv > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     else {
      c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
      if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	  ( ( UBi - xiv > sfeps ) && ( RCi > ceps ) ) )
       return( false );
      }
     }
    }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   if( dUB.empty() ) {
    for( ; i < get_NStaticArcs() ; ++i )
     if( ! is_deleted( i ) ) {
      c_CNumber Ci = lfo->get_coefficient( i );
      CNumber RCi = RC[ i ];
      if( Ci )
       RCi /= Ci;
      if( ( (dxi++)->get_value() > feps ) && ( RCi < - ceps ) )
       return( false );
      }
    }
   else {
    auto dubi = dUB.begin();

    for( ; i < get_NArcs() ; ++i )
     if( ! is_deleted( i ) ) {
      c_CNumber Ci = lfo->get_coefficient( i );
      CNumber RCi = RC[ i ];
      if( Ci )
       RCi /= Ci;
      c_FNumber dxiv = (dxi++)->get_value();
      c_FNumber UBi = (dubi++)->get_rhs();
      if( UBi >= Inf<RowConstraint::RHSValue>() ) {
       if( ( dxiv > feps ) && ( RCi < - ceps ) )
	return( false );
       }
      else {
       c_FNumber sfeps = ( UBi == 0 ? feps : feps * UBi );
       if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	   ( ( UBi - dxiv > sfeps ) && ( RCi > ceps ) ) )
	return( false );
       }
      }
    }
   }
  }
 else {
  // do it using the physical representation- - - - - - - - - - - - - - - - -
  Index i = 0;

  // static part
  for(  ; i < get_NStaticArcs() ; ++i ) {
   c_CNumber Ci = get_C( i );
   CNumber RCi = RC[ i ];
   if( Ci != 0 )
    RCi /= C[ i ];

   c_FNumber Ui = get_U( i );
   c_FNumber xiv = x[ i ].get_value();

   if( Ui >= Inf<FNumber>() ) {
    if( ( xiv > feps ) && ( RCi < - ceps ) )
     return( false );
    }
   else {
    c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
    if( ( ( xiv > sfeps ) && ( RCi < - ceps ) ) ||
	( ( Ui - xiv > sfeps ) && ( RCi > ceps ) ) )
     return( false );
    }
   }

  // dynamic part
  if( HasDynamicX() ) {
   auto dxi = dx.begin();

   for( ; i < get_NArcs() ; ++i , ++dxi )
    if( ! is_deleted( i ) ) {
     c_FNumber dxiv = dxi->get_value();
     c_CNumber Ci = get_C( i );
     CNumber RCi = RC[ i ];
     if( Ci != 0 )
      RCi /= C[ i ];

     c_FNumber Ui = get_U( i );
     if( Ui >= Inf<FNumber>() ) {
      if( ( dxiv > feps ) && ( RCi < - ceps ) )
       return( false );
      }
     else {
      c_FNumber sfeps = ( Ui == 0 ? feps : feps * Ui );
      if( ( ( dxiv > sfeps ) && ( RCi < - ceps ) ) ||
	  ( ( Ui - dxiv > sfeps ) && ( RCi > ceps ) ) )
       return( false );
      }
     }
   }
  }

 return( true );

 }  // end( CapacitatedFacilityLocationBlock::complementary_slackness )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::is_feasible( bool useabstract , Configuration *fsbc )
{
 FNumber eps = 0;
 auto tfsbc = dynamic_cast<SimpleConfiguration<FNumber> *>( fsbc );

 if( ( ! tfsbc ) && f_BlockConfig &&
     f_BlockConfig->f_is_feasible_Configuration )
  tfsbc = dynamic_cast<SimpleConfiguration<FNumber> *>(
                        f_BlockConfig->f_is_feasible_Configuration );
 if( tfsbc )
  eps = tfsbc->f_value;

 return( flow_feasible( eps , useabstract ) &&
	 bound_feasible( eps , useabstract ) );

 }  //  end( CapacitatedFacilityLocationBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::is_optimal( bool useabstract , Configuration *optc )
{
 CNumber ceps = 0;
 FNumber feps = 0;
 if( optc ) {
  if( auto toptc =
      dynamic_cast< SimpleConfiguration< std::pair< CNumber , FNumber > > * >(
								    optc ) ) {
   ceps = toptc->f_value.first;
   feps = toptc->f_value.second;
   }
  else {
   auto ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >( optc );

   if( ( ! ttoptc ) && f_BlockConfig &&
       f_BlockConfig->f_is_optimal_Configuration )
    ttoptc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                           f_BlockConfig->f_is_optimal_Configuration );
   if( ttoptc )
    ceps = ttoptc->f_value;

   if( f_BlockConfig && f_BlockConfig->f_is_feasible_Configuration ) {
    auto fsbc = dynamic_cast< SimpleConfiguration< FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration );
    if( fsbc )
     feps = fsbc->f_value;
    }
   }
  }
 else
  if( f_BlockConfig ) {
   if( f_BlockConfig->f_is_optimal_Configuration )
    if( auto csbc = dynamic_cast< SimpleConfiguration< CNumber > * >(
                              f_BlockConfig->f_is_optimal_Configuration ) )
     ceps = csbc->f_value;

   if( f_BlockConfig->f_is_feasible_Configuration )
    if( auto fsbc = dynamic_cast< SimpleConfiguration<FNumber > * >(
                              f_BlockConfig->f_is_feasible_Configuration ) )
     feps = fsbc->f_value;
   }

 return( flow_feasible( feps , useabstract ) &&
	 bound_feasible( feps , useabstract ) &&
	 dual_feasible( ceps , useabstract ) &&
	 complementary_slackness( ceps , feps , useabstract ) );

 }  //  end( CapacitatedFacilityLocationBlock::is_optimal )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/

Block * CapacitatedFacilityLocationBlock::get_R3_Block( Configuration *r3bc , Block * base  ,
				Block * father )
{
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 CapacitatedFacilityLocationBlock *MCFB;
 if( base ) {
  MCFB = dynamic_cast< CapacitatedFacilityLocationBlock * >( base );
  if( ! MCFB )
   throw( std::invalid_argument( "base is not a CapacitatedFacilityLocationBlock" ) );
  }
 else
  MCFB = new CapacitatedFacilityLocationBlock( father );

 MCFB->load( get_NNodes() , get_NArcs() , EN , SN , U , C , B ,
	     get_NNodes() - get_NStaticNodes() ,
	     get_NArcs() - get_NStaticArcs() ,
	     get_MaxNNodes() - get_NStaticNodes() ,
	     get_MaxNArcs() - get_NStaticArcs() );
 
 return( MCFB );

 }  // end( CapacitatedFacilityLocationBlock::get_R3_Block )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::map_back_solution( Block *R3B , Configuration *r3bc ,
				               Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast<CapacitatedFacilityLocationBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a CapacitatedFacilityLocationBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
   tsolc = dynamic_cast<SimpleConfiguration<int> *>(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map back primal solution - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 2 ) && ( AR & HasVar ) ) {  // ... if any
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! xi->is_fixed() )
     xi->set_value( r3bxi->get_value() );
 
  // dynamic part
  // note that if MCFB->dx is longer than this->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicX() ) {
   auto dxi = dx.begin();
   for( auto r3bdxi = MCFB->dx.begin() ;
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! dxi->is_fixed() )
     dxi->set_value( r3bdxi->get_value() );

   for( ; dxi != dx.end() ; ++dxi )
    if( ! dxi->is_fixed() )
     dxi->set_value();
   }
  }

 // if required, map back dual solution - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any

  // map back the potentials- - - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (ei++)->set_dual( (r3bei++)->get_dual() );
 
  // dynamic part
  // note that if MCFB->dE is longer than this->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( HasDynamicE() ) {
   auto dei = dE.begin();

   for( auto r3bdei = MCFB->dE.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (dei++)->set_dual( (r3bdei++)->get_dual() );
 
   while( dei != dE.end() )
    (dei++)->set_dual( 0 );
   }

  // map back the reduced costs - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( HasStaticX() && ( ! UB.empty() ) ) {
   if( MCFB->UB.empty() ) {
    for( auto & cnst : UB )
     cnst.set_dual();
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (drci++)->set_dual( (r3bdrci++)->get_dual() );
   }

  // dynamic part
  if( HasDynamicX() && ( ! dUB.empty() ) ) {
   if( MCFB->dUB.empty() ) {
    for( auto & cnst : dUB )
     cnst.set_dual();
    }
   else {
    auto drci = dUB.begin();
    for( auto r3bdrci = MCFB->dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ;
	 ++drci , ++r3bdrci )
     drci->set_dual( r3bdrci->get_dual() );
 
    while( drci != dUB.end() )
     (drci++)->set_dual();
    }
   }
  }
 }  // end( CapacitatedFacilityLocationBlock::map_back_solution )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::map_forward_solution( Block *R3B , Configuration *r3bc ,
				                  Configuration *solc )
{
 // process Configuration - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto MCFB = dynamic_cast<CapacitatedFacilityLocationBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a CapacitatedFacilityLocationBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast<SimpleConfiguration<int> *>(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 // if required, map forward primal solution- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any
  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static flow size" ) );

  // static part
  if( MCFB->HasStaticX() )
   for( auto xi = x.begin() , r3bxi = MCFB->x.begin() ; xi != x.end() ;
	++xi , ++r3bxi )
    if( ! r3bxi->is_fixed() )
     r3bxi->set_value( xi->get_value() );
 
  // dynamic part
  // note that if this->dx is longer than MCFB->dx the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicX() ) {
   auto r3bdxi = MCFB->dx.begin();
   for( auto dxi = dx.begin();
	( dxi != dx.end() ) && ( r3bdxi != MCFB->dx.end() ) ;
	++dxi , ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value( dxi->get_value() );

   for( ; r3bdxi != MCFB->dx.end() ; ++r3bdxi )
    if( ! r3bdxi->is_fixed() )
     r3bdxi->set_value();
   }
  }

 // if required, map forward dual solution- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wsol != 1 ) && ( AR & HasFlw ) ) {  // ... if any
  // map forward the potentials - - - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticNodes() != get_NStaticNodes() )
   throw( std::invalid_argument( "incompatible static potential size" ) );

  // static part
  if( MCFB->HasStaticE() )
   for( auto ei = E.begin() , r3bei = MCFB->E.begin() ; ei != E.end() ; )
    (r3bei++)->set_dual( (ei++)->get_dual() );
 
  // dynamic part
  // note that if this->dE is longer than MCFB->dE the last part is
  // ignored, while if the converse happens it is filled with zeros
  if( MCFB->HasDynamicE() ) {
   auto r3bdei = MCFB->dE.begin();

   for( auto dei = dE.begin() ;
	( dei != dE.end() ) && ( r3bdei != MCFB->dE.end() ) ; )
    (r3bdei++)->set_dual( (dei++)->get_dual() );
 
   while( r3bdei != MCFB->dE.end() )
    (r3bdei++)->set_dual( 0 );
   }

  // map forward the reduced costs- - - - - - - - - - - - - - - - - - - - - -

  if( MCFB->get_NStaticArcs() != get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible static reduced cost size" ) );

  // static part
  if( MCFB->HasStaticX() && ( ! MCFB->UB.empty() ) ) {
   if( UB.empty() ) {
    for( auto & cnst : MCFB->UB )
     cnst.set_dual( 0 );
    }
   else
    for( auto drci = UB.begin() , r3bdrci = MCFB->UB.begin() ;
	 drci != UB.end() ; )
      (r3bdrci++)->set_dual( (drci++)->get_dual() );
   }

  // dynamic part
  if( MCFB->HasDynamicX() && ( ! MCFB->dUB.empty() ) ) {
   if( UB.empty() ) {
    for( auto & cnst : MCFB->dUB )
     cnst.set_dual( 0 );
    }
   else {
    auto r3bdrci = MCFB->dUB.begin();

    for( auto drci = dUB.begin() ;
	 ( drci != dUB.end() ) && ( r3bdrci != MCFB->dUB.end() ) ; )
     (r3bdrci++)->set_dual( (drci++)->get_dual() );
 
    while( r3bdrci != MCFB->dUB.end() )
     (r3bdrci++)->set_dual( 0 );
    }
   }
  }
 }  // end( CapacitatedFacilityLocationBlock::map_forward_solution )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::map_forward_Modification( Block *R3B , c_p_Mod mod ,
					 Configuration *r3bc ,
					 ModParam issuePMod ,
					 ModParam issueAMod )
{
 if( mod->concerns_Block() )  // an abstract Modification
  return( false );            // none of my business
 
 auto MCFB = dynamic_cast< CapacitatedFacilityLocationBlock * >( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a CapacitatedFacilityLocationBlock" ) );
 if( r3bc != nullptr )
  throw( std::invalid_argument( "non-nullptr R3B Configuration" ) );

 /* When a GroupModification is processed, if no channel is provided, then
    one is opened. This only happens "at root", after which in guts_of_mfM()
    whenever a GroupModification is processed, then the channel is nested.
    Indeed, if the "root" Modification is not a GroupModification, then there
    cannot be any GroupModification in it. */

 ModParam iPM = issuePMod;
 ModParam iPA = make_par( std::min( ModParam( eNoBlck ) ,
				    par2mod( issueAMod ) ) ,
			  par2chnl( issueAMod ) );

 /* Use a Lambda to define a "guts" of the method that can be called
    recursively without having to pass "local globals". Note the trick of
    defining the std::function object and "passing" it to the lambda,
    which allows recursive calls. Note the need to explicitly capture
    "this" to use fields/methods of the class. */

 std::function< bool( c_p_Mod )> guts_of_mfM;
 guts_of_mfM = [ this , & guts_of_mfM , & MCFB , & iPM , & iPA ]( c_p_Mod mod
								  ) {
  // process Modification- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* This requires to patiently sift through the possible Modification types
     to find what this Modification exactly is, and call the appropriate
     method of either MCFB, for a "physical Modification", or of the "abstract
     representation" of MCFB for an "abstract Modification". */

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
     (a Modification changin nothing from the "empty" state is not issued). */

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
   * CapacitatedFacilityLocationBlockSbstMod only provides read-only access to the vector. */

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

 bool ok = true;  // final return value

 if( const auto tmod = dynamic_cast< GroupModification * const >( mod ) ) {
  // if the channels are the default ones, open new ones
  if( ! par2chnl( issuePMod ) )
   iPM = make_par( par2mod( issuePMod ) , MCFB->open_channel() );
  if( ! par2chnl( issueAMod ) )
   iPA = make_par( par2mod( issueAMod ) , MCFB->open_channel() );

  for( const auto & submod : tmod->sub_Modifications() )  // for each sub-Mod
   if( ! guts_of_mfM( submod.get() ) )                    // make the call
    ok = false;

  // now close the opened channels, if any
  if( ! par2chnl( issuePMod ) )
   MCFB->close_channel( par2chnl( iPM ) );
  if( ! par2chnl( issueAMod ) )
   MCFB->close_channel( par2chnl( iPA ) );
  }
 else                             // any other Modification
  ok = guts_of_mfM( mod );        // just make the call

 return( ok );

 }  // end( CapacitatedFacilityLocationBlock::map_forward_Modification )

/*--------------------------------------------------------------------------*/

bool CapacitatedFacilityLocationBlock::map_back_Modification( Block *R3B , c_p_Mod mod ,
				      Configuration *r3bc ,
				      ModParam issuePMod ,
				      ModParam issueAMod )
{
 /* Fantastically dirty trick: because the two objects are copies, mapping
    back a Modification to this from R3B is the same as mapping forward a
    Modification from R3B to this. */

 auto MCFB = dynamic_cast<CapacitatedFacilityLocationBlock *>( R3B );
 if( ! MCFB )
  throw( std::invalid_argument( "R3B is not a CapacitatedFacilityLocationBlock" ) );

 return( MCFB->map_forward_Modification( this , mod , r3bc , issuePMod ,
					 issueAMod ) );

 }  // end( CapacitatedFacilityLocationBlock::map_back_Modification )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * CapacitatedFacilityLocationBlock::get_Solution( Configuration *solc , bool emptys )
{
 int wsol = 0;
 auto tsolc = dynamic_cast<SimpleConfiguration<int> *>( solc );

 if( ( ! tsolc ) && f_BlockConfig && f_BlockConfig->f_solution_Configuration )
  tsolc = dynamic_cast<SimpleConfiguration<int> *>(
                                    f_BlockConfig->f_solution_Configuration );
 if( tsolc )
  wsol = tsolc->f_value;

 auto *sol = new MCFSolution();

 if( wsol != 2 )
  sol->v_x.resize( get_NArcs() );

 if( wsol != 1 )
  sol->v_pi.resize( get_NNodes() );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( CapacitatedFacilityLocationBlock::get_Solution )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_x( Vec_FNumber & FSol , Range rng )
{
 auto FSi = FSol.begin();
 for( ; rng.first < std::min( rng.second , get_NStaticArcs() ) ; )
  *(FSi++) = x[ rng.first++ ].get_value();

 if( HasDynamicX() ) {
  auto dxi = dx.begin();
  for( ; rng.first++ < std::min( rng.second , get_NArcs() ) ; )
   *(FSi++) = (*(dxi++)).get_value();
  }
 }  // end( CapacitatedFacilityLocationBlock::get_x( interval ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_x( Vec_FNumber & FSol , c_Subset & nms )
{
 if( ! ( AR & HasVar ) )
  throw( std::logic_error( "flow Variable not available" ) );

 auto FSi = FSol.begin();
 auto nmsi = nms.begin();

 if( HasDynamicX() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
   *(FSi++) = x[ *(nmsi++) ].get_value();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticArcs();
  for( auto dxi = dx.begin() ; nmsi != nms.end() ; ++dxi )
   if( *nmsi == i++ ) {
    *(FSi++) = (*dxi).get_value();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() ) 
   *(FSi++) = x[ *(nmsi++) ].get_value();

 }  // end( CapacitatedFacilityLocationBlock::get_x( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_pi( Vec_CNumber & PSol , Range rng )
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 auto PSi = PSol.begin();
 Index i = rng.first;
 for( ; i < std::min( rng.second , get_NStaticNodes() ) ; )
  *(PSi++) = E[ i++ ].get_dual();

 if( HasDynamicE() ) {
  auto dei = std::next( dE.begin() , rng.first >= get_NStaticNodes() ?
			             i - get_NStaticNodes() : 0 );
   for( ; i++ < std::min( rng.second , get_NNodes() ) ; )
   *(PSi++) = (*(dei++)).get_dual();
  }
 }  // end( CapacitatedFacilityLocationBlock::get_pi( interval ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_pi( Vec_CNumber & PSol , c_Subset & nms )
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "potentials unavailable if Constraint aren't" ) );

 auto PSi = PSol.begin();
 auto nmsi = nms.begin();

 if( HasDynamicE() ) {
  while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticNodes() ) )
   *(PSi++) = E[ *(nmsi++) ].get_dual();

  if( nmsi == nms.end() )
   return;

  auto i = get_NStaticNodes();
  for( auto dei = dE.begin() ; nmsi != nms.end() ; ++dei )
   if( *nmsi == i++ ) {
    *(PSi++) = (*dei).get_dual();
    nmsi++;
    }
  }
 else
  while( nmsi != nms.end() )
   *(PSi++) = E[ *(nmsi++) ].get_dual();

 }  // end( CapacitatedFacilityLocationBlock::get_pi( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_rc( Vec_CNumber & RC , Range rng )
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 auto RCi = RC.begin();

 if( AR & HasBnd ) {
  Index i = rng.first;
  if( HasStaticX() )
   for( ; i < std::min( rng.second , get_NStaticArcs() ) ; )
    *(RCi++) = UB[ i++ ].get_dual();

  if( HasDynamicX() ) {
   auto dubi = std::next( dUB.begin() , rng.first >= get_NStaticArcs() ?
			                i - get_NStaticArcs() : 0 );
   for( ; i++ < std::min( rng.second , get_NArcs() ) ; )
    *(RCi++) = (*(dubi++)).get_dual();
   }
  }
 else
  for( ; rng.first < std::min( rng.second , get_NArcs() ) ; ++rng.first )
   *(RCi++) = get_C( rng.first ) + get_pi( SN[ rng.first ] - 1 )
                                 - get_pi( EN[ rng.first ] - 1 );

 }  // end( CapacitatedFacilityLocationBlock::get_rc( interval ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::get_rc( Vec_CNumber & RC , c_Subset & nms )
{
 if( ! ( AR & HasFlw ) )
  throw( std::logic_error( "reduced costs unavailable if Constraint aren't" )
	 );

 auto RCi = RC.begin();
 auto nmsi = nms.begin();

 if( AR & HasBnd ) {
  if( HasDynamicX() ) {
   while( ( nmsi != nms.end() ) && ( *nmsi < get_NStaticArcs() ) )
     *(RCi++) = UB[ *(nmsi++) ].get_dual();

   if( nmsi == nms.end() )
    return;

   auto i = get_NStaticArcs();
   for( auto dubi = dUB.begin() ; nmsi != nms.end() ; ++dubi )
    if( *nmsi == i++ ) {
     *(RCi++) = (*dubi).get_dual();
     nmsi++;
     }
   }
  else
   while( nmsi != nms.end() ) 
    *(RCi++) = UB[ *(nmsi++) ].get_dual();
  }
 else
  for( ; nmsi != nms.end() ; ++nmsi ) 
   *(RCi++) = get_C( *nmsi ) + get_pi( SN[ *nmsi ] - 1 )
                             - get_pi( EN[ *nmsi ] - 1 );
 
 }  // end( CapacitatedFacilityLocationBlock::get_rc( subset ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::set_x( c_Vec_FNumber_it fstrt , c_Vec_FNumber_it fstop ,
		      c_Index strt )
{
 if( ! ( AR & HasVar ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( std::distance( fstrt , fstop ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticX() )
  for( auto xi = x.begin() + strt ;
       ( fstrt != fstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (xi++)->set_value( *(fstrt++) );

 if( ( fstrt == fstop ) || ( ! HasDynamicX() ) )
  return;

 auto dxi = dx.begin();
 if( i > get_NStaticArcs() )
  dxi = std::next( dxi , i - get_NStaticArcs() );

 while( fstrt != fstop )
  (dxi++)->set_value( *(fstrt++) );

 }  // end( CapacitatedFacilityLocationBlock::set_x( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::set_pi( c_Vec_CNumber_it pstrt , c_Vec_CNumber_it pstop ,
		       c_Index strt )
{
 if( ! ( AR & HasFlw ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( std::distance( pstrt , pstop ) + strt > get_NNodes() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticE() )
  for( auto ei = E.begin() + strt ;
       ( pstrt != pstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (ei++)->set_dual( *(pstrt++) );

 if( ( pstrt == pstop ) || ( ! HasDynamicE() ) )
  return;

 auto dei = dE.begin();
 if( i > get_NStaticArcs() )
  dei = std::next( dei , i - get_NStaticArcs() );

 while( pstrt != pstop )
  (dei++)->set_dual( *(pstrt++) );

 }  // end( CapacitatedFacilityLocationBlock::set_pi( range ) )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::set_rc( c_Vec_CNumber_it rcstrt , c_Vec_CNumber_it rcstop ,
		       c_Index strt )
{
 if( ! ( AR & HasBnd ) )  // nowhere to put the value in
  return;                 // cowardly (and silently) return

 if( std::distance( rcstrt , rcstop ) + strt > get_NArcs() )
  throw( std::invalid_argument( "too many values provided" ) );

 Index i = strt;

 if( HasStaticX() )
  for( auto ubi = UB.begin() + strt ;
       ( rcstrt != rcstop ) && ( i < get_NStaticArcs() ) ; ++i )
   (ubi++)->set_dual( *(rcstrt++) );

 if( ( rcstrt == rcstop ) || ( ! HasDynamicX() ) )
  return;

 auto dubi = dUB.begin();
 if( i > get_NStaticArcs() )
  dubi = std::next( dubi , i - get_NStaticArcs() );

 while( rcstrt != rcstop )
  (dubi++)->set_dual( *(rcstrt++) );

 }  // end( CapacitatedFacilityLocationBlock::set_rc( range ) )

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 //!! std::cout << *mod << std::endl;

 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod , chnl );
 }

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE CapacitatedFacilityLocationBlock ---------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::serialize( netCDF::NcGroup & group ) const
{
 // call the method of Block- - - - - - - - - - - - - - - - - - - - - - - - -

 Block::serialize( group );

 // now the CapacitatedFacilityLocationBlock data - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim nn = group.addDim( "NNodes" , get_NNodes() );
 netCDF::NcDim na = group.addDim( "NArcs" , get_NArcs() );

 if( get_NNodes() > get_NStaticNodes() )
  group.addDim( "DynNNodes" , get_NNodes() - get_NStaticNodes() );

 if( get_NArcs() > get_NStaticArcs() )
  group.addDim( "DynNArcs" , get_NArcs() - get_NStaticArcs() );

 if( get_MaxNNodes() > get_NStaticNodes() )
  group.addDim( "MaxDynNNodes" , get_MaxNNodes() - get_NStaticNodes() );

 if( get_MaxNArcs() > get_NStaticArcs() )
  group.addDim( "MaxDynNArcs" , get_MaxNArcs() - get_NStaticArcs() );

 ( group.addVar( "SN" , netCDF::NcUint64() , na ) ).putVar( SN.data() );

 ( group.addVar( "EN" , netCDF::NcUint64() , na ) ).putVar( EN.data() );

 if( ! C.empty() )
  ( group.addVar( "C" , netCDF::NcDouble() , na ) ).putVar( C.data() );

 if( ! U.empty() )
  ( group.addVar( "U" , netCDF::NcDouble() , na ) ).putVar( U.data() );

 if( ! B.empty() )
  ( group.addVar( "B" , netCDF::NcDouble() , nn ) ).putVar( B.data() );

 }  // end( CapacitatedFacilityLocationBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::chg_costs( c_Vec_CNumber_it NCost , Range rng ,
			  c_ModParam issueMod , c_ModParam issueAMod )
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
    themselves from Variable that are going to be deleted anyway. Then
    deletes all the "abstract representation", if any. */

 // clear the bound constraints
 for( auto & cnst : UB )
  cnst.clear();
 for( auto & cnst : dUB )
  cnst.clear();

 // clear the flow conservation constraints
 for( auto & cnst : E )
  cnst.clear();
 for( auto & cnst : dE )
  cnst.clear();

 // then delete them all
 dUB.clear();
 UB.clear();
 dE.clear();
 E.clear();

 // clear the objective function
 c.clear();

 // delete all Variable
 dx.clear();
 x.clear();

 // explicitly reset all Constraint and Variable
 // this is done for the case where this method is called prior to re-loading
 // a new instance: if not, the new representation would be added to the
 // (no longer 
 reset_static_constraints();
 reset_static_variables();
 reset_dynamic_constraints();
 reset_dynamic_variables();
 reset_objective();

 AR = 0;

 }  // end( CapacitatedFacilityLocationBlock::guts_of_destructor )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
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

 throw( std::invalid_argument( "unsupported Modification to CapacitatedFacilityLocationBlock" ) );

 }  // end( CapacitatedFacilityLocationBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void CapacitatedFacilityLocationBlock::compute_conditional_bounds( void )
{
 f_cond_lower = f_cond_upper = 0;

 auto tC = C.begin();
 auto tU = U.begin();

 for( ; tC < C.end() ; ++tC , ++tU ) {
  if( *tC == 0 )
   continue;

  if( *tC < 0 ) {
   if( *tU == Inf<FNumber>() ) {
    f_cond_lower = - Inf<double>();
    break;
    }
   else
    f_cond_lower += *tC * (*tU);
   }
  else
   if( *tU == Inf<FNumber>() ) {
    f_cond_upper = Inf<double>();
    break;
    }
   else
    f_cond_upper += *tC * (*tU);
   }

 if( f_cond_lower > - Inf<double>() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC < 0 ) {
    if( *tU == Inf<FNumber>() ) {
     f_cond_lower = - Inf<double>();
     break;
     }
    else
     f_cond_lower += *tC * (*tU);
    }
  }

 if( f_cond_upper < Inf<double>() ) {
  for( ; tC < C.end() ; ++tC , ++tU )
   if( *tC > 0 ) {
    if( *tU == Inf<FNumber>() ) {
     f_cond_upper = Inf<double>();
     break;
     }
    else
     f_cond_upper += *tC * (*tU);
    }
  }
 }  // end( CapacitatedFacilityLocationBlock::compute_conditional_bounds )

/*--------------------------------------------------------------------------*/

ModParam CapacitatedFacilityLocationBlock::make_amod_param( ModParam issueAMod , Index num )
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

void CapacitatedFacilityLocationBlock::unmake_amod_param( ModParam oldiAM , ModParam newiAM ,
				  Index num )
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

 // check variables, these are always there - - - - - - - - - - - - - - - - -
 if( x.size() != get_NStaticArcs() )
  std::cerr << "x.size() != NStaticArcs" << std::endl;

 if( dx.size() < get_NArcs() - get_NStaticArcs() )
  std::cerr << "dx.size() too small" << std::endl;

 // check flow constraints- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasFlw ) { 
  if( E.size() != get_NStaticNodes() )
   std::cerr << "E.size() != NStaticNodes" << std::endl;

  if( dE.size() < get_NNodes() - get_NStaticNodes() )
   std::cerr << "dE.size() too small" << std::endl;
  
  Subset NRIncid( get_NNodes() , 0 );

  Index objbnd = 0;          // number of actives between obj and bound
  if( AR & HasBnd )
   ++objbnd;
  if( AR & HasObj )
   ++objbnd;
  Index expnc = 2 + objbnd;  // total number of active stuff per variable

  // static arcs
  Index a = 0;
  for( auto xi = x.begin() ; a < get_NStaticArcs() ; ++a , ++xi ) {
   if( is_deleted( a ) ) {
    std::cerr << "static arc " << a << " deleted" << std::endl;
    continue;
    }
   
   if( xi->get_num_active() != expnc )
    std::cerr << "arc " << a << " active in " << xi->get_num_active()
	      << " constraints" << std::endl;

   Index sna = SN[ a ];
   if( ! sna )
    std::cerr << "SN[ " << a << " ] == 0" << std::endl;
   else
    --sna;
   if( sna >= get_NArcs() )
    std::cerr << "SN[ " << a << " ] == " << sna << " >= |A!" << std::endl;

   if( ( SN[ a ] > 0 ) && ( sna < get_NArcs() ) ) {
    ++NRIncid[ sna ];
    auto snc = get_lfc( i2p_e( sna ) );
    auto sni = snc->is_active( &(*xi) );
    if( sni >= snc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for SN[ a ] == "
	       << sna << std::endl;
    }
    
   Index ena = EN[ a ];
   if( ! ena )
    std::cerr << "EN[ " << a << " ] == 0" << std::endl;
   else
    --ena;
   if( ena >= get_NArcs() )
    std::cerr << "EN[ " << a << " ] == " << ena << " >= |A!" << std::endl;

   if( ( EN[ a ] > 0 ) && ( ena < get_NArcs() ) ) {
    ++NRIncid[ ena ];
    auto enc = get_lfc( i2p_e( ena ) );
    auto eni = enc->is_active( &(*xi) );
    if( eni >= enc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for EN[ a ] == "
	       << ena << std::endl;
    }
   }

  // dynamic arcs
  for( auto xi = dx.begin() ; a < get_NArcs() ; ++a , ++xi ) {
   if( is_deleted( a ) ) {
    if( xi->get_num_active() != objbnd )
     std::cerr << "deleted arc " << a << " active in "
	       << xi->get_num_active() << " constraints" << std::endl;
    continue;
    }
   
   if( xi->get_num_active() != expnc )
    std::cerr << "arc " << a << " active in " << xi->get_num_active()
	      << " constraints" << std::endl;

   Index sna = SN[ a ];
   if( ! sna )
    std::cerr << "SN[ " << a << " ] == 0 " << std::endl;
   else
    --sna;
   if( sna >= get_NArcs() )
    std::cerr << "SN[ " << a << " ] == " << sna << " >= |A!" << std::endl;

   if( ( SN[ a ] > 0 ) && ( sna < get_NArcs() ) ) {
    ++NRIncid[ sna ];
    auto snc = get_lfc( i2p_e( sna ) );
    auto sni = snc->is_active( &(*xi) );
    if( sni >= snc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for SN[ a ] == "
	       << sna << std::endl;
    }

   Index ena = EN[ a ];
   if( ! ena )
    std::cerr << "EN[ " << a << " ] == 0 " << std::endl;
   else
    --ena;
   if( ena >= get_NArcs() )
    std::cerr << "EN[ " << a << " ] == " << ena << " >= |A!" << std::endl;

   if( ( EN[ a ] > 0 ) && ( ena < get_NArcs() ) ) {
    ++NRIncid[ ena ];
    auto enc = get_lfc( i2p_e( ena ) );
    auto eni = enc->is_active( &(*xi) );
    if( eni >= enc->get_num_active_var() )
     std::cerr << "static arc " << a
	       << " absent in flow constraint for EN[ a ] == "
	       << ena << std::endl;
    }
   }

  if( HasDynamicX() && is_deleted( get_NArcs() - 1 ) )
   std::cerr << "last dynamic arc is deleted" << std::endl;

  // static nodes
  Index n = 0;
  for( auto ni = E.begin() ; n < get_NStaticNodes() ; ++n , ++ni ) {
   auto lni = get_lfc( &(*ni) );
   if( lni->get_num_active_var() != NRIncid[ n ] )
    std::cerr << "active variables in static flow constraint " << n
	      << " == " << lni->get_num_active_var()
	      << " do not match with incident arcs " << NRIncid[ n ]
	      << std::endl;
   }
   
  // dynamic nodes
  for( auto ni = dE.begin() ; n < get_NNodes() ; ++n , ++ni ) {
   auto lni = get_lfc( &(*ni) );
   if( lni->get_num_active_var() != NRIncid[ n ] )
    std::cerr << "active variables in dynamic flow constraint " << n
	      << " == " << lni->get_num_active_var()
	      << " do not match with incident arcs " << NRIncid[ n ]
	      << std::endl;
   }
  }  // end( if( AR & HasFlw ) )

 // check bound constraints - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasBnd ) {
  // static bounds
  Index a = 0;
  auto UBi = UB.begin();
  for( auto xi = x.begin() ; xi != x.end() ; ++a , ++xi , ++UBi ) {
   if( UBi->is_active( &(*xi) ) >= UBi->get_num_active_var() )
    std::cerr << "static arc " << a << " absent in bound constraint"
	      << std::endl;
   }

 // dynamic bounds
  auto dUBi = dUB.begin();
  for( auto xi = dx.begin() ; xi != dx.end() ; ++a , ++xi , ++dUBi ) {
   if( is_deleted( a ) )
    continue;

   if( dUBi->is_active( &(*xi) ) >= dUBi->get_num_active_var() )
    std::cerr << "dynamic arc " << a << " absent in bound constraint"
	      << std::endl;
   }
  }  // end( if( AR & HasBnd ) )

 // check objective - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( AR & HasObj ) {
  auto lfo = get_lfo();
  if( lfo->get_num_active_var() != get_NArcs() )
   std::cerr << "objective has " << lfo->get_num_active_var()
	     << " variables while |A| = " << get_NArcs() << std::endl;

  Index a = 0;
  for( auto & xi : x ) {
   if( lfo->is_active( & xi ) >= get_NArcs() )
    std::cerr << "static arc " << a << " absent from objective " << std::endl;
   ++a;
   }

  for( auto & xi : dx ) {
   if( lfo->is_active( & xi ) >= get_NArcs() )
    std::cerr << "dynamic arc " << a << " absent from objective "
	      << std::endl;
   ++a;
   }
  }  // end( if( AR & HasObj ) )
 }  // end( CapacitatedFacilityLocationBlock::CheckAbsVSPhys )

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF MCFSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void MCFSolution::deserialize( const netCDF::NcGroup & group )
{
 std::vector<size_t> start = { 0 };

 netCDF::NcDim na = group.getDim( "NumArcs" );
 if( na.isNull() )
  v_x.clear();
 else {
  netCDF::NcVar fs = group.getVar( "FlowSolution" );
  if( fs.isNull() )
   v_x.clear();
  else {
   v_x.resize( na.getSize() );
   fs.getVar( v_x.data() );
   }
  }

 netCDF::NcDim nn = group.getDim( "NumNodes" );
 if( nn.isNull() )
  v_pi.clear();
 else {
  netCDF::NcVar ps = group.getVar( "Potentials" );
  if( ps.isNull() )
   v_pi.clear();
  else {
   v_pi.resize( nn.getSize() );
   ps.getVar( v_pi.data() );
   }
  }
 }  // end( MCFSolution::deserialize )

/*--------------------------------------------------------------------------*/

void MCFSolution::read( const Block * const block )
{
 auto MCFB = dynamic_cast< const CapacitatedFacilityLocationBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a CapacitatedFacilityLocationBlock" ) );

 // read flows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( v_x.size() < MCFB->get_NArcs() )
  v_x.resize( MCFB->get_NArcs() );

 auto vxi = v_x.begin();

 // static part
 for( auto & xi : MCFB->x )
  *(vxi++) = xi.get_value();

 // dynamic part
 for( auto & xi : MCFB->dx )
  *(vxi++) = xi.get_value();

 // read potentials - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MCFB->E.empty() && MCFB->dE.empty() )  // no potentials available
  return;

 if( v_pi.size() < MCFB->get_NNodes() )
  v_pi.resize( MCFB->get_NNodes() );

 auto vpii = v_pi.begin();
 
 // static part
 for( auto & ei : MCFB->E )
  *(vpii++) = ei.get_dual();

 // dynamic part
 for( auto & ei : MCFB->dE )
  *(vpii++) = ei.get_dual();

 }  // end( MCFSolution::read )

/*--------------------------------------------------------------------------*/

void MCFSolution::write( Block * const block ) 
{
 auto MCFB = dynamic_cast< CapacitatedFacilityLocationBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument( "block is not a CapacitatedFacilityLocationBlock" ) );

 // write flows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_x.empty() ) {
  if( v_x.size() < MCFB->get_NStaticArcs() )
   throw( std::invalid_argument( "incompatible flow size" ) );

  auto vxi = v_x.begin();

  // static part
  for( auto & xi : MCFB->x )
   xi.set_value( *(vxi++) );

  // dynamic part
  for( auto dxi = MCFB->dx.begin() ;
       ( dxi != MCFB->dx.end() ) && ( vxi != v_x.end() ) ; )
   (*(dxi++)).set_value( *(vxi++) );
  }

 // write potentials- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( v_pi.empty() )  // no potentials to write
  return;

 if( MCFB->E.empty() && MCFB->dE.empty() )  // no Constraint to write to
  return;

 if( v_pi.size() < MCFB->get_NStaticNodes() )
  throw( std::invalid_argument( "incompatible potential size" ) );

 auto vpii = v_pi.begin();

 // static part
 for( auto & ei : MCFB->E )
  ei.set_dual( *(vpii++) );

 // dynamic part
 for( auto dei = MCFB->dE.begin() ;
      ( dei != MCFB->dE.end() ) && ( vpii != v_pi.end() ) ; )
  (*(dei++)).set_dual( *(vpii++) );

 // write reduced costs (if any)- - - - - - - - - - - - - - - - - - - - - - -
  
 if( MCFB->UB.empty() && MCFB->dUB.empty() )  // no bounds to write to
  return;

 CapacitatedFacilityLocationBlock::Index i = 0;

 // static part
 for( auto ubi = MCFB->UB.begin() ; ubi != MCFB->UB.end() ; ++i )
  (ubi++)->set_dual( MCFB->get_C( i ) + v_pi[ MCFB->SN[ i ] - 1 ]
		                      - v_pi[ MCFB->EN[ i ] - 1 ] );
 // dynamic part
 for( auto dubi = MCFB->dUB.begin() ;
      ( dubi != MCFB->dUB.end() ) && ( i < MCFB->get_NArcs() ) ; ++i )
  (dubi++)->set_dual( MCFB->get_C( i ) + v_pi[ MCFB->SN[ i ] - 1 ]
		                       - v_pi[ MCFB->EN[ i ] - 1 ] );
 }  // end( MCFSolution::write )

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
