/*--------------------------------------------------------------------------*/
/*---------------------------- File txt2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing CapacitatedFacilityLocationBlock instance
 * files, be them netCDF ones of text ones in the "standard" ORLib one which
 * is somehow "customer-oriented), out of text files in three different
 * formats, one of which is the the ORLib one implemented in
 * CapacitatedFacilityLocationBlock::load( istream ) and the other two
 * slightly different versions of "facility-oriented" formats found around
 * the web.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include <fstream>

#include <CapacitatedFacilityLocationBlock.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

using Index = Block::Index;

using Demand = CapacitatedFacilityLocationBlock::Demand;

using DVector = CapacitatedFacilityLocationBlock::DVector;

using Cost = CapacitatedFacilityLocationBlock::Cost;

using CVector = CapacitatedFacilityLocationBlock::CVector;

using CMatrix = CapacitatedFacilityLocationBlock::CMatrix;


/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 if( argc < 4 ) {
  cerr << "Usage: " << argv[ 0 ] << " frmt txt_in file_out" << endl
       << "        frmt: 0 = ORLib, 1 = demands-first, 2 demand-last" << endl
       << "        if file_out ends in .dmx a netCDF file is created" << endl
       << "        otherwise a text file is created" << endl;
  return( 1 );
  }

 Index frmt;
 Str2Sthg( argv[ 1 ] , frmt );
 if( frmt > 2 ) {
  cerr << "Error: unknown file format " << argv[ 1 ] << endl;
  return( 1 );  
  }

 // open input file in text format
 ifstream ProbFile( argv[ 2 ] );
 if( ! ProbFile.is_open() ) {
  cerr << "Error: cannot open file " << argv[ 2 ] << endl;
  return( 1 );
  }

 // create the CapacitatedFacilityLocationBlock
 CapacitatedFacilityLocationBlock CFLB;

 if( ! frmt )  // standard ORLib format- - - - - - - - - - - - - - - - - - - -
  ProbFile >> CFLB;  // just use the laod() method
 else {  // any of the two facility-oriented formats - - - - - - - - - - - - -
  // the two facility oriented formats are similar, starting with
  //
  // m: number of potential locations for facilities
  // n: number of customers.
  //
  // then, for demands-first
  //
  // D_1 D_2 D_3 ... D_n  (customers demand)
  // Q_1 Q_2 Q_3 ... Q_m  (facilities capacity)
  // F_1 F_2 F_3 ... f_m  (facilities cost)
  //
  // while for demands-last
  //
  // for i = 1 ... m
  //     Q_i  F_i
  // D_1 D_2 D_3 ... D_n
  //
  // then, for both formats we have
  //
  // C_{11} C_{12} C_{13} ... C_{1n} 
  // C_{21} C_{22} C_{23} ... C_{2n} 
  //   :     :      :     ...  :
  // C_{m1} C_{m2} c_{m3} ... C_{mn} 
  //
  // i.e., the transportation costs arranged facility-wise, but:
  //
  // - for demands-first this is the unitary transportation cost
  // - for demands-last this is the total cost of serving customer j out of i
  //
  // note that in both cases NO COMMENTS ARE ALLOWED

  Index m;
  ProbFile >> m;
  if( ProbFile.fail() ) {
   cerr << "Error reading number of facilities" << endl;
   return( 1 );
   }

  Index n;
  ProbFile >> n;
  if( ProbFile.fail() ) {
   cerr << "Error reading number of customers" << endl;
   return( 1 );
   }

  DVector D( n );
  DVector Q( m );
  CVector F( m );

  if( frmt == 1 ) {  // demands-first
   // first read demands (doh!)
   for( Index j = 0 ; j < n ; ++j ) {
    ProbFile >> D[ j ];
    if( ProbFile.fail() ) {
     cerr << "Error reading demand " << j << endl;
     return( 1 );
     }
    }

   // then read capacities
   for( Index i = 0 ; i < m ; ++i ) {
    ProbFile >> Q[ i ];
    if( ProbFile.fail() ) {
     cerr << "Error reading capacity " << i << endl;
     return( 1 );
     }
    }

   // then read fixed costs
   for( Index i = 0 ; i < m ; ++i ) {
    ProbFile >> F[ i ];
    if( ProbFile.fail() ) {
     cerr << "Error reading fixed cost " << i << endl;
     return( 1 );
     }
    }
   }
  else {             // demands-last
   // first read pairs ( capacity , fixed cost )
   for( Index i = 0 ; i < m ; ++i ) {
    ProbFile >> Q[ i ];
    if( ProbFile.fail() ) {
     cerr << "Error reading capacity " << i << endl;
     return( 1 );
     }
    ProbFile >> F[ i ];
    if( ProbFile.fail() ) {
     cerr << "Error reading fixed cost " << i << endl;
     return( 1 );
     }
    }

   // last read demands (doh!)
   for( Index j = 0 ; j < n ; ++j ) {
    ProbFile >> D[ j ];
    if( ProbFile.fail() ) {
     cerr << "Error reading demand " << j << endl;
     return( 1 );
     }
    }
   }

  // now read transportation costs
  CMatrix C;
  C.resize( boost::extents[ m ][ n ] );

  for( Index i = 0 ; i < m ; ++i )
   for( Index j = 0 ; j < n ; ++j ) {
    ProbFile >> C[ i ][ j ];
    if( ProbFile.fail() ) {
     cerr << "Error reading transportation cost ( " << i << ", " << j
	  << " )" << endl;
     return( 1 );
     }

    if( frmt == 1 )  // in the demands-first case these are unitary costs
     C[ i ][ j ] *= D[ j ]; // make them total costs
    }

  // finally, load everything into the CapacitatedFacilityLocationBlock
  CFLB.load( m , n , move( Q ) , move( F ) , move( D ) , std::move( C ) );
  }

 ProbFile.close();  // input done, close the file

 // now either serialize or << the CapacitatedFacilityLocationBlock
 string name( argv[ 3 ] );
 if( name.size() <= 4 ) {
  cerr << "Error: suffix missing oun output filename" << name << endl;
  return( 1 );
  } 
 
 if( name.compare( name.size() - 4 , 4 , ".txt" ) ) {
  ofstream OutFile( name , ofstream::trunc );
  OutFile << CFLB;
  }
 else
  CFLB.Block::serialize( name , int( eBlockFile ) );
  // why the Block:: is needed completely evades me, but clang++ seems to
  // think it is

 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File txt2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

