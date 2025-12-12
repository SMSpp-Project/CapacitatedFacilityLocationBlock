/*--------------------------------------------------------------------------*/
/*---------------------------- File txt2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing instance files, be them netCDF ones of text
 * ones, of the Capacitated Facility Location problem by reading them into
 * the CapacitatedFacilityLocationBlock (in whatever format it supports) and
 * either print()-ing to a text file (in whatever format this is supported)
 * or deserialize()-int into a netCDF file.
 *
 * While this main() is written for CapacitatedFacilityLocationBlock, in
 * fact it does not even include CapacitatedFacilityLocationBlock.h: by
 * having the name used in the factory call (new_Block), say, be provided
 * by in the command line it would work with any other kind of :Block.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include "Block.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 if( argc < 4 ) {
  std::cerr << "Usage: " << argv[ 0 ]
	    << " file_in frmt_in file_out [frmt_out]" << std::endl
	    << "        if file_in ends in .nc4 a netCDF file is assumed"
	    << std::endl
	    << "        else frmt is the parameter of load()"
	    << std::endl
	    << "        if file_out ends in .nc4 a netCDF file is created"
	    << std::endl
	    << "        else frmt_out (if any) is the parameter of print()"
	    << std::endl;
  return( 1 );
  }

 char frmt_in = argv[ 2 ][ 0 ];
 char frmt_out = argc > 4 ? argv[ 4 ][ 0 ] : 'C';

 // either deserialize() or load() the [CapacitatedFacilityLocation]Block
 Block * CFLB;

 std::string iname( argv[ 1 ] );
 if( ( iname.size() > 4 ) &&
     ( iname.substr( iname.size() - 4 , 4 ) == ".nc4" ) )  // a netCDF file
  CFLB = Block::deserialize( iname );
 else {                                                    // a text file
  // construct the [CapacitatedFacilityLocation]Block via the factory
  CFLB = Block::new_Block( "CapacitatedFacilityLocationBlock" );
  CFLB->load( iname , frmt_in );  // now load() it
  }

 // either serialize() or print() the [CapacitatedFacilityLocation]Block
 std::string oname( argv[ 3 ] );
 if( ( oname.size() > 4 ) &&
     ( oname.substr( oname.size() - 4 , 4 ) == ".nc4" ) ) {
  CFLB->Block::serialize( oname , eBlockFile );
  // why the Block:: is needed completely evades me, but clang++ seems to
  // think it is
  }
 else
  CFLB->print( oname , frmt_out );

 // cleanup
 delete CFLB;

 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File txt2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
