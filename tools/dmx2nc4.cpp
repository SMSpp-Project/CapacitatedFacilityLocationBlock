/*--------------------------------------------------------------------------*/
/*---------------------------- File dmx2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing MCFBlock netCDF files out of DIMACS ones.
 * It just creates one MCFBlock and loads it from a DIMACS file, then
 * de-serializes it to a netCDF file.
 *
 * Optionally it hacks into the netCDF file to change the number of static
 * and dynamic nodes and arcs, as well as the maximum number of nodes and
 * arcs (remember that all nodes and arcs loaded out of a DIMACS file are
 * treated as static ones).
 *
 * \version 0.10
 *
 * \date 11 - 04 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include <MCFBlock.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void Str2Sthg( const char* const str , T &sthg )
{
 std::istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 double dyn_node = 0;
 double dyn_arc = 0;
 double max_node = 0;
 double max_arc = 0;
 switch( argc ) {
  case( 6 ): Str2Sthg( argv[ 5 ] , max_arc );
  case( 5 ): Str2Sthg( argv[ 4 ] , max_node );
  case( 4 ): Str2Sthg( argv[ 3 ] , dyn_arc );
  case( 3 ): Str2Sthg( argv[ 2 ] , dyn_node );
  case( 2 ): break;
  default:
   std::cerr << "Usage: " << argv[ 0 ]
	     << " DIMACS_in [ dyn_node dyn_arc max_node max_arc ]"
	     << std::endl
	     << "        if DIMACS_in ends in .dmx the suffix is removed"
	     << std::endl
	     << "        and DIMACS_in[-suffix].nc4 is created"
	     << std::endl
	     << "        any number -1 <= . < 0 is treated as a fraction"
	     << std::endl
	     << "        any positive number as the actual number"
	     << std::endl;
  return( 1 );
  }

 // open input file in DIMACS format
 std::ifstream ProbFile( argv[ 1 ] );
 if( ! ProbFile.is_open() ) {
  std::cerr << "Error: cannot open file " << argv[ 1 ] << std::endl;
  return( 1 );
  }

 // load the MCFBlock from the file
 MCFBlock *MCFB = dynamic_cast< MCFBlock * >( Block::new_Block( "MCFBlock" ) );
 if( ! MCFB ) {
  std::cerr << "Failed to initialize MCFBlock" << std::endl;
  return( 1 );
  }

 ProbFile >> *MCFB;
 ProbFile.close();

 // serialize the MCFBlock (in a netCDF BlockFile)
 std::string name( argv[ 1 ] );
 if( name.size() > 4 ) {
  std::string sffx = name.substr( name.size() - 4 , 4 );
  std::string dmx( ".dmx" );

  if( std::equal( sffx.begin() , sffx.end() , dmx.begin() ,
		  []( auto a , auto b ){
		   return( std::tolower( a ) == std::tolower( b ) ); } ) )
   name.erase( name.size() - 4 , 4 );
  }

 name.append( ".nc4" );
 
 // first open the file
 netCDF::NcFile f( name , netCDF::NcFile::replace );

 // put in the type attribute
 f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );

 // serialize the MCFBlock into Block_0
 MCFB->Block::serialize( f , eBlockFile );

 // if it has to be fully static, all is done
 if( ( dyn_node == 0 ) && ( dyn_arc == 0 ) &&
     ( max_node == 0 ) && ( max_arc == 0 ) )
  return( 0 );

 netCDF::NcGroup bg = f.getGroup( "Block_0" );

 if( bg.isNull() ) {
  std::cout << "Error: can't find Block_0" << std::endl;
  return( 1 );
  }

 if( dyn_node != 0 ) {
  MCFBlock::Index dn = MCFB->get_NNodes();
  dn = dyn_node > 0 ? std::min( dn , MCFBlock::Index( dyn_node ) )
                    : MCFBlock::Index( - dyn_node * dn );
  bg.addDim( "DynNNodes" , dn );
  }

 if( dyn_arc != 0 ) {
  MCFBlock::Index da = MCFB->get_NArcs();
  da = dyn_arc > 0 ? std::min( da , MCFBlock::Index( dyn_arc ) )
                   : MCFBlock::Index( - dyn_arc * da );
  bg.addDim( "DynNArcs" , da );
  }

 if( max_node != 0 ) {
  MCFBlock::Index mn = MCFB->get_NNodes();
  mn = max_node > 0 ? std::min( mn , MCFBlock::Index( max_node ) )
                    : MCFBlock::Index( - max_node * mn );
  bg.addDim( "MaxDynNNodes" , mn );
  }

 if( max_arc != 0 ) {
  MCFBlock::Index ma = MCFB->get_NArcs();
  ma = max_arc > 0 ? std::min( ma , MCFBlock::Index( max_arc ) )
                   : MCFBlock::Index( - max_arc * ma );
  bg.addDim( "MaxDynNArcs" , ma );
  }
 
 // delete the MCFBlock
 delete MCFB;

 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File dmx2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

