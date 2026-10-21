/*--------------------------------------------------------------------------*/
/*--------------- File CapacitatedFacilityLocationBlock.cppm ---------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.capacitatedfacilitylocationblock, which
 * exports every name that the headers of CapacitatedFacilityLocationBlock
 * declare at namespace scope, together with the modules it depends on; a
 * translation unit can then write
 *
 *     import smspp.capacitatedfacilitylocationblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoit Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "CapacitatedFacilityLocationBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.capacitatedfacilitylocationblock;

export import smspp.binaryknapsackblock;
export import smspp.mcfblock;

export namespace SMSpp_di_unipi_it {

 // CapacitatedFacilityLocationBlock.h
 using SMSpp_di_unipi_it::CapacitatedFacilityLocationBlock;
 using SMSpp_di_unipi_it::CapacitatedFacilityLocationBlockMod;
 using SMSpp_di_unipi_it::CapacitatedFacilityLocationBlockRngdMod;
 using SMSpp_di_unipi_it::CapacitatedFacilityLocationBlockSbstMod;
 using SMSpp_di_unipi_it::CapacitatedFacilityLocationSolution;
}

/*--------------------------------------------------------------------------*/
/*------------- End File CapacitatedFacilityLocationBlock.cppm -------------*/
/*--------------------------------------------------------------------------*/
