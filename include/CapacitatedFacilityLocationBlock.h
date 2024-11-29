/*--------------------------------------------------------------------------*/
/*-------------- File CapacitatedFacilityLocationBlock.h -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the concrete class CapacitatedFacilityLocationBlock, which
 * implements the Block concept [see Block.h] for the "basic version" of the
 * Capacitated Facility Location (CFL) problem, a.k.a. the Capacitated
 * Warehouse Location (CWL) problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CapacitatedFacilityLocationBlock
 #define __CapacitatedFacilityLocationBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BinaryKnapsackBlock.h"

#include "LinearFunction.h"

#include "MCFBlock.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class CapacitatedFacilityLocationBlock;
 // forward declaration of CapacitatedFacilityLocationBlock

 class CapacitatedFacilityLocationSolution;
 // forward declaration of CapacitatedFacilityLocationSolution

 class BinaryKnapsackBlockMod;
 // forward declaration of BinaryKnapsackBlockMod

 class MCFBlock;     // forward declaration of MCFBlock

 class MCFBlockMod;  // forward declaration of MCFBlockMod

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CapacitatedFacilityLocationBlock_CLASSES
 *            Classes in CapacitatedFacilityLocationBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------- CLASS CapacitatedFacilityLocationBlock --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the CFL / CWL problem
/** The CapacitatedFacilityLocationBlock class implements the Block concept
 * [see Block.h] for a "pretty basic version" of the Capacitated Facility
 * Location (CFL) problem, a.k.a. the Capacitated Warehouse Location (CWL)
 * problem.
 *
 * The data of the problem consist of a set of m Locations where to install
 * Facilities / Warehouses whose role is to serve a set of n customers,
 * usually (but not necessarily) with n >> m. For ease of notation we define
 * the set I = { 1 , ... m } of Facilities / Warehouses and the set J =
 * { 1 , ... , n } of customers. Each Facility [ / Warehouse ] i \in I has a
 * capacity Q[ i ] corresponding to the maximum amount that it can deliver of
 * the unique commodity that the customers require; however, this can be done
 * only if the Facility [ ... ] is open, which has a fixed cost F[ i ].
 * Underlying the problem there is a complete bipartite graph G = ( I x J ,
 * A ) with m + n nodes and m * n = (directed) arcs corresponding to
 * (directed) routes from each facility i \in I to each customer j \in J.
 * Each arc ( i , j ) has a linear cost coefficient C[ i , j ] corresponding
 * to the cost of serving all the demand D[ j ] of the user j for the unique
 * commodity from facility i.
 *
 * Introducing flow variables X[ i , j ] that represent the percentage of the
 * total demand of customer j served by facility i and binary variables
 * Y[ i ] that represent (in the obvious way) if the facility i is opened, a
 * "natural formulation" of the problem is:
 * \f[
 *  \min \sum_{ i \in I } \sum_{ j \in J } C[ i , j ] X[ i , j ] +
 *       \sum_{ i \in I } F[ i ] Y[ i ]
 * \f]
 * \f[
 *  \sum_{ i \in I } X[ i , j ] = 1                        \quad j \in J  (1)
 * \f]
 * \f[
 *  \sum_{ j \in J } D[ j ] X[ i , j ] \leq Q[ i ] Y[ i ]  \quad i \in I  (2)
 * \f]
 * \f[
 *  Y[ i ] \in \{ 0 , 1 \}                                 \quad i \in I  (3)
 * \f]
 * \f[
 *  0 \leq X[ i , j ] \leq 1                     \quad j \in J , i \in I  (4)
 * \f]
 * The n equations (1) impose that all the demand of each customer is
 * satisfied. The m inequalities (2) impose that a facility is only used to
 * serve the customers if it has been opened, and in this case only up to its
 * capacity. (3) and (4) are the bounds and integrality restrictions on the
 * variables. Note that the above is the *splittable* version of the problem
 * whereby a customer can be served by multiple facilities; with the "minor"
 * change in the model where (4) is replaced with
 * \f[
 *  X[ i , j ] \in \{ 0 , 1 \}                  \quad j \in J , i \in I  (4')
 * \f]
 * one can represent the *unsplittable* version where each customer need be
 * served by exactly one facility.
 *
 * This "basic" version of the problem is typically complicated in reality
 * by many further considerations: multiple commodities, service to customers
 * during multiple time instants along a given time horizon, with varying
 * demand possibly subject to uncertainties, more complex installation and/or
 * handling costs for facilities (if a facility incurs a unitary cost H[ j ]
 * for handling one unit of commodity this is considered included in the
 * transportation cost C[ i , j ] here, but this is possible only if the
 * handling cost is linear, which it may not be).
 *
 * This class only represent the "basic version" and it is primarily intended
 * as a "didactic" implementation for showing some of the features of SMS++,
 * among which:
 *
 * - CapacitatedFacilityLocationBlock supports a number of different
 *   formulations of the problem, among which ones suitable for the use of
 *   decomposition approaches; see the comments to
 *   generate_abstract_variables();
 *
 * - CapacitatedFacilityLocationBlock supports reformulations/relaxations of
 *   the problem via the "R3Block" mechanism; see the comments to
 *   get_R3_Block().
 *
 * The current implementation of the class is only an initial version, and
 * features are expected to be added along time.
 *
 * As such, the class currently1 only supports only a subset of the possible
 * operations on the data of the problem. This starts with the fact that
 *
 *     THE SIZE OF THE PROBLEM IS STATIC AND CAN NEVER BE CHANGED
 *     (save if the problem is completely re-loaded)
 *
 * which implies that there are no dynamic Variable and Constraint, which in
 * turn dramatically simplifies a lot of the underlying logic. Furthermore:
 *
 * - Changing customers' demands via the abstract representation is not
 *   allowed in the Standard Formulation and the Knapsack Formulation, since
 *   the same demand is replicated in multiple constraints; one could ask
 *   that all the changes happen at the same time and that the corresponding
 *   Modification are bunched together in a GroupModification, which is
 *   possible but complex and not implemented yet. The change is instead
 *   possible in the Flow Formulation where demands are node deficits.
 *
 * - Changing the splittable/unsplittable form of the problem, i.e., the
 *   integrality of all variables x[ i ][ j ], via the abstract
 *   representation is never allowed.
 *
 * - map_[forward/back]_[Modification/Solution]() are fully implemented for   
 *   both types of R3Block, except "back Modification" that is not
 *   implemented for the MCF R3Block.
 *
 * - The unsplittable version of the problem cannot be represented when
 *   using the Flow Formulation (this is inherent and unlikely to ever
 *   change).
 *
 * - The Modification corresponding to changing the type of the problem are
 *   ignored in map_forward_Modification() to a MCF R3Block, since this
 *   represents a continuous relaxation of the original CFL and therefore it
 *   is identical in the splittable and unsplittable case. */

class CapacitatedFacilityLocationBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 using Demand = double;  ///< demands of customers = capacity of facilities

 using DVector = std::vector< Demand >;  ///< a vector of demands

 using c_DVector = const DVector;       ///< a const vector of demands

 using DV_it = DVector::iterator;       ///< iterator into a DVector

 using c_DV_it = DVector::const_iterator;
 ///< const_iterator into a DVector

 using Cost = double;                  ///< transportation / fixed costs

 using CVector = std::vector< Cost >;  ///< a vector of costs

 using c_CVector = const CVector;      ///< a const vector of costs

 using CV_it = CVector::iterator;      ///< iterator into a CVector

 using c_CV_it = CVector::const_iterator;
 ///< const_iterator into a CVector

 using CMatrix = boost::multi_array< Cost , 2 >;
 ///< a 2-dimensional matrix of costs

 using c_CMatrix = const CMatrix;
 ///< a const 2-dimensional matrix of costs

 using IntSolution = std::vector< bool >;   ///< an integer (binary) solution

 using IS_it = IntSolution::iterator;       ///< iterator into a IntSolution

 using c_IS_it = IntSolution::const_iterator;
 ///< const_iterator into a IntSolution

 using CntSolution = std::vector< double >;  ///< a continuous solution

 using CS_it = CntSolution::iterator;        ///< iterator into a CntSolution

 using c_CS_it = CntSolution::const_iterator;
 ///< const_iterator into a CntSolution

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend CapacitatedFacilityLocationSolution;
 ///< make CapacitatedFacilityLocationSolution friend

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTANTS ---------------------------------*/
/*--------------------------------------------------------------------------*/

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor, taking a pointer to the father (generic) Block
 /** Constructor. It accepts a pointer to the father Block, which can be of
  * any type, defaulting to nullptr so that this can also be used as the void
  * constructor. */

 explicit CapacitatedFacilityLocationBlock( Block *father = nullptr )
  : Block( father ) , f_n_facilities( 0 ) , f_n_customers( 0 ) ,
    f_unsplittable( false ) , AR( 0 ) , f_cond_lower( dNaN ) ,
    f_cond_upper( dNaN ) , f_mod_skip( false ) {}

/*--------------------------------------------------------------------------*/
 /// destructor; deletes the abstract representation, if any

 virtual ~CapacitatedFacilityLocationBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the CFL instance from memory, moving the input data
 /** Loads the CFL instance from memory. The parameters are what you expect:
  *
  * - m is the number of facilities
  *
  * - n is the number of customers
  *
  * - Q is the std::vector< Demand > of the facility capacities, that must
  *   have exactly size() == m
  *
  * - F is the std::vector< FCost > of the facility opening costs, that must
  *   have exactly size() == m
  *
  * - D is the std::vector< Demand > of the customers demand, that must have
  *   exactly size() == n
  *
  * - C is the boost::multi_array< TCost , 2 > matrix of transportation costs,
  *   arranged facility-wise; this means that C[ i ] for i \in I is a n-vector
  *   so that C[ i ][ j ] is the total transportation cost between facility i
  *   and customer j, i.e., the cost of serving all the demand D[ j ] of
  *   customer j out of facility i.
  *
  * - unsplt is the bool that, if true [default], denotes that the
  *   unsplittable version of of the problem needs be solved.
  *
  * As the && tells, all the data becomes property of the
  * CapacitatedFacilityLocationBlock.
  *
  * Like load( std::istream & ), if there is any Solver attached to this
  * CapacitatedFacilityLocationBlock then a NBModification (the "nuclear
  * option") is issued. */

 void load( Index m , Index n , DVector && Q , CVector && F ,
	    DVector && D , CMatrix && C , bool unsplt = false );

/*--------------------------------------------------------------------------*/
 /// loads the CFL instance from memory, copying the input data
 /** Loads the CFL instance from memory. The parameters are the same as
  * the other form of load() except that (by being const & rather than &&)
  * the vectors/matrices are copied rather than moved.  */

 void load( Index m , Index n , c_DVector & Q , c_CVector & F ,
	    c_DVector & D , c_CMatrix & C , bool unsplt = false ) {
  load( m , n , DVector( Q ) , CVector( F ) , DVector( D ) , CMatrix( C ),
	unsplt );
  }

/*--------------------------------------------------------------------------*/
 /// loads the CFL instance from file in standard format
 /** Loads a CapacitatedFacilityLocationBlock out of a std::istream (which is
  * what operator>> is dispatched to). The std::istream is assumed to contain
  * the description of a CFL instance, which according to \p frmt can be in
  * three different formats:
  *
  * The default (frmt == 0 or frmt == 'C') is the ORLib format", which is
  * somehow customer-oriented and is the following:
  *
  * number of potential facility locations (m)
  * number of customers (n)
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     capacity of facility i, if opened
  *     fixed cost to open facility i
  *
  * for each customer j (j = 1, ..., n):
  *     demand of customer j
  *     for each potential facility location i (i = 1, ..., m): 
  *         cost of allocating all of the demand of j to facility i
  *
  * frmt == 'F' is the facility oriented, demands-first format:
  *
  * number of potential facility locations (m)
  * number of customers (n)
  *
  * for each customer j (j = 1, ..., n):
  *     demand of customer j
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     capacity of facility i, if opened
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     capacity of facility i, if opened
  *     fixed cost to open facility i
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     for each customer j (j = 1, ..., n):
  *         unitary transportation cost from facility i to customer j
  *
  * frmt == 'L' is the facility oriented, demands-last format:
  *
  * number of potential facility locations (m)
  * number of customers (n)
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     capacity of facility i, if opened
  *     fixed cost to open facility i
  *
  * for each customer j (j = 1, ..., n):
  *     demand of customer j
  *
  * for each potential facility location i (i = 1, ..., m): 
  *     for each customer j (j = 1, ..., n):
  *         cost of allocating all of the demand of j to facility i
  *
  * In all the cases (possibly unlikely what the original formats assumed),
  * comments (starting with '#' and taking up to the following newline) can
  * be placed anywhere in the file and are skipped.
  *
  * Like load( memory ), if there is any Solver attached to this
  * CapacitatedFacilityLocationBlock then a NBModification (the "nuclear
  * option") is issued.
  *
  * Note that the text input format does not allow to specify whether the
  * splittable or unsplittable version of the problem is solved; by default
  * the splittable one is assumed. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a CapacitatedFacilityLocationBlock. Besides what is managed by the
  * serialize() method of the base Block class, the group should contain the
  * following:
  *
  * - the dimension "NFacilities" containing the number of facilities
  *
  * - the dimension "NCustomers" containing the number of customers
  *
  * - the variable "FacilityCapacity", of type double and indexed over the
  *   dimension "NFacilities"; the j-th entry of the variable is assumed to
  *   contain the capacity of the j-th facility
  *
  * - the variable "FacilityCost", of type double and indexed over the
  *   dimension "NFacilities"; the j-th entry of the variable is assumed to
  *   contain the opening cost of the j-th facility
  *
  * - the variable "CustomerDemand", of type double and indexed over the
  *   dimension "NCustomers"; the i-th entry of the variable is assumed to
  *   contain the demand of the i-th customer
  *
  * - the variable "TransportationCost", of type double and indexed over
  *   both the dimensions "NFacilities" and "NCustomers"; the entry ( i , j )
  *   is assumed to contain the *total* cost of serving customer i from
  *   facility j
  *
  * All the dimensions and variables are mandatory. However, the optional
  *
  * - variable "FacilityFix", of type char and indexed over the dimension
  *   "NFacilities"; the j-th entry of the variable is assumed to contain
  *   0 (yFree) if the facility is not fixed, 1 (yFix0) if the facility is
  *   closed, and 2 (yFix1) if the facility is fixed open
  *
  * may also be present to represent if any of the facilities is either
  * closed or fixed open (if not, all facilities are considered free). Also,
  * the optional
  *
  * - dimension "UnSplittable" can be present; if so, and it contains a
  *   nonzero value, then the problem is considered an unsplittable one
  *   (each customer must be served by one and only one facility).
  *
  * Otherwise (the dimension is not there or it contains zero) then the
  * problem is considered an splittable one (customer can be served by any
  * number of facilities). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the Capacitated Facility Location
 /** Method that generates the abstract Variable of the Capacitated Facility
  * Location problem, meanwhile deciding which of the different formulations
  * of the problem is produced as the "abstract representation" of the
  * CapacitatedFacilityLocationBlock. The different possible formulations are
  * represented by a single int value "wf" that is obtained as follows:
  *
  * - if either stvv is not nullptr and it is a SimpleConfiguration< int >,
  *   or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_variables_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< int >, then wf is the f_value of the
  *   SimpleConfiguration< int >
  *
  * - otherwise, wf is 0
  *
  * The list of supported formulations is:
  *
  * - wf & 3 == 0 is the "natural formulation" (NF) always comprising
  *   \f[
  *    \min \sum_{ i \in I } \sum_{ j \in J } C[ i , j ] X[ i , j ] +
  *         \sum_{ i \in I } F[ i ] Y[ i ]
  *   \f]
  *   \f[
  *    \sum_{ i \in I } X[ i , j ] = 1                           j \in J  (1)
  *   \f]
  *   \f[
  *    \sum_{ j \in J } D[ j ] X[ i , j ] \leq Q[ i ] Y[ i ]     i \in I  (2)
  *   \f]
  *   \f[
  *     Y[ i ] \in \{ 0 , 1 \}                                   i \in I  (3)
  *   \f]
  *   \f[
  *    0 \leq X[ i , j ] \leq 1                        j \in J , i \in I  (4)
  *   \f]
  *   If wf & 4 (wf == 4) the "splittable constraints" (4) are replaced with
  *   the "unsplittable constraints"
  *   \f[
  *    X[ i , j ] \in \{ 0 , 1 \}                      j \in J , i \in I  (4')
  *   \f]
  *   This means that the CapacitatedFacilityLocationBlock has no sub-Block,
  *   and the two groups of static variables X[] and Y[]:
  *
  *   = "x", a boost::multi_array< ColVariable , 2 > with sizes
  *     f_n_facilities and f_n_customers, which are of type kPosUnitary
  *     (is_positive() == is_unitary() == true) if wf == 0, and of type
  *     kBinary (in addition, is_unitary() == true) if wf == 4
  *
  *   = "y", a std::vector< ColVariable > of size f_n_facilities, which are
  *     all of type kBinary (is_integer() == is_positive() == is_unitary()
  *     == true)
  *
  *   and no dynamic variables
  *
  * - wf & 3 == 1: the Lagrange-friendly "knapsack formulation" (KF). In
  *   this case, CapacitatedFacilityLocationBlock "grows" f_n_facilities
  *   sub-Block, each of type BinaryKnapsackBlock and with f_n_customers + 1
  *   variables. sub-Block i corresponds to facility i: the first
  *   f_n_customers variables correspond to the transportation variables
  *   X[ j , i ] between i and all the customers (in the natural order),
  *   while the last variable correspond to the design variable Y[ i ].
  *   That is, the i-th knapsack problem is
  *   \f[
  *    \min \sum_{ j \in J } C[ i , j ] X[ i , j ] + F[ i ] Y[ i ]
  *   \f]
   *  \f[
  *    \sum_{ j \in J } D[ j ] X[ i , j ] - Q[ i ] Y[ i ] \leq 0
  *   \f]
  *   \f[
  *     Y[ i ] \in \{ 0 , 1 \}
  *   \f]
  *   \f[
  *    0 \leq X[ i , j ] \leq 1                        j \in J
  *   \f]
  *   Then, wf & 4 (wf == 5) the previous "splittable constraints" are
  *   replaced by the "unsplittable constraints"
  *   \f[
  *    X[ i , j ] \in \{ 0 , 1 \}                      j \in J
  *   \f]
  *   That is, Y[ i ] is always kBinary (is_integer() == is_positive() ==
  *   is_unitary() == true), whereas X[ j ] are kBinary for wf == 5 and
  *   kPosUnitary (is_positive() == is_unitary() == true) if wf == 1.
  *   The linking constraints (1) are the only static group of Constraint
  *   in the CapacitatedFacilityLocationBlock.
  *
  * - wf & 3 >= 2: the Benders-friendly "flow formulation" (FF). In this case,
  *   CapacitatedFacilityLocationBlock "grows" two sub-Block. The first one
  *   only has f_n_facilities kBinary variables corresponding with the
  *   design ones Y[ i ]. The second is instead a MCFBlock representing the
  *   continuous relaxation of the problem as produced by get_R3_Block()
  *   with wr3b == ( wf & 3 ) - 1, except the costs of the "facility arcs"
  *   are set to 0. Then, the CapacitatedFacilityLocationBlock contains the
  *   linking constraints
  *   \f[
  *     arc_flow[ i ] \leq Q[ i ] Y[ i ]                       i \in I
  *   \f]
  *   where arc_flow[ i ] is the flow on the "facility arc" corresponding to
  *   facility i in the MCFBlock. The difference between wf & 3 == 2 (i.e.,
  *   wr3b == 1) and wf & 3 == 3 (i.e., wr3b == 2) is that in the former
  *   case the reformulation is "exact" (the problem is completely equivalent
  *   to that of all the other formulations), while in the second it is
  *   "approximate" in that wr3b == 2 causes the addition of extra high-cost
  *   "slack arcs" that ensure that the instance is always feasible even if
  *   the aggregate demand is larger than the aggregate facility capacity
  *   (see get_R3_Block() for details). In this case, setting wf & 4 true
  *   is not supported in that the flows in the MCFBlock are scaled and there
  *   is no (simple) way to include the required integrality constraints. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the CFL
 /** Method that generates the abstract constraint of the CFL. The actual
  * form of these constraints depend on the formulation, which is decided
  * during the call to generate_abstract_variables(); see the comments to
  * that method for details. Yet, generate_abstract_constraints() allows to
  * only partly generate the abstract constraints according to the single
  * int value "wc" that is obtained as follows:
  *
  * - if either stcc is not nullptr and it is a SimpleConfiguration< int >,
  *   or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_constraints_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< int >, then wc is the f_value of the
  *   SimpleConfiguration< int >
  *
  * - otherwise, wc is 3
  *
  * The meaning of wc is bit-wise: the first bit being 1 means that the
  * customers satisfaction constraints are constructed, while the second bit
  * being 1 means that the capacity constraints are constructed.
  * Note that "constraints X constructed" has different meanings according
  * to which formulation is used, as decided by generate_abstract_variables():
  *
  * - If the "natural formulation" (SF) is used, there are the two explicit
  *   groups of "linear constraints" (FRowConstraint with a LinearFunction)
  *   of the "natural formulation", i.e.,
  *
  *   = "sat", a std::vector< FRowConstraint > of size f_n_customers
  *     imposing the satisfaction of customers' demands (1)
  *
  *   = "cap", a std::vector< FRowConstraint > of size f_n_facilities
  *     imposing the maximum capacity of facilities as well as the logical
  *     constraints that a facility can only be used to serve any customer
  *     if it is open (2)
  *
  *   and no dynamic constraints. "sat" is constructed unless ( wc & 1 ) ==
  *   true, and "cap" is constructed unless ( wc & 2 ) == true.
  *
  * - If the "knapsack formulation" (KF) is used, then there is only one
  *   explicit groups of "linear constraints", the "sat" one with a
  *   std::vector< FRowConstraint > of size f_n_facilities imposing the
  *   satisfaction of customers' demands (1), which is constructed unless
  *   ( wc & 1 ) == true; the capacity constraints are inside the
  *   BinaryKnapsackBlock sub-Block, and the corresponding constraints are
  *   constructed unless ( wc & 2 ) == true.
  *
  * - If the "flow formulation" (FF) is used, then there is only one
  *   explicit groups of "linear constraints", the "cap" one with a
  *   std::vector< FRowConstraint > of size f_n_customers imposing the
  *   linking between the Y[] variables in the first sub-Block and the
  *   (appropriate) arc flow variables in the MCFBlock sub-Block; this
  *   is constructed unless ( wc & 2 ) == true, while the constraints in
  *   the MCFBlock sub-Block (which impose the satisfaction of customers'
  *   demands, although they also are a part of the capacity ones) are
  *   constructed unless ( wc & 1 ) == true.
  *
  * Finally, if the third bit of ws is 1, then an appropriately arranged
  * group of dynamic Constraint is added that support the separation of
  * the "strong linking" constraints x_{ij} \leq y_i. If this is done,
  * then generate_dynamic_constraints() implements this separation. Note
  * that, whatever the formulation, che corresponding "strong" group of
  * dynamic Constraint is always added to the "root"
  * CapacitatedFacilityLocationBlock */
 
 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the dynamic Constraints of the CapacitatedFacilityLocationBlock
 /** The CapacitatedFacilityLocationBlock (in whatever formulation) only has
  * a single group of dynamic Constraint: the "strong forcing" constraints
  * x_{ij} <= y_i. These can only be dynamically separated if "declared" as
  * being part of the formulation in generate_abstract_constraints() (see
  * the comments there).
  *
  * If not nullptr, dycc must be a SimpleConfiguration< std::pair< int ,
  * double > > *. The first parameter is the maximum number of constraints
  * to generate at each separation round; a negative number (default) means
  * "infinite". If the number of constraints to be inserted exceeds the
  * number of those that could be separated, the ones with larger violations
  * are preferred. The second parameter (default 1e-4) is the minimal absolute
  * violation required to declare a constraint violated (since both variables
  * have no coefficients and are in [ 0 , 1 ] an absolute threshold, for once,
  * makes sense). */

 void generate_dynamic_constraints( Configuration * dycc = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the CFL
 /** Method that generates the objective of the CFL. The actual form of the
  * objective obviously depends on the variables, and therefore on the
  * formulation which is decided during the call to
  * generate_abstract_variables(), as follows:
  *
  * - If the "natural formulation" (SF) is used, there is a single "dense"
  *   "linear objective" (FRealObjective with a LinearFunction) having
  *   first the terms F[ i ] Y[ i ] in order of i, and then all the terms
  *   C[ i , j ] X[ i , j ] in order of i and then of j.
  *
  * - If the "knapsack formulation" (KF) is used, then all the objective is
  *   expressed in terms of the objectives of the f_n_facilities
  *   BinaryKnapsackBlock sub-Block.
  *
  * - If the "flow formulation" (FF) is used, then there is a single "dense"
  *   "linear objective" (FRealObjective with a LinearFunction) having
  *   the terms F[ i ] Y[ i ] (in order of i) in the first sub-Block where
  *   the Y[] variables are defined; the remaining part of the objective
  *   is represented by the "linear objective" of the MCFBlock, that has
  *   cost C[ i , j ] / D[ j ] (since the flow on the arc represent the
  *   actual amount of commodity shipped along the arc, as opposed to the
  *   fraction of the demand D[ j ]) in the "transportation arc" ( i , j )
  *   (from facility i to customer j), but zero costs on the "facility arcs"
  *   (from the super-source to facilities).
  *
  * There are no options (save those above relative to the employed
  * formulation which must have already been set before the method is
  * called), hence objc is ignored. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-- Methods for reading the data of the CapacitatedFacilityLocationBlock --*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the CapacitatedFacilityLocationBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// getting the current sense of the Objective, which is always minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }
  
/*--------------------------------------------------------------------------*/
 /// getting upper bounds on the value of the Objective
 /** An upper bound on the optimal value of the problem is computed as
  * \f$ \sum_{ i \in I } : F[ i ] > 0 F[ i ] \f$ plus, for each customer
  * j \in J, the term \f$ C[ i , j ] \f$ corresponding to the maximum
  * \f$ C[ i , j ] \f$ among all possible i \in I. */

 [[nodiscard]] double get_valid_upper_bound( bool conditional = false )
  override final {
  if( ! conditional )
   return( Inf< double >() );
   
  if( std::isnan( f_cond_upper ) )
   compute_conditional_bounds();

  return( f_cond_upper );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// getting a global valid lower bound on the value of the Objective
 /** An upper bound on the optimal value of the problem is computed as
  * \f$ \sum_{ i \in I } : F[ i ] < 0 F[ i ] \f$ plus, for each customer
  * j \in J, the term \f$ C[ i , j ] \f$ corresponding to the minimum
  * \f$ C[ i , j ] \f$ among all possible i \in I. */

 [[nodiscard]] double get_valid_lower_bound( bool conditional = false )
  override final {
  if( std::isnan( f_cond_lower ) )
   compute_conditional_bounds();

  return( f_cond_lower );
  }

/*--------------------------------------------------------------------------*/
 /// get the number of facilities

 [[nodiscard]] Index get_NFacilities( void ) const {
  return( f_n_facilities );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of customers

 [[nodiscard]] Index get_NCustomers( void ) const {
  return( f_n_customers );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the splittable/unsplittable status

 [[nodiscard]] bool get_UnSplittable( void ) const {
  return( f_unsplittable );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of facility capacities

 [[nodiscard]] c_DVector & get_Capacities( void ) const {
  return( v_capacity );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the capacity of facility i (0 <= i < get_NFacilities())

 [[nodiscard]] Demand get_Capacity( Index i ) const {
  return( v_capacity[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of facility fixed costs

 [[nodiscard]] c_CVector & get_Fixed_Costs( void ) const {
  return( v_f_cost );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the fixed cost of facility i (0 <= i < get_NFacilities())

 [[nodiscard]] Cost get_Fixed_Cost( Index i ) const {
  return( v_f_cost[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of customers' demands

 [[nodiscard]] c_DVector & get_Demands( void ) const {
  return( v_demand );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the demands of customer i (0 <= i < get_NCustomers())

 [[nodiscard]] Demand get_Demand( Index i ) const {
  return( v_demand[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the matrix of transportation costs
 /** Returns a const reference to the matrix of transportation costs, i.e.,
  * element [ i ][ j ] is the *unitary* transportation cost between facility
  * i and customer j. */

 [[nodiscard]] c_CMatrix & get_Transportation_Costs( void ) const {
  return( v_t_cost );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the transportation cost for a given pair ( facility , customer )

 [[nodiscard]] Cost get_Transportation_Cost( Index facility , Index customer )
  const {
  return( v_t_cost[ facility ][ customer ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i return true if the corresponding facility is fixed

 [[nodiscard]] bool is_fixed( Index i ) const {
  return( v_fxd[ i ] != 0 );  // 0 == yFree
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i return true if the corresponding facility is closed

 [[nodiscard]] bool is_closed( Index i ) const {
  return( v_fxd[ i ] == 1 );  // 1 == yFxd0
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// given an index i return true if the corresponding facility is fixed open

 [[nodiscard]] bool is_fixed_open( Index i ) const {
  return( v_fxd[ i ] == 2 );  // 2 == yFxd1
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the design
  * (Y) and transportation (X) Variable of the
  * CapacitatedFacilityLocationBlock is approximately feasible. This clearly
  * requires the Variable of the CapacitatedFacilityLocationBlock to have
  * been defined, i.e., that generate_abstract_variables() has been called
  * prior to this method, and then of course that a solution has been
  * written there.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single double value, representing the *relative* tolerance for
  * satisfaction of both the customer satisfaction constraint and the
  * facility capacity ones. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< double >, then
  *   it if fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< double >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 1e-10. */

 bool is_feasible( bool useabstract = false ,
		   Configuration * fsbc = nullptr ) override;

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) customer feasible
 /** Returns true if the solution encoded in the current value of the design
  * (Y) and transportation (X) Variable of the
  * CapacitatedFacilityLocationBlock (approximately) satisfies the customer
  * satisfaction constraints. This clearly requires the Variable to have been
  * defined, i.e., that generate_abstract_variables() has been called prior
  * to this method, and then of course that a solution has been written
  * there. The parameter feps is the relative accuracy defining
  * "approximately". The parameter "useabstract" has the same meaning as in
  * is_feasible(). */

 bool customer_feasible( double eps , bool useabstract = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the current solution is (approximately) facility feasible
 /** Returns true if the solution encoded in the current value of the design
  * (Y) and transportation (X) Variable of the
  * CapacitatedFacilityLocationBlock (approximately) satisfies the facility
  * capacity constraints. This clearly requires the Variable to have been
  * defined, i.e., that generate_abstract_variables() has been called prior
  * to this method, and then of course that a solution has been written
  * there. The parameter feps is the relative accuracy defining
  * "approximately". The parameter "useabstract" has the same meaning as in
  * is_feasible(). */

 bool facility_feasible( double eps , bool useabstract = false );

/** @} ---------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
 *  @{ */

 /// gets an R3 Block of CapacitatedFacilityLocationBlock 
 /** Gets an R3 Block of the CapacitatedFacilityLocationBlock. The specific
  * type of R3 Block is encoded in r3bc, which (if not nullptr) must be a
  * SimpleConfiguration< int >, whose f_value can be:
  *
  * - 0 (which is also assumed if r3bc == nullptr), indicating the standard
  *   copy (a CapacitatedFacilityLocationBlock identical to this); in this
  *   case, if \p base is passed then it must be a
  *   CapacitatedFacilityLocationBlock (otherwise exception is thrown)
  *
  * - 1 and 2: a MCFBlock containing the "flow relaxation" of the
  *   CapacitatedFacilityLocationBlock is produced, the difference between
  *   the two values is that with 2 the "flow relaxation" is arranged in
  *   such a way as to always being feasible even of the "normal" arcs of
  *   the graph are closed (or their capacity is reduced). If \p base is
  *   passed then it must be a MCFBlock (otherwise exception is thrown).
  *
  *   In both cases, the graph has f_n_facilities + f_n_customers + 1 nodes,
  *   with the following arrangement:
  *
  *   = nodes 1 ... f_n_facilities: facilities, deficit == 0
  *
  *   = nodes f_n_facilities + 1 ... f_n_facilities + f_n_customers: 
  *     customers, deficit == customer demand
  *
  *   = node f_n_facilities + f_n_customers + 1: super source,
  *     deficit == - sum of all demands
  *
  *   and *at least* f_n_facilities * ( f_n_customers + 1 ) with the
  *   following arrangement:
  *
  *   = arcs 0 ...f_n_facilities - 1: from source to facilities, capacity ==
  *     facility capacity, cost == fixed cost / facility capacity
  *
  *   = arcs f_n_facilities ... f_n_facilities * ( f_n_customers + 1 ) - 1:
  *     from facilities to customers, arranged facility-wise (first all the
  *     arcs of the first facilities, then all the arcs of the second, ...),
  *     capacity == infinite (Inf< MCFBlock::FNumber >()), cost ==
  *     unitary transportation cost between facility and customer
  *
  *   When the value is 2, f_n_customers "artificial" arcs are also added
  *   from the super source to each customer, with capacity == infinite
  *   (Inf< (MCFBlock::FNumber >()) and a very large cost (something like
  *   100 * ( max facility cost + max transportation cost from any
  *   facility to the customer). Note that
  *
  *       THE COST OF THESE ARCS ARE NOT UPDATED WHEN THE ORIGINAL COSTS
  *       CHANGE
  *
  *   as it is assumed that they are "large enough for good".
  *
  * - any other value: not supported (exception is thrown) */

 [[nodiscard]] Block * get_R3_Block( Configuration * r3bc = nullptr ,
				     Block * base = nullptr ,
				     Block * father = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// maps back the solution from a R3Block
 /** Takes the solution stored in R3B, which is supposed to be a R3Block of
  * the type indicated by the [Simple]Configuration[< int >] r3bc with the
  * same encoding as in get_R3_Block(), and moves it in the current Block.
  *
  * If not nullptr, solc is assumed to be a SimpleConfiguration< int > whose
  * f_value encodes bit-wise which part of the solution is mapped back:
  *
  * - bit 0 (+1): the design solution (y) is mapped back
  *
  * - bit 1 (+2): the transportation solution (x) is mapped back
  *
  * If solc == nullptr, the value of 3 (map back everything) is assumed. */

 void map_back_solution( Block * R3B , Configuration * r3bc = nullptr ,
			 Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps forward the solution to a R3Block
 /** Takes the solution stored in the current Block and moves it in R3B,
  * which is supposed to be a R3Block of the type indicated by the
  * [Simple]Configuration[< int >] r3bc with the same encoding as in
  * get_R3_Block(), 
  *
  * If not nullptr, solc is assumed to be a SimpleConfiguration< int > whose
  * f_value encodes bit-wise which part of the solution is mapped forward:
  *
  * - bit 0 (+1): the design solution (y) is mapped forward
  *
  * - bit 1 (+2): the transportation solution (x) is mapped forward
  *
  * If solc == nullptr, the value of 3 (map forward everything) is assumed.
  */

 void map_forward_solution( Block * R3B , Configuration * r3bc = nullptr ,
			    Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// maps forward Modification to a R3Block
 /** Maps forward Modification to a R3Block; the r3bc has the same meaning as
  * in get_R3_Block() and specifies what kind of R3Block R3B is.
  *
  * IMPORTANT NOTE: map_forward_Modification() only maps "physical"
  * Modification. The point is that if any part of the "abstract
  * representation" of CapacitatedFacilityLocationBlock is changed, the
  * corresponding "abstract" Modification is intercepted in add_Modification()
  * and a "physical" Modification is also issued. Hence, for any change in
  * CapacitatedFacilityLocationBlock there will always be both Modification
  * "in flight", and therefore there is no need (and good reasons not) to map
  * both.
  *
  * In particular, the method handles the following Modification:
  *
  * - GroupModification
  *
  * - CapacitatedFacilityLocationBlockRngdMod
  *
  * - CapacitatedFacilityLocationBlockSbstMod
  *
  * - NBModification
  *
  * Any other Modification is ignored (and false is returned).
  *
  * Note that for GroupModification, true is returned only if all the
  * inner Modification of the GroupModification return true.
  *
  * Note that if the issueAMod param is eModBlck, then it is "downgraded" to
  * eNoBlck: the method directly does "physical" changes, hence there is no
  * reason for it to issue "abstract" Modification with concerns_Block() ==
  * true. */

 bool map_forward_Modification( Block * R3B , c_p_Mod mod ,
				Configuration * r3bc = nullptr ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eModBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps back Modification from a R3Block
 /** Maps back Modification from a R3Block; the r3bc has the same meaning as
  * in get_R3_Block() and specifies what kind of R3Block R3B is.
  *
  * Currently, this is only implemented for the "copy" R3Block (dur to the
  * fantastically dirty trick whereby, since the two objects are copies,
  * mapping back a Modification to this from R3B is the same as mapping
  * forward a Modification from R3B to this); it will always return false
  * otherwise. */

 bool map_back_Modification( Block * R3B , c_p_Mod mod ,
			     Configuration * r3bc = nullptr ,
			     ModParam issuePMod = eNoBlck ,
			     ModParam issueAMod = eModBlck ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns the current CapacitatedFacilityLocationBlockSolution
 /** Returns a CapacitatedFacilityLocationBlockSolution representing the
  * current solution status of this CapacitatedFacilityLocationBlock. What
  * part of the solution is saved depends on the integer value ws, obtained
  * as follows:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then
  *   ws == solc->f_value:
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, ws is its f_value
  *
  * - otherwise ws == 0.
  *
  * The encoding of ws is:
  *
  *   = 1 means "only save the design part of the solution (y)"
  *
  *   = 2 means "only save the transportation part of the solution (x)"
  *
  *   = everything else (e.g., 0) means "save everything";
  *
  * Note that CapacitatedFacilityLocationBlockSolution may not contain the
  * required solution if the Variable have not been constructed yet: this
  * throws an exception, unless emptys == true, in which case the
  * CapacitatedFacilityLocationBlockSolution object is only prepped for
  * getting a solution, but it is not really getting one now.
  *
  * Note that, although the method clearly returns a
  * CapacitatedFacilityLocationBlockSolution, formally the return type is
  * Solution *. This is because it is not possible to forward declare
  * CapacitatedFacilityLocationBlockSolution as a derived class from
  * Solution, nor to define CapacitatedFacilityLocationBlockSolution before
  * CapacitatedFacilityLocationBlock because the former uses some type
  * information declared in the latter. */ 

 [[nodiscard]]  Solution * get_Solution( Configuration * solc = nullptr ,
					 bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// gets (a pointer to) the y[] variable corresponding to facility i
 /** Gets a pointer to the ColVariable representing whether or not facility
  * i is constructed; returns nullptr if the abstract representation has not
  * been constructed (yet). */

 [[nodiscard]] ColVariable * get_y( Index i ) const {
  if( ! ( AR & 8 ) )  // 8 == HasVar
   return( nullptr );

  if( i >= f_n_facilities )
   throw( std::logic_error( "get_y: invalid facility index" ) );

  if( ( AR & 3 ) == 1 )  // 3 = FormMsk , 1 = KskForm
   return( static_cast< BinaryKnapsackBlock * >(
				   v_Block[ i ] )->get_Var( f_n_customers ) );

  // note the need for the const_cast as all fields of the class are const
  // inside of a const method (this is const)
  return( const_cast< ColVariable * >( & v_y[ i ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// gets (a pointer to) the x[ i ][ j ] variable
 /** Gets a pointer to the ColVariable representing how much of the demand
  * of the customer j is served by facility i; returns nullptr if the
  * abstract representation has not been constructed (yet).
  *
  * Important note: in the Standard and Knapsack Formulation the variable
  *                 is in [ 0 , 1 ] and represents the fraction of the
  *                 overall demand of j served by i, but in the Flow
  *                 Formulation the variable actually contains the total
  *                 amount of demand of j served by i, which menas that to
  *                 get the actual value in [ 0 , 1 ] one must divide its
  *                 value by the demand of j. */

 [[nodiscard]] ColVariable * get_x( Index i , Index j ) const {
  if( ! ( AR & 8 ) )  // 8 == HasVar
   return( nullptr );

  if( i >= f_n_facilities )
   throw( std::logic_error( "get_x: invalid facility index" ) );
  if( j >= f_n_customers )
   throw( std::logic_error( "get_x: invalid customer index" ) );

  if( ( AR & 3 ) == 0 )  // 3 = FormMsk, 0 = StdForm
   return( const_cast< ColVariable * >( & v_x[ i ][ j ] ) );
   // note the need for the const_cast as all fields of the class are const
   // inside of a const method like this one is

  if( ( AR & 3 ) == 1 )  // 3 = FormMsk , 1 = KskForm
   return( static_cast< BinaryKnapsackBlock * >(
				              v_Block[ i ] )->get_Var( j ) );

  return( static_cast< MCFBlock * >( v_Block[ 1 ] )->i2p_x(
			          f_n_facilities + i * f_n_customers + j ) );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the (integer) facility solution
 /** Method to get the facility solution: the components of the facility
  * solution vector in the range [ rng.first , rng.second ) are written in
  * the IntSolution in the positions starting from where the iterator Sol
  * points, in the obvious order. Since Sol is an iterator to an IntSolution,
  * this is the integer (binary) version of the solution. */

 void get_facility_solution( IS_it Sol , Range rng = INFRange ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets an arbitrary subset of the (integer) facility solution
 /** Method to get the facility solution: the components of the facility
  * solution vector whose indices are specified in nms[] are written in the
  * IntSolution in the positions starting from where the iterator Sol points,
  * in the obvious order. Since Sol is an iterator to an IntSolution, this is
  * the integer (binary) version of the solution. Note that we don't require
  * nms[] to be ordered, just the entries written in the (sub)vector
  * will be in whatever order nms[] is. */

 void get_facility_solution( IS_it Sol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets a contiguous interval of the (continuous) facility solution
 /** Method to get the facility solution: the components of the facility
  * solution vector in the range [ rng.first , rng.second ) are written in
  * the CntSolution in the positions starting from where the iterator Sol
  * points, in the obvious order. Since Sol is an iterator to a CntSolution,
  * this can be a fractional solution, e.g., as produced by a Solver that can
  * only solve a continuous relaxation of the problem. */

 void get_facility_solution( CS_it Sol , Range rng = INFRange ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets an arbitrary subset of the (continuous) facility solution
 /** Method to get the facility solution: the components of the facility
  * solution vector whose indices are specified in nms[] are written in the
  * IntSolution in the positions starting from where the iterator Sol points,
  * in the obvious order.  Note that we don't require nms[] to be ordered,
  * just the entries written in the (sub)vector will be in whatever order
  * nms[] is. Since Sol is an iterator to a CntSolution, this can be a
  * fractional solution, e.g., as produced by a Solver that can only solve a
  * continuous relaxation of the problem. */

 void get_facility_solution( CS_it Sol , c_Subset & nms ) const;

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the (continuous) transportation solution
 /** Method to get the transportation solution: the components of the
  * transportation solution vector in the range [ rng.first , rng.second )
  * are written in the CntSolution in the positions starting from where the
  * iterator Sol points, in the obvious order. The two-dimensional
  * transportation solution is considered "flattened" into a one-dimensional
  * vector, arranged facility-wise: first the components corresponding to
  * the first facility (in order of customer), then these corresponding to
  * the second facility ...
  *
  * Since Sol is an iterator to a CntSolution, this is the fractional version
  * of the solution, i.e., the "true" solution if the splittable version of
  * the problem is solved (get_Unsplittable() == false), or the integer
  * rounding of the true continuous solution otherwise. */

 void get_transportation_solution( CS_it Sol , Range rng = INFRange ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// gets an arbitrary subset of the (continuous) transportation solution
 /** Method to get the transportation solution: the components of the 
  * transportation solution vector whose indices are specified in nms[] are
  * written in the CntSolution in the positions starting from where the
  * iterator Sol points, in the obvious order. The two-dimensional
  * transportation solution is considered "flattened" into a one-dimensional
  * vector, arranged facility-wise: first the components corresponding to
  * the first facility (in order of customer), then these corresponding to
  * the second facility ... nms[] must be ordered in increasing sense.
  *
  * Since Sol is an iterator to a CntSolution, this is the fractional version
  * of the solution, i.e., the "true" solution if the splittable version of
  * the problem is solved (get_Unsplittable() == false), or the integer
  * rounding of the true continuous solution otherwise. */

 void get_transportation_solution( CS_it Sol , c_Subset & nms ) const;

/*--------------------------------------------------------------------------*/
 /// returns the objective value of the current solution

 [[nodiscard]] RealObjective::OFValue get_objective_value( void );

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the (integer) facility solution
 /** Method to set the facility solution: the components of the facility
  * solution vector in the range [ rng.first , rng.second ) are read from
  * the IntSolution in the positions starting from where the iterator Sol
  * points, in the obvious order. Since Sol is an iterator to an IntSolution,
  * this is the integer (binary) version of the solution. */

 void set_facility_solution( c_IS_it Sol , Range rng = INFRange );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets an arbitrary subset of the (integer) facility solution
 /** Method to set the facility solution: the components of the facility
  * solution vector whose indices are specified in nms[] are read from the
  * IntSolution in the positions starting from where the iterator Sol points,
  * in the obvious order. Since Sol is an iterator to an IntSolution, this is
  * the integer (binary) version of the solution. Note that we don't require
  * nms[] to be ordered, just the entries read from the (sub)vector need be in
  * whatever order nms[] is. */

 void set_facility_solution( c_IS_it Sol , c_Subset & nms );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets a contiguous interval of the (continuous) facility solution
 /** Method to set the facility solution: the components of the facility
  * solution vector in the range [ rng.first , rng.second ) are read from
  * the CntSolution in the positions starting from where the iterator Sol
  * points, in the obvious order. Since Sol is an iterator to a CntSolution,
  * this can be a fractional solution, e.g., as produced by a Solver that can
  * only solve a continuous relaxation of the problem. */

 void set_facility_solution( c_CS_it Sol , Range rng = INFRange );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets an arbitrary subset of the (continuous) facility solution
 /** Method to set the facility solution: the components of the facility
  * solution vector whose indices are specified in nms[] are read from the
  * IntSolution in the positions starting from where the iterator Sol points,
  * in the obvious order. nms[] must be ordered in increasing sense. Since
  * Sol is an iterator to a CntSolution, this can be a  fractional solution,
  * e.g., as produced by a Solver that can only solve a continuous relaxation
  * of the problem. */

 void set_facility_solution( c_CS_it Sol , c_Subset & nms );

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the (continuous) transportation solution
 /** Method to set the transportation solution: the components of the
  * transportation solution vector in the range [ rng.first , rng.second )
  * are read from the CntSolution in the positions starting from where the
  * iterator Sol points, in the obvious order. The two-dimensional
  * transportation solution is considered "flattened" into a one-dimensional
  * vector, arranged facility-wise: first the components corresponding to
  * the first facility (in order of customer), then these corresponding to
  * the second facility ...
  *
  * Since Sol is an iterator to a CntSolution, this is the fractional version
  * of the solution, i.e., the "true" solution if the splittable version of
  * the problem is solved (get_Unsplittable() == false), or the integer
  * rounding of the true continuous solution otherwise. */

 void set_transportation_solution( c_CS_it Sol , Range rng = INFRange );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets an arbitrary subset of the (continuous) transportation solution
 /** Method to set the transportation solution: the components of the 
  * transportation solution vector whose indices are specified in nms[] are
  * read from the CntSolution in the positions starting from where the
  * iterator Sol points, in the obvious order. The two-dimensional
  * transportation solution is considered "flattened" into a one-dimensional
  * vector, arranged facility-wise: first the components corresponding to
  * the first facility (in order of customer), then these corresponding to
  * the second facility ... nms[] must be ordered in increasing sense.
  *
  * Since Sol is an iterator to a CntSolution, this is the fractional version
  * of the solution, i.e., the "true" solution if the splittable version of
  * the problem is solved (get_Unsplittable() == false), or the integer
  * rounding of the true continuous solution otherwise. */

 void set_transportation_solution( c_CS_it Sol , c_Subset & nms );

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// adding a new Modification to the CapacitatedFacilityLocationBlock
 /** Method for handling Modification.
  *
  * The version of CapacitatedFacilityLocationBlock has to intercept any
  * "abstract Modification" that modifies the "abstract representation" of
  * the CapacitatedFacilityLocationBlock, and "translate" them into both
  * changes of the actual data structures and corresponding "physical
  * Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the
  * CapacitatedFacilityLocationBlock, and this must not be repeated. In
  * particular, this would be an issue if the Modification would be
  * [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * CapacitatedFacilityLocationBlock), which would mean that the "physical
  * Modification" would be issued twice.
  *
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - C05FunctionModRngd and C05FunctionModSbst changing coefficients coming
  *   from the (LinearFunction into the FRow)Objective, but *not* from the
  *   (LinearFunction into the FRow)Constraint;
  *
  * - RowConstraintMod changing the RHS of the bound constraints and both
  *   sides at once of the flow conservation ones, but not any other
  *   combination; and note that the RHS of the bound constraints may not
  *   be changeable at all if they have not been constructed, in which
  *   case there cannot be any Modification to handle here;
  *
  * - VariableMod fixing and un-fixing a flow ColVariable; however, note
  *   that *fixing is only permitted if the value() of the ColVariable is
  *   zero*, because that corresponds to closing the arc, exception being
  *   thrown otherwise.
  *
  * Any other Modification reaching the CapacitatedFacilityLocationBlock
  * will lead to exception being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*-------- PRINTING & SAVING THE CapacitatedFacilityLocationBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the CapacitatedFacilityLocationBlock
 *  @{ */

 /// print the CapacitatedFacilityLocationBlock on an ostream
 /** Prints information about the CapacitatedFacilityLocationBlock with the
  * given verbosity level. The "complete" level ('C') outputs the
  * CapacitatedFacilityLocationBlock in the standard ORLib text format,
  * any other level prints only basic information.
  *
  * The verbosity levels corresponding to the other two formats (facility
  * oriented, demands-first and demands-last) are not implemented yet. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * CapacitatedFacilityLocationBlock. See
  * CapacitatedFacilityLocationBlock::deserialize( netCDF::NcGroup & ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the Capacitated Facility Location instance
 *
 * All the methods in this section have two parameters issueMod and issueAMod
 * which control if and how the, respectively, "physical Modification" and
 * "abstract Modification" corresponding to the change have to be issued, and
 * where (to which channel). The format of the parameters is that of
 * Observer::make_par(), except that the value eModBlck is ignored and
 * treated it as if it were eNoBlck [see Observer::issue_pmod()]. This is
 * because it makes no sense to issue an "abstract" Modification with
 * concerns_Block() == true, since the changes in the
 * CapacitatedFacilityLocationBlock have surely been done already, and this
 * is just not possible for a "physical" Modification.
 *
 * Note: the methods accept the eDryRun value for the issueAMod parameter for
 * the "abstract" representation. This allows to re-use them within
 * CapacitatedFacilityLocationBlock itself when reacting to abstract
 * Modification, where the  "abstract" representation has been changed
 * already. However, the eDryRun value is not allowed (it is ignored) for
 * the issuePMod parameter for the "physical" representation, as there is no
 * reasonable use for this. Basically, this makes eDryRun equivalent to
 * eNoMod.
 * @{ */

 /// change the facility_costs of a contiguous interval
 /** Method to change the costs of a subset of facility with "contiguous
  * names". That is, *( NCost + h ) becomes the new cost of the facility
  * rng.first + h, for all 0 <= h < rmg.second - rng.first. Note that any 
  * rng.second >= get_NFacilities() means "up until the end".
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void chg_facility_costs( c_CV_it NCost , Range rng = INFRange ,
                          ModParam issueMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the facility_costs of an arbitrary subset of facilities
 /** Method to change the costs of an arbitrary subset of facilities. That is,
  * *( NCost + h ) becomes the new cost of facility nms[ h ] for all 0 <= h <
  * NCost.size(). \p ordered tells if \p nms is already ordered in increasing
  * sense. As the the && tells, \p nms is "consumed" by the method, typically
  * being shipped to the appropriate CapacitatedFacilityLocationBlockSbstMod
  * that is issued. */

 void chg_facility_costs( c_CV_it NCost ,
                          Subset && nms , bool ordered = false ,
                          ModParam issueMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the cost of the given facility

 void chg_facility_cost( Cost NCost , Index i ,
                         ModParam issueMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the transportation costs of a contiguous interval
 /** Method to change the transportation costs for a subset of the pairs
  * ( facility , customer ) with "contiguous names". The matrix of
  * transportation costs is considered "flattened" into a vector in row-major
  * format, i.e., first all costs for the first facility, then all costs for
  * the second facility, ...; in other words, the pair ( i , j ) is mapped in
  * the element k = i * get_NCustomers() + j of the vector. Given this, 
  * *( NCost + h ) becomes the new  transportation cost of the pair with
  * "name" k = rng.first + h, for all 0 <= h < rmg.second - rng.first. Note
  * that any  rng.second >= get_NFacilities() * get_NCustomers() means "up
  * until the end".
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void chg_transportation_costs( c_CV_it NCost , Range rng = INFRange ,
                                ModParam issueMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the transportation costs of an arbitrary subset
 /** Method to change the transportation costs for an arbitrary subset of the
  * pairs ( facility , customer ). The matrix of transportation costs is
  * considered "flattened" into a vector in row-major format; see the comments
  * to the Range version for details.
  * As the the && tells, \p nms is "consumed" by the method, typically being
  * shipped to the appropriate CapacitatedFacilityLocationBlockSbstMod that
  * is issued. \p order tells if \p nms is already ordered in increasing
  * sense. */

 void chg_transportation_costs( c_CV_it NCost ,
                                Subset && nms , bool ordered = false ,
                                ModParam issueMod = eNoBlck ,
                                ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the transportation of the given pair ( facility , customer )
 /** Method to change the transportation costs for the given pair
  * ( facility , customer ) \p, i.e., facility = p / get_NCustomers() and
  * customer = p % get_NCustomers(). */

 void chg_transportation_cost( Cost NCost , Index p ,
                               ModParam issueMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the capacities of a contiguous interval of facilities
 /** Method to change the capacities of a subset of facility with "contiguous
  * names". That is, *( NCap + h ) becomes the new capacity of the facility
  * rng.first + h, for all 0 <= h < rmg.second - rng.first. Note that any 
  * rng.second >= get_NFacilities() means "up until the end".
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void chg_facility_capacities( c_DV_it NCap , Range rng = INFRange ,
                               ModParam issueMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the capacities of an arbitrary subset of facilities
 /** Method to change the capacities of an arbitrary subset of facilities.
  * That is, *( NCop + h ) becomes the capacity of facility nms[ h ] for all
  * 0 <= h < NCap.size(). \p ordered tells if \p nms is already ordered in
  * increasing sense. As the the && tells, \p nms is "consumed" by the method,
  * typically being shipped to the appropriate
  * CapacitatedFacilityLocationBlockSbstMod that is issued. */

 void chg_facility_capacities( c_DV_it NCap ,
                               Subset && nms , bool ordered = false ,
                               ModParam issueMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the capacity of the given facility

 void chg_facility_capacity( Demand NCap , Index i ,
                             ModParam issueMod = eNoBlck ,
                             ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the demands of a contiguous interval of customers
 /** Method to change the demands of a subset of customers with "contiguous
  * names". That is, *( NDem + h ) becomes the new demand of the customer
  * rng.first + h, for all 0 <= h < rmg.second - rng.first. Note that any 
  * rng.second >= get_NCustomers() means "up until the end".
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void chg_customer_demands( c_DV_it NDem , Range rng = INFRange ,
                            ModParam issueMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// change the demands of an arbitrary subset of customers
 /** Method to change the demands of an arbitrary subset of customers. That
  * is, *( NDem + h ) becomes the demand of customer nms[ h ] for all 0 <= h
  * < NDem.size(). \p ordered tells if \p nms is already ordered in increasing
  * sense. As the the && tells, \p nms is "consumed" by the method, typically
  * being shipped to the appropriate
  * CapacitatedFacilityLocationBlockSbstMod that is issued. */

 void chg_customer_demands( c_DV_it NDem ,
                            Subset && nms , bool ordered = false ,
                            ModParam issueMod = eNoBlck ,
                            ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the demand of the given customer

 void chg_customer_demand( Demand NDem , Index j ,
                           ModParam issueMod = eNoBlck ,
                           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// closes a contiguous interval of facilities
 /** Method to close a subset of facility with "contiguous names", i.e., all
  * facilities rng.first <= i < rmg.second - rng.first. Note that any
  * rng.second >= get_NFacilities() means "up until the end". Closing an
  * already closed facility or a fixed-open one does nothing, i.e., closing
  * a facility does not override the fixed-open status.
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void close_facilities( Range rng = INFRange ,
                        ModParam issueMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes an arbitrary subset of facilities
 /** Method to close an arbitrary subset of facilities, i.e., all those whose
  * names are found in \p nms. Closing an already closed facility or a
  * fixed-open one does nothing, i.e., closing a facility does not override
  * the fixed-open status. \p ordered tells if \p nms is already ordered in
  * increasing sense. As the && tells, \p nms is "consumed" by the method,
  * typically being shipped to an appropriate
  * CapacitatedFacilityLocationBlockSbstMod object. */

 void close_facilities( Subset && nms , bool ordered = false ,
                        ModParam issueMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes the given facility
 /** Closes the given facility, unless i is already closed or fixed-open in
  * which case it one does nothing; that is, closing a facility does not
  * override the fixed-open status. */

 void close_facility( Index i , ModParam issueMod = eNoBlck ,
                      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// re-opens a contiguous interval of facilities
 /** Method to re-open a subset of facility with "contiguous names", i.e.,
  * all facilities rng.first <= i < rmg.second - rng.first that had
  * previously been closed are now open again. Note that any rng.second >=
  * get_NFacilities() means "up until the end". Re-opening a facility that
  * had not been previously closed does nothing, while re-openomg a precently
  * fixed-open facility overrides the fixed-open status.
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void open_facilities( Range rng = INFRange ,
                       ModParam issueMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens an arbitrary subset of facilities
 /** Method to re-opens an arbitrary subset of facilities, i.e., all those
  * whose names are found in \p nms. Re-opening a facility that had not been
  * previously closed does nothing, while re-openomg a precently fixed-open
  * facility overrides the fixed-open status. \p ordered tells if \p nms is
  * already ordered in increasing  sense. As the && tells, \p nms is
  * "consumed" by the method, typically being shipped to an appropriate
  * CapacitatedFacilityLocationBlockSbstMod object. */

 void open_facilities( Subset && nms , bool ordered = false ,
                       ModParam issueMod = eNoBlck ,
                       ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens the given facility

 void open_facility( Index i , ModParam issueMod = eNoBlck ,
                     ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// fix open a contiguous interval of facilities
 /** Method to fix open a subset of facility with "contiguous names", i.e.,
  * all facilities rng.first <= i < rmg.second - rng.first are now considered
  * to be open already. Note that their construction cost is added to the
  * Objective value, but it is not optimised upon since the decision is taken
  * already. Note that any rng.second >= get_NFacilities() means "up until 
  * the end". Fixing open an already fixed-open or a closed facility does
  * nothing, i.e., fixing-open a facility does not override its closed status.
  *
  * If issueMod says so then a "physical"
  * CapacitatedFacilityLocationBlockRngdMod is issued. */

 void fix_open_facilities( Range rng = INFRange ,
			   ModParam issueMod = eNoBlck ,
			   ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// fix open an arbitrary subset of facilities
 /** Method to fix open an arbitrary subset of facilities, i.e., all those
  * whose names are found in \p nms are now considered to be open already.
  * Note that their construction cost is added to the Objective value, but it
  * is not optimised upon since the decision is taken  already. Fixing open
  * an already fixed-open or a closed facility does nothing, i.e., fixing-open
  * a facility does not override its closed status. \p ordered tells if
  * \p nms is already ordered in increasing  sense. As the && tells, \p nms
  * is "consumed" by the method, typically being shipped to an appropriate
  * CapacitatedFacilityLocationBlockSbstMod object. */

 void fix_open_facilities( Subset && nms , bool ordered = false ,
			   ModParam issueMod = eNoBlck ,
			   ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// fix opens the given facility
 /** Fixes open a facility, unless it is already fixed-open or a closed in
  * which case it does nothing; that is, fixing-open a facility does not
  * override its closed status. */ 

 void fix_open_facility( Index i , ModParam issueMod = eNoBlck ,
			           ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the splittable/unsplittable status of the problem
 /** If \p unsplt == true [default], sets the problem to be the unsplittable
  * version, where each customer need be served by exactly one facility;
  * otherwise sets the problem to be the splittable version, where each
  * customer can be served by any number of facilities).
  *
  * The unsplittable version of the problem cannot be represented when
  * using the Flow Formulation (this is inherent since flow variables are
  * scaled), so calling chg_UnSplittable( true  ) will result in an exception
  * been thrown. */

 void chg_UnSplittable( bool unsplt = true ,
                        ModParam issueMod = eNoBlck ,
                        ModParam issueAMod = eNoBlck );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FRIENDS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 // problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index f_n_facilities;  ///< the number of facilities
 Index f_n_customers;   ///< the number of customers

 DVector v_capacity;    ///< vector of facility capacities
 CVector v_f_cost;      ///< vector of facility fixed costs
 DVector v_demand;      ///< vector of customers demands
 CMatrix v_t_cost;      ///< matrix of transportation costs
                        /**< The matrix of transportation costs is arranged
			 * facility-wise, i.e., v_t_cost[ i ][ j ] is the
  * *total* transportation cost between facility i and customer j. */

 std::vector< unsigned char > v_fxd;  ///< vector saying how the y are fixed
                                      /**< f_fxd[ i ] indicates if y_i is
				       * fixed, with the following encoding:
  * 0 = not fixed, 1 = fixed to 0 (closed), 2 = fixed to 1 (fixed_open) */

 // abstract representation stuff - - - - - - - - - - - - - - - - - - - - - -

 unsigned char AR;       ///< bit-wise coded: what abstract is there
                         /**< The char field AR keeps track of which part of
			  * the abstract representation has been constructed
  * already, as well as *which formulation* is used.
  * The second part is coded in the first three bits of AR, as follows.
  * The first part is coded in the first three bits of AR, as follows:
  * The first two bits encode the "large-scale shape" of the formulation:
  *
  * - ( AR & FormMsk ) == StdForm: the "standard" formulation is used
  * - ( AR & FormMsk ) == KskForm: the "knapsack" formulation is used
  * - ( AR & FormMsk ) == FlwForm: the "flow" formulation is used
  *
  * Then, the third bit ( AR & UnSpltF ) is 1 if the formulation is
  * unsplittable (the X[] are integer).
  *
  * The following bits encode which parts of the abstract formulation have
  * been constructed:
  *
  * - AR & HasVar:      if the Variable have been constructed
  * - AR & HasObj:      if the Objective has been constructed
  * - AR & HasSatCns:   if the customer satisfaction Constraints have been
  *                     constructed
  * - AR & HasCapCns:   if the capacity Constraints have been constructed
  * - AR & HasStrngCns: if the strong Linking Constraints have been
  *                     prepared for being separated */

 bool f_unsplittable;    ///< if customers can only be served once

 double f_cond_lower;    ///< conditional lower bound, can be infinite
 double f_cond_upper;    ///< conditional upper bound, can be infinite

 boost::multi_array< ColVariable , 2 > v_x;  ///< the flow variables
                                             /**< x is a bi-dimensional array
				              * of ColVariable representing
  * transportation; that is, v_x[ i ][ j ] is the fraction of demand of
  * customer j served by facility i. */

 std::vector< ColVariable > v_y;  ///< the design variables
                                  /**< y is the vector of ColVariable
				   * representing facility opening. */

 std::vector< FRowConstraint> v_sat;  ///< the customer satisfaction constrs.

 std::vector< FRowConstraint> v_cap;  ///< the facility capacity constraints

 std::vector< std::list< FRowConstraint > > v_sfc;
 ///< the strong forcing constraints
 /**< v_sfc is an array of FRowConstraint for the dynamic separation of
  * "strong forcing" constraints x_{ij} \leq y_i; v_sfc[ i ] are all the
  * strong forcing" constraints corresponding to y_i. */

 FRealObjective f_obj;                ///< the (linear) objective function

 // Modification handling - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool f_mod_skip;  ///< if the Modification is a "abstract-physical" one
                   /**< CapacitatedFacilityLocationBlock uses the "physical"
		    * Modification of the inner Block as "abstract"
  * Modification. This creates an issue whereby when a  "physical"
  * Modification of an inner Block is processed, it is not known whether it
  * has been issued due to a change that the CapacitatedFacilityLocationBlock
  * already knows about, and therefore can (and must) be ignored, or it has 
  * been issued due to a change that the CapacitatedFacilityLocationBlock
  * does not know about, and therefore it must be processed. The solution to
  * this issue is that CapacitatedFacilityLocationBlock ensures that
  * f_mod_skip == true if and only if this is what is happening. That is,
  * f_mod_skip is false by default, it is set to true right before calling
  * the chg_* methods of the sub-Block (within which "physical" Modification
  * are issued and CapacitatedFacilityLocationBlock::add_Modification() is
  * called) and put back to false immediately after this happened. This would
  * be dangerous in case multiple changes would be happening at the same time,
  * with some of them "known already" by CapacitatedFacilityLocationBlock and
  * others not, but this is not supposed to happen since the mechanism is only
  * activated inside the chg_* methods of CapacitatedFacilityLocationBlock,
  * which are only supposed to be called when it is lock()-ed. */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// register CapacitatedFacilityLocationBlock methods in the method factories
/** Although in general private methods should not be commented, this one is
 * because it does the registration of the following
 * CapacitatedFacilityLocationBlock methods:
 *
 * - chg_facility_costs() (both range and subset version)
 *
 * - chg_transportation_costs() (both range and subset version)
 *
 * - chg_facility_capacities() (both range and subset version)
 *
 * - chg_customers_demands() (both range and subset version)
 *
 * - close_facilities() (both range and subset version)
 *
 * - open_facilities() (both range and subset version)
 *
 * into the corresponding method factories.
 */

 static void static_initialization( void )
 {
  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Range >(
   "CapacitatedFacilityLocationBlock::chg_facility_costs",
   & CapacitatedFacilityLocationBlock::chg_facility_costs );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Subset && ,
   bool >(
   "CapacitatedFacilityLocationBlock::chg_facility_costs" ,
   & CapacitatedFacilityLocationBlock::chg_facility_costs );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Range >(
   "CapacitatedFacilityLocationBlock::chg_transportation_costs" ,
   &CapacitatedFacilityLocationBlock::chg_transportation_costs );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Subset && ,
   bool >(
   "CapacitatedFacilityLocationBlock::chg_facility_capacities",
   & CapacitatedFacilityLocationBlock::chg_facility_capacities );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Range >(
   "CapacitatedFacilityLocationBlock::chg_facility_capacities",
   & CapacitatedFacilityLocationBlock::chg_facility_capacities );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Subset && ,
   bool >(
   "CapacitatedFacilityLocationBlock::chg_customer_demands",
   &CapacitatedFacilityLocationBlock::chg_customer_demands );

  register_method< CapacitatedFacilityLocationBlock , MF_dbl_it , Range >(
   "CapacitatedFacilityLocationBlock::chg_customer_demands",
   & CapacitatedFacilityLocationBlock::chg_customer_demands );

   register_method< CapacitatedFacilityLocationBlock , Range >(
   "CapacitatedFacilityLocationBlock::close_facilities",
   & CapacitatedFacilityLocationBlock::close_facilities );

  register_method< CapacitatedFacilityLocationBlock , Subset && , bool >(
   "CapacitatedFacilityLocationBlock::close_facilities" ,
   & CapacitatedFacilityLocationBlock::open_facilities );

  register_method< CapacitatedFacilityLocationBlock , Range >(
   "CapacitatedFacilityLocationBlock::open_facilities",
   & CapacitatedFacilityLocationBlock::open_facilities );

  register_method< CapacitatedFacilityLocationBlock , Subset && , bool >(
   "CapacitatedFacilityLocationBlock::open_facilities" ,
   & CapacitatedFacilityLocationBlock::open_facilities );

  register_method< CapacitatedFacilityLocationBlock , Range >(
   "CapacitatedFacilityLocationBlock::fix_open_facilities",
   & CapacitatedFacilityLocationBlock::fix_open_facilities );

  register_method< CapacitatedFacilityLocationBlock , Subset && , bool >(
   "CapacitatedFacilityLocationBlock::fix_open_facilities" ,
   & CapacitatedFacilityLocationBlock::fix_open_facilities );
   }

/*--------------------------------------------------------------------------*/

 LinearFunction * get_lfo( void ) {
  #ifdef NDEBUG
   return( static_cast< LinearFunction * >( f_obj.get_function() ) );
  #else
   auto lfo = dynamic_cast< LinearFunction * >( f_obj.get_function() );
   assert( lfo );
   return( lfo );
  #endif
  }

/*--------------------------------------------------------------------------*/
 
 void guts_of_destructor( void );

 void guts_of_get_R3B_MCF( MCFBlock * mcfb , int wR3B );

 void guts_of_chg_tcost_MCF( MCFBlock * mcfb , Range rng ,
			     ModParam issueMod , ModParam issueAMod );

 void guts_of_chg_tcost_MCF( MCFBlock * mcfb , c_Subset & nms , bool ordered ,
			     ModParam issueMod , ModParam issueAMod );
 
 void guts_of_chg_dem_MCF( MCFBlock * mcfb , Range rng ,
			   ModParam issueMod , ModParam issueAMod );

 void guts_of_chg_dem_MCF( MCFBlock * mcfb , c_Subset & nms , bool ordered ,
			   ModParam issueMod , ModParam issueAMod );

 void guts_of_add_ModificationSF( c_p_Mod mod , ChnlName chnl );

 void guts_of_add_ModificationKFG( const GroupModification * mod ,
				   ChnlName chnl );

 void guts_of_add_ModificationKFP( const BinaryKnapsackBlockMod * mod ,
				   ChnlName chnl );

 void guts_of_add_ModificationFFA( c_p_Mod mod , ChnlName chnl );

 void guts_of_add_ModificationFFG( const GroupModification * mod ,
				   ChnlName chnl );

 void guts_of_add_ModificationFFP( const MCFBlockMod * mod , ChnlName chnl );

 bool guts_of_map_f_Mod_copy(
			CapacitatedFacilityLocationBlock * R3B , c_p_Mod mod ,
			ModParam issuePMod , ModParam issueAMod );

 bool guts_of_guts_of_map_f_Mod_copy(
			CapacitatedFacilityLocationBlock * R3B , c_p_Mod mod ,
			ModParam issuePMod , ModParam issueAMod );

 bool guts_of_map_f_Mod_MCF( MCFBlock * R3B , c_p_Mod mod ,
			     ModParam issuePMod , ModParam issueAMod );

 bool guts_of_guts_of_map_f_Mod_MCF( MCFBlock * R3B , c_p_Mod mod ,
				    ModParam issuePMod , ModParam issueAMod );

 void compute_conditional_bounds( void );

 template< class T >
 void get_y( typename std::vector< T >::iterator Sol , Range rng ) const;

 template< typename T >
 void get_y( typename std::vector< T >::iterator Sol , c_Subset nms ) const;

 template< typename T >
 void get_x( typename std::vector< T >::iterator Sol , Range rng ) const;

 template< typename T >
 void get_x( typename std::vector< T >::iterator Sol , c_Subset nms ) const;

 template< typename T >
 void set_y( typename std::vector< T >::const_iterator Sol , Range rng );

 template< typename T >
 void set_y( typename std::vector< T >::const_iterator Sol , c_Subset nms );

 template< typename T >
 void set_x( typename std::vector< T >::const_iterator Sol , Range rng );

 template< typename T >
 void set_x( typename std::vector< T >::const_iterator Sol , c_Subset nms );
 
 ModParam make_amod_param( ModParam issueAMod , Index num );

 void unmake_amod_param( ModParam oldiAM , ModParam newiAM , Index num );

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

 void CheckAbsVSPhys( void );
 
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;
 // insert CapacitatedFacilityLocationBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( CapacitatedFacilityLocationBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------- CLASS CapacitatedFacilityLocationBlockMod -----------------*/
/*--------------------------------------------------------------------------*/
/// Modification for changes to a CapacitatedFacilityLocationBlock
/** Derived class from Modification to describe changes to a
 * CapacitatedFacilityLocationBlock. This is "almost abstract", since it is
 * only directly used to communicate changes "without data", i.e., the
 * problem type to be set to splittable/unsplittable, while most of the
 * changes involve specifying which (part of the) data has changed, and this
 * is demanded to derived classes (which do this in different ways). Note
 * that it is derived from Modification rather than, say, BlockMod (which has
 * the same structure) because this is a class of "physical Modification".
 * This means that a CapacitatedFacilityLocationBlockMod refers to changes in
 * the "physical representation" of the CapacitatedFacilityLocationBlock; the
 * corresponding changes in the "abstract representation" (if any) are dealt
 * with by means of "abstract Modification", i.e., derived classes from
 * AModification (as is BlockMod, which is why
 * CapacitatedFacilityLocationBlockMod is not derived from BlockMod). */

class CapacitatedFacilityLocationBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of CapacitatedFacilityLocationBlockMod
 
 enum MCFB_mod_type {
  eChgFCost = 0 ,   ///< change the facility (design) costs
  eChgTCost     ,   ///< change the transportation costs
  eChgCap       ,   ///< change the facility capacities
  eChgDem       ,   ///< change the customers demands
  eCloseF       ,   ///< close facilities
  eOpenF        ,   ///< re-open facilities
  eBuyF         ,   ///< fix open facilities
  eChgUnSplt    ,   ///< change problem type to unsplittable
  eChgSplt          ///< change problem type to splittable
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the CapacitatedFacilityLocationBlock and the type

 CapacitatedFacilityLocationBlockMod(
		       CapacitatedFacilityLocationBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~CapacitatedFacilityLocationBlockMod() = default;
 ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /** returns the [CapacitatedFacilityLocationBlock]Block to which the
  * CapacitatedFacilityLocationBlockMod refers */

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the CapacitatedFacilityLocationBlockMod

 void print( std::ostream &output ) const override {
  output << "CapacitatedFacilityLocationBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgFCost ): output << "change the facility (design) costs "; break;
   case( eChgTCost ): output << "change the transportation costs "; break;
   case( eChgCap ):   output << "change the facility capacities "; break;
   case( eChgDem ):   output << "change the customers demands "; break;
   case( eCloseF ):   output << "close facilities "; break;
   case( eOpenF ):    output << "re-open facilities "; break;
   default:           output << "fix open facilities ";
   }
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 CapacitatedFacilityLocationBlock *f_Block;
 /**< pointer to the CapacitatedFacilityLocationBlock to which the
  * CapacitatedFacilityLocationBlockMod refers */

 int f_type;   ///< type of modification

/*--------------------------------------------------------------------------*/

 };  // end( class( CapacitatedFacilityLocationBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------ CLASS CapacitatedFacilityLocationBlockRngdMod ---------------*/
/*--------------------------------------------------------------------------*/
/// derived from CapacitatedFacilityLocationBlockMod for ranged modifications
/** Derived class from CapacitatedFacilityLocationBlockMod to describe ranged
 * modifications to a CapacitatedFacilityLocationBlock, i.e., modifications
 * that apply to a Range of either the facilities or the customers. */

class CapacitatedFacilityLocationBlockRngdMod
 : public CapacitatedFacilityLocationBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the CFLnBlock, the type, and the range

 CapacitatedFacilityLocationBlockRngdMod(
				  CapacitatedFacilityLocationBlock * fblock ,
				  int type , Block::Range rng )
  : CapacitatedFacilityLocationBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~CapacitatedFacilityLocationBlockRngdMod() = default;
 ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the CapacitatedFacilityLocationBlockRngdMod

 void print( std::ostream &output ) const override {
  CapacitatedFacilityLocationBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( CapacitatedFacilityLocationBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------- CLASS CapacitatedFacilityLocationBlockSbstMod -----------------*/
/*--------------------------------------------------------------------------*/
/// derived from CapacitatedFacilityLocationBlockMod for subset modifications
/** Derived class from Modification to describe "subset" modifications to a
 * CapacitatedFacilityLocationBlock, i.e., modifications that apply to an
 * arbitrary Subset of either the facilities or the customers. */

class CapacitatedFacilityLocationBlockSbstMod
 : public CapacitatedFacilityLocationBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:


/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the Block *, the type, and the subset
 /**< Constructor: takes the CapacitatedFacilityLocationBlock *, the type,
  * and the subset. As the the && tells, nms is "consumed" by the constructor
  * and its resources become property of the
  * CapacitatedFacilityLocationBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 CapacitatedFacilityLocationBlockSbstMod(
				 CapacitatedFacilityLocationBlock * fblock ,
				 int type , Block::Subset && nms )
  : CapacitatedFacilityLocationBlockMod( fblock , type ) ,
    f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~CapacitatedFacilityLocationBlockSbstMod() = default;
 ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the CapacitatedFacilityLocationBlockSbstMod

 void print( std::ostream &output ) const override {
  CapacitatedFacilityLocationBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( CapacitatedFacilityLocationBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*------------- CLASS CapacitatedFacilityLocationSolution ------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a CapacitatedFacilityLocationBlock
/** The CapacitatedFacilityLocationSolution class, derived from Solution,
 * represents a solution of a CapacitatedFacilityLocationBlock, i.e.:
 *
 * - an m-vector of double for the facility solution
 *
 * - an  (m * n)-vector of double for the transportation solution
 *
 * where m is the number of facilities and n is the number of customers. */

class CapacitatedFacilityLocationSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend CapacitatedFacilityLocationBlock;
 ///< make CapacitatedFacilityLocationBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING MCFSolution ----------------*/

 explicit CapacitatedFacilityLocationSolution( void ) { }
 /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~CapacitatedFacilityLocationSolution() = default;
 ///< destructor: it is virtual, and empty

/* METHODS DESCRIBING THE BEHAVIOR OF A CapacitatedFacilityLocationSolution */

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a CapacitatedFacilityLocationSolution into a netCDF::NcGroup
 /** Serialize a CapacitatedFacilityLocationSolution into a netCDF::NcGroup,
  * with the following format:
  *
  * - The dimension "NFacilities" containing the number of facilities. The
  *   dimension is optional and it is only necessary if the facility
  *   solution is present.
  *
  * - The dimension "TransportationDim" containing the total number of
  *   elements in a transportation solution ( facilities x customers ).
  *   The dimension is optional and it is only necessary if the
  *   transportation solution is present.
  *
  * - The variable "FacilitySolution", of type double and indexed over the
  *   dimension NFacilities. The variable is optional, if it is not specified
  *   then the CapacitatedFacilityLocationSolution object does not contain
  *   any facility solution.
  *
  * - The variable "TransportationSolution", of type double and indexed over
  *   both the dimension NFacilities and NCustomers. The variable is
  *   optional, if it is not specified then the
  *   CapacitatedFacilityLocationSolution object does not contain any
  *   transportation solution. */
 
 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 CapacitatedFacilityLocationSolution * scale( double factor )
  const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 CapacitatedFacilityLocationSolution * clone( bool empty = false )
  const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
  output << "CapacitatedFacilityLocationSolution [" << this << "]: ";
  if( ! v_y.empty() )
   output  << "F";
  if( ! v_x.empty() )
   output  << "T";
  output << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 CapacitatedFacilityLocationBlock::CntSolution v_y;
 ///< the facility solution

 CapacitatedFacilityLocationBlock::CntSolution v_x;
 ///< the transportation solution

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;
 // insert CapacitatedFacilityLocationSolution in the Solution factory
 
 };  // end( class( CapacitatedFacilityLocationSolution ) )

/** @} end( group( CapacitatedFacilityLocationBlock_CLASSES ) ) ------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* CapacitatedFacilityLocationBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------- End File CapacitatedFacilityLocationBlock.h ----------------*/
/*--------------------------------------------------------------------------*/
