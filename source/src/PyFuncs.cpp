// ============================================================================
// Include files 
// ============================================================================
// STD& STL 
// ============================================================================
#include <limits>
#include <climits>
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyFuncs.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// local
// ============================================================================
#include "CallPython.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from namespace Ostap::Functions 
 *  @date   2018-03-31
 *  @author Vanya Belyaev
 */
// ============================================================================
namespace 
{
  // ===========================================================================
  static_assert ( std::numeric_limits<float>::is_specialized , 
                  "std::numeric_limits<float> is not specialized" ) ;
  // ==========================================================================
  constexpr double s_max =   std::numeric_limits<float>::max() ;  
  constexpr double s_min = - std::numeric_limits<float>::max() ;
  // ==========================================================================
  static_assert ( s_max > 0 , "std::numeric_limits<float>::max is too small" );
  static_assert ( s_min < 0 , "std::numeric_limits<float>::max is too small" );
  // ==========================================================================
  static char s_method[] = "evaluate" ;
  // =========================================================================
} //                                           the end of anonnnymous namespace 
// ============================================================================
/*  constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncTree::PyFuncTree 
( const TTree* tree )
  : Ostap::IFuncTree () 
  , m_tree ( tree )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncTree::~PyFuncTree()  {}
// ============================================================================
Ostap::Functions::PyFuncTree*
Ostap::Functions::PyFuncTree::clone( const char* /* name */ ) const
{
  Ostap::Assert ( false ,
                  "Method `clone` must be overriden!"    ,  
                  "Ostap::Functions::PyFuncTree"         ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  // fake ....
  return new PyFuncTree ( *this ) ;
}
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncTree::operator() ( const TTree* t ) const
{ 
  /// redefine the current  tree 
  if ( nullptr != t ) { m_tree = t ; }
  return evaluate () ;
}
// ============================================================================
// function that needs to be redefiend in python 
// ============================================================================
double Ostap::Functions::PyFuncTree::evaluate () const
{
  Ostap::Assert ( false ,
                  "Method `evaluate` must be overriden!" , 
                  "Ostap::Functions::PyFuncTree"         ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return -1000 ;
}
// ============================================================================
/* constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncData::PyFuncData 
( const RooAbsData* data )
  : Ostap::IFuncData () 
  , m_data ( data )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncData::~PyFuncData()  {}
// ============================================================================
Ostap::Functions::PyFuncData*
Ostap::Functions::PyFuncData::clone ( const char* /* name */ ) const
{
  Ostap::Assert ( false ,
                  "Method `clone` must be overriden!"    ,  
                  "Ostap::Functions::PyFuncTree"         ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return new PyFuncData ( *this ) ;
}
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncData::operator() ( const RooAbsData* d ) const
{
  /// redefine the current  tree 
  if ( nullptr != d ) { m_data = d ; }
  //
  return evaluate () ;
}  
// ============================================================================
// function that needs to be redefiend in python 
// ============================================================================
double Ostap::Functions::PyFuncData::evaluate () const
{
  Ostap::Assert ( false ,
                  "Method `evaluate` must be overriden!" , 
                  "Ostap::Functions::PyFuncData"         ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return -1000 ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

