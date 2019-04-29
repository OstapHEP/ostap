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
// ============================================================================
// local
// ============================================================================
#include "CallPython.h"
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
/** constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncTree::PyFuncTree 
( PyObject* self , const TTree* tree )
  : Ostap::IFuncTree () 
  , m_tree ( tree )
  , m_self ( self ) 
{
  Py_INCREF( m_self ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncTree::~PyFuncTree() { Py_XDECREF ( m_self ) ; }
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncTree::operator() ( const TTree* t ) const
{
  /// redefine the current  tree 
  if ( nullptr != t ) { m_tree = t ; }
  return call_method ( m_self , s_method ) ;
}
// ============================================================================
/** constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncData::PyFuncData 
( PyObject* self , const RooAbsData* data )
  : Ostap::IFuncData () 
  , m_data ( data )
  , m_self ( self ) 
{
  Py_INCREF( m_self ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncData::~PyFuncData() { Py_XDECREF ( m_self ) ; }
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncData::operator() ( const RooAbsData* d ) const
{
  /// redefine the current  tree 
  if ( nullptr != d ) { m_data = d ; }
  return call_method ( m_self , s_method ) ;
}  
// ============================================================================

// ============================================================================
// The END 
// ============================================================================

