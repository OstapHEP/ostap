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
#include "Exception.h"
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
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
// ============================================================================
/*  constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncTree::PyFuncTree 
( PyObject*    self  , 
  const TTree* tree )
  : Ostap::IFuncTree () 
  , m_tree ( tree )
  , m_self ( self ) 
{ if ( nullptr != m_self ) { Py_INCREF( m_self ) ; } }
// ============================================================================
#else 
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
#endif 
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncTree::~PyFuncTree() 
{
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT 
  Py_XDECREF ( m_self ) ;
#endif 
}
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncTree::operator() ( const TTree* t ) const
{ 
  /// redefine the current  tree 
  if ( nullptr != t ) { m_tree = t ; }
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  //
  Ostap::Assert ( m_self                   , 
                  "self*  points to NULL"  , 
                  "PyFuncTree::operator()" , 
                  Ostap::StatusCode(400)   ) ;
  Ostap::Assert ( m_tree                   , 
                  "TTree* points to NULL"  , 
                  "PyFuncTree::operator()" , 
                  Ostap::StatusCode(401)   ) ;
  return call_method ( m_self , s_method ) ;
  //
#else 
  //
  return evaluate () ;
  //
#endif 
}
// ============================================================================
// function that needs to be redefiend in python 
// ============================================================================
double Ostap::Functions::PyFuncTree::evaluate () const { return -1000 ; }
// ============================================================================



// ============================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
// ============================================================================
/* constructor
 *  @param self python objects
 *  @param tree pointer to the tree
 */
// ============================================================================
Ostap::Functions::PyFuncData::PyFuncData 
( PyObject* self , 
  const RooAbsData* data )
  : Ostap::IFuncData () 
  , m_data ( data )
  , m_self ( self ) 
{ if ( 0 != m_self ) { Py_INCREF( m_self ) ;} }
// ============================================================================
#else 
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
#endif
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::PyFuncData::~PyFuncData() 
{
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  Py_XDECREF ( m_self ) ;
#endif 
}
// ============================================================================
// the basic 
// ============================================================================
double Ostap::Functions::PyFuncData::operator() ( const RooAbsData* d ) const
{
  /// redefine the current  tree 
  if ( nullptr != d ) { m_data = d ; }
  //
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  //
  Ostap::Assert ( m_self                   , 
                  "self*  points to NULL"  , 
                  "PyFuncData::operator()" , 
                  Ostap::StatusCode(400)   ) ;
  Ostap::Assert ( m_data                   , 
                  "RooabsData* points to NULL" , 
                  "PyFuncData::operator()" , 
                  Ostap::StatusCode(401)   ) ;
  return call_method ( m_self , s_method ) ;
  //
#else 
  //
  return evaluate () ;
  //
#endif 
}  
// ============================================================================
// function that needs to be redefiend in python 
// ============================================================================
double Ostap::Functions::PyFuncData::evaluate () const { return -1000 ; }
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================

