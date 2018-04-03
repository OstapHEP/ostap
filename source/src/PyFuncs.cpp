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
#include "Ostap/StatusCode.h"
#include "Ostap/PyFuncs.h"
// ============================================================================
// local
// ============================================================================
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
  double call_python ( PyObject* self )  
  {
    // check arguments
    if ( !self )
    {
      throw Ostap::Exception( "Python exception: invalid ``self''",
                              "call_python"          ,
                              Ostap::StatusCode(400) ) ;
    }    
    // call Python
    PyObject* r = PyObject_CallMethod ( self, s_method , nullptr );    
    // error/exception ?
    if ( !r ) 
    {
      PyErr_Print();
      throw Ostap::Exception( "Python exception: invalid ``result''" ,
                              "call_python"          ,
                              Ostap::StatusCode(500) ) ;
    }
    // float or integer 
    if      ( PyFloat_Check( r ) )  // floating value? 
    {
      const double result = PyFloat_AS_DOUBLE( r );
      Py_DECREF( r );
      return result ;                                    // RETURN 
    } 
    else if ( PyInt_Check ( r ) )   // integer value ? 
    {
      const double result = PyInt_AS_LONG( r );
      Py_DECREF( r );
      return result ;                                    //  RETURN
    } 
    // ?
    const double result = PyFloat_AsDouble ( r ) ;
    if ( PyErr_Occurred() ) 
    {
      PyErr_Print();
      Py_DECREF ( r );      
      throw Ostap::Exception( "Python exception: invalid conversion" ,
                              "call_python"          ,
                              Ostap::StatusCode(600) ) ;
    }
    //
    Py_DECREF ( r ) ; 
    return result   ;                                     // RETURN
  }
  // ==========================================================================
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
  return call_python ( m_self ) ;
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
  return call_python ( m_self ) ;
}


  
// ======================================================================

// ============================================================================
// The END 
// ============================================================================

