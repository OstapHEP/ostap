// ============================================================================
// Include files 
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyCallable.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "TPython.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "CallPython.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::PyCallable
 *  @date 2019-09-25 
 *  @author Vanya Belyaev
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const char s_ERROR1  [] = "Invalid callable"                          ; 
  const char s_ERROR2  [] = "Cannot fill tuple"                         ; 
  const char s_METHOD1 [] = "Ostap::Functions::PyCallable::constructor" ; 
  const char s_METHOD2 [] = "Ostap::Functions::PyCallable::evaluate"    ; 
  // ==========================================================================
}
// ============================================================================
// constructor from the callable object 
// ============================================================================
Ostap::Functions::PyCallable::PyCallable 
( PyObject* callable , const bool /* ok */ ) 
  : m_callable  ( callable ) 
  , m_arguments ( PyTuple_New ( 1 ) )  
{
  //
  Ostap::Assert ( m_callable && PyCallable_Check ( m_callable ) , 
                  s_ERROR1 , s_METHOD1 , 510  ) ;
  //
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::PyCallable::PyCallable 
( const Ostap::Functions::PyCallable&  right ) 
  : m_callable  ( right.m_callable  ) 
  , m_arguments ( PyTuple_New ( 1 ) )  
{
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// Move constructor 
// ============================================================================
Ostap::Functions::PyCallable::PyCallable 
( Ostap::Functions::PyCallable&&  right ) 
  : m_callable  ( right.m_callable  ) 
  , m_arguments ( right.m_arguments ) 
{
  right.m_callable  = nullptr ;
  right.m_arguments = nullptr ;
}
// ============================================================================
// Destructor 
// ============================================================================
Ostap::Functions::PyCallable::~PyCallable() 
{
  if  ( m_callable  ) { Py_DECREF ( m_callable  ) ; m_callable  = nullptr ; }
  if  ( m_arguments ) { Py_DECREF ( m_arguments ) ; m_arguments = nullptr ; }
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
double Ostap::Functions::PyCallable::evaluate ( const  double x ) const 
{
  // 
  if  ( nullptr == m_callable || !PyCallable_Check( m_callable  ) ) 
  {
    PyErr_Print() ;
    Ostap::throwException ( s_ERROR1 , s_METHOD2 , 511 ) ;
  }
  //
  PyObject* px =  PyFloat_FromDouble ( x ) ;
  if ( 0 != PyTuple_SetItem ( m_arguments , 0 , px ) ) 
  {
    PyErr_Print () ;
    Ostap::throwException ( s_ERROR2  , s_METHOD2 , 512 ) ;
  }
  //
  PyObject* result = PyObject_CallObject ( m_callable , m_arguments ) ;
  //
  return result_to_double ( result , s_METHOD2 ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
