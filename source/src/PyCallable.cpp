// ============================================================================
// Include files 
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "TPython.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyCallable.h"
// ============================================================================
// Local
// ============================================================================
#include "CallPython.h"
#include "status_codes.h"
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
  /// defautl (invalid) value 
  const double s_DEFAULT {-1000 } ;   /// default (invalid) value 
  // ==========================================================================
}
// ============================================================================
// constructor from the callable object 
// ============================================================================
Ostap::Functions::PyCallable::PyCallable 
( PyObject* callable , const bool /* ok */ ) 
  : m_callable  ( callable ) 
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
{
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// Move constructor 
// ============================================================================
Ostap::Functions::PyCallable::PyCallable 
( Ostap::Functions::PyCallable&&  right ) 
  : m_callable  ( right.m_callable  ) 
{
  right.m_callable  = nullptr ;
}
// ============================================================================
// Destructor 
// ============================================================================
Ostap::Functions::PyCallable::~PyCallable() 
{
  if  ( m_callable  ) { Py_DECREF ( m_callable  ) ; m_callable  = nullptr ; }
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
      Ostap::Assert ( false     ,
                      s_ERROR1  ,
                      s_METHOD2 ,
                      INVALID_CALLABLE , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* arguments = PyTuple_New ( 1 ) ;
  if ( !arguments )
    {
      PyErr_Print () ;      
      Ostap::Assert ( false ,
                      "Can't create PyTuple"             ,
                      "PyCallable::evaluate"             ,
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* px =  PyFloat_FromDouble ( x ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 0 , px ) ) 
    {
      PyErr_Print () ;
      Py_DECREF ( arguments ) ; arguments = nullptr ;
      Ostap::Assert ( false     ,
                      s_ERROR2  ,
                      s_METHOD2 ,
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ; 
    }
  //
  PyObject* result = PyObject_CallObject ( m_callable , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , s_METHOD2 ) ;
}

// ============================================================================
// constructor from the callable object 
// ============================================================================
Ostap::Functions::PyCallable2::PyCallable2 
( PyObject* callable , const bool /* ok */ ) 
  : m_callable  ( callable ) 
{
  //
  Ostap::Assert ( m_callable && PyCallable_Check ( m_callable ) , 
                  s_ERROR1                               ,
                  s_METHOD1                              ,
                  INVALID_CALLABLE , __FILE__ , __LINE__ ) ; 
  //
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::PyCallable2::PyCallable2 
( const Ostap::Functions::PyCallable2&  right ) 
  : m_callable  ( right.m_callable  ) 
{
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// Move constructor 
// ============================================================================
Ostap::Functions::PyCallable2::PyCallable2 
( Ostap::Functions::PyCallable2&&  right ) 
  : m_callable  ( right.m_callable  ) 
{
  right.m_callable  = nullptr ;
}
// ============================================================================
// Destructor 
// ============================================================================
Ostap::Functions::PyCallable2::~PyCallable2() 
{
  if  ( m_callable  ) { Py_DECREF ( m_callable  ) ; m_callable  = nullptr ; }
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
double Ostap::Functions::PyCallable2::evaluate 
(  const  double x ,
   const  double y ) const 
{
  // 
  if  ( nullptr == m_callable || !PyCallable_Check( m_callable  ) ) 
  {
    PyErr_Print() ;
    Ostap::Assert  ( false                                  ,
                     s_ERROR1                               ,
                     s_METHOD2                              ,
                     INVALID_CALLABLE , __FILE__ , __LINE__ ) ;
    return s_DEFAULT ;
  }
  //
  PyObject* arguments = PyTuple_New ( 2 ) ;
  if ( !arguments )
    {
      PyErr_Print () ;
      Ostap::Assert ( false                               ,
                      "Can't create PyTuple"              ,
                      "PyCallable2::evaluate"             ,
                      ERROR_PYTHON  , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* px =  PyFloat_FromDouble ( x ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 0 , px ) ) 
    {
      PyErr_Print () ;
      Py_DECREF ( arguments ) ; arguments = nullptr ;      
      Ostap::Assert ( false                              ,
                      s_ERROR2                           ,
                      s_METHOD2                          ,
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ; 
    }
  //
  PyObject* py =  PyFloat_FromDouble ( y ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 1 , py ) ) 
    {
      PyErr_Print () ;
      Py_DECREF ( arguments ) ; arguments = nullptr ;
      Ostap::Assert ( false                              ,
                      s_ERROR2                           ,
                      s_METHOD2                          , 
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* result = PyObject_CallObject ( m_callable , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , s_METHOD2 ) ;
}

// ============================================================================
// constructor from the callable object 
// ============================================================================
Ostap::Functions::PyCallable3::PyCallable3 
( PyObject* callable , const bool /* ok */ ) 
  : m_callable  ( callable ) 
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
Ostap::Functions::PyCallable3::PyCallable3 
( const Ostap::Functions::PyCallable3&  right ) 
  : m_callable  ( right.m_callable  ) 
{
  Py_XINCREF ( m_callable ) ;
}
// ============================================================================
// Move constructor 
// ============================================================================
Ostap::Functions::PyCallable3::PyCallable3 
( Ostap::Functions::PyCallable3&&  right ) 
  : m_callable  ( right.m_callable  ) 
{
  right.m_callable  = nullptr ;
}
// ============================================================================
// Destructor 
// ============================================================================
Ostap::Functions::PyCallable3::~PyCallable3() 
{
  if  ( m_callable  ) { Py_DECREF ( m_callable  ) ; m_callable  = nullptr ; }
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
double Ostap::Functions::PyCallable3::evaluate 
(  const  double x ,
   const  double y ,
   const  double z ) const 
{
  // 
  if  ( nullptr == m_callable || !PyCallable_Check( m_callable  ) ) 
    {
      PyErr_Print() ;    
      Ostap::Assert ( false                                  ,
                      s_ERROR1                               ,
                      s_METHOD2                              ,
                      INVALID_CALLABLE , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* arguments = PyTuple_New ( 3 ) ;
  if ( !arguments )
    {
      PyErr_Print () ;
      Ostap::Assert ( false , 
                      "Can't create PyTuple"    ,
                      "PyCallable3::evaluate"   ,
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* px =  PyFloat_FromDouble ( x ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 0 , px ) ) 
  {
    PyErr_Print () ;
    Py_DECREF ( arguments ) ; arguments = nullptr ;    
    Ostap::Assert ( false                              ,
                    s_ERROR2                           ,
                    s_METHOD2                          ,
                    ERROR_PYTHON , __FILE__ , __LINE__ ) ;
    return s_DEFAULT ;
  }
  //
  PyObject* py =  PyFloat_FromDouble ( y ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 1 , py ) ) 
    {
      PyErr_Print () ;
      Py_DECREF ( arguments ) ; arguments = nullptr ;
      Ostap::Assert ( false                              ,
                      s_ERROR2                           ,
                      s_METHOD2                          ,
                      ERROR_PYTHON , __FILE__ , __LINE__ ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* pz =  PyFloat_FromDouble ( z ) ;
  if ( 0 != PyTuple_SetItem ( arguments , 2 , pz ) ) 
  {
    PyErr_Print () ;
    Py_DECREF ( arguments ) ; arguments = nullptr ;
    Ostap::Assert ( false                              ,
                    s_ERROR2                           ,
                    s_METHOD2                          ,
                    ERROR_PYTHON , __FILE__ , __LINE__ ) ;
    return s_DEFAULT ;
  }
  //
  PyObject* result = PyObject_CallObject ( m_callable , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , s_METHOD2 ) ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
