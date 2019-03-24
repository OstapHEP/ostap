// ============================================================================
#ifndef CALLPYTHON_H 
#define CALLPYTHON_H 1
// ============================================================================
// Incldue files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
//   Local
// ============================================================================
#include "Exception.h"
// ============================================================================
namespace 
{
  // ==========================================================================
  /// void -> double  
  double call_python ( PyObject* self , char* method )  
  {
    // check arguments
    Ostap::Assert ( self                          ,
                    "CallPython:invalid ``self''" ,
                    "CallPython:call_python"      ,
                    Ostap::StatusCode(400)        ) ;
    // call Python
    PyObject* r = PyObject_CallMethod ( self , method , nullptr );    
    // error/exception ?
    if ( !r ) 
    {
      PyErr_Print();
      Ostap::throwException ( "CallPython:invalid ``result''" ,
                              "CallPython:call_python"          ,
                              Ostap::StatusCode(500) ) ;
    }
    // float or integer ?
    if      ( PyFloat_Check( r ) )  // floating value? 
    {
      const double result = PyFloat_AS_DOUBLE( r );
      Py_DECREF( r );
      return result ;                                    // RETURN 
    }
    //
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3 
    //
    else if ( PyInt_Check ( r ) )   // integer value ? 
    {
      const double result = PyInt_AS_LONG( r );
      Py_DECREF( r );
      return result ;                                    //  RETURN
    } 
    //
#else
    //
    else if ( PyLong_Check ( r ) )   // integer value ? 
    {
      int overflow = 0 ;
      const long result = PyLong_AsLongAndOverflow ( r , &overflow );
      if      ( -1 == overflow || 1 ==   overflow ) 
      {
        Py_DECREF ( r );      
        throwException ( "CallPython:long overflow" ,
                         "CallPython:call_python"          ,
                         Ostap::StatusCode(600) ) ;
      }
      else if ( -1 == result && PyErr_Occurred() ) 
      {
        PyErr_Print();
        Py_DECREF ( r );      
        throwException ( "CallPython:invalid conversion" ,
                         "CallPython:call_python"          ,
                         Ostap::StatusCode(700) ) ;
      }
      Py_DECREF( r );
      return result ;                                    //  RETURN
    } 
    //
#endif 
    // ?
    const double result = PyFloat_AsDouble ( r ) ;
    if ( PyErr_Occurred() ) 
    {
      PyErr_Print();
      Py_DECREF ( r );      
      throwException ( "CallPython:invalid conversion" ,
                       "CallPython:call_python"          ,
                       Ostap::StatusCode(800) ) ;
    }
    //
    Py_DECREF ( r ) ; 
    return result   ;                                     // RETURN
  }
  // ==========================================================================
}
// ============================================================================
#endif // CALLPYTHON_H
