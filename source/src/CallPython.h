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
  /// convert Python to double
  double result_to_double ( PyObject* r , const char* tag  )
  {
    // ========================================================================
    if  ( !r ) 
    {
      PyErr_Print () ;
      Ostap::throwException ( "CallPython:invalid `result'"  , tag , Ostap::StatusCode(500) ) ;
    }
    // ========================================================================
    // float or integer ?
    // ========================================================================
    if      ( PyFloat_Check ( r ) )  // floating value? 
    {
      const double result = PyFloat_AS_DOUBLE ( r ) ;
      Py_DECREF ( r ) ;
      return result ;                                    // RETURN 
    }
    // ========================================================================
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3 
    // ========================================================================
    // Int/long in python2 
    // ========================================================================
    else if ( PyInt_Check ( r ) )   // integer value ? 
    {
      const double result = PyInt_AS_LONG ( r ) ;
      Py_DECREF ( r ) ;
      return result ;                                    //  RETURN
    } 
    //
#else
    // ========================================================================
    // int in python3 
    // ========================================================================
    else if ( PyLong_Check ( r ) )   // integer value ? 
    {
      int overflow = 0 ;
      const long result = PyLong_AsLongAndOverflow ( r , &overflow );
      if      ( -1 == overflow || 1 == overflow ) 
      {
        Py_DECREF ( r ) ;      
        Ostap::throwException ( "CallPython:long overflow"      , tag , Ostap::StatusCode(600) ) ;
      }
      else if ( -1 == result && PyErr_Occurred() ) 
      {
        PyErr_Print();
        Py_DECREF ( r ) ;      
        Ostap::throwException ( "CallPython:invalid conversion" , tag , Ostap::StatusCode(700) ) ;
      }
      Py_DECREF ( r ) ;
      return result ;                                    //  RETURN
    } 
#endif 
    // ========================================================================
    const double result = PyFloat_AsDouble ( r ) ;
    if ( PyErr_Occurred() ) 
    {
      PyErr_Print();
      Py_DECREF ( r ) ;      
      Ostap::throwException ( "CallPython:invalid conversion" , tag , Ostap::StatusCode(800) ) ;
    }
    //
    Py_DECREF ( r ) ; 
    return result   ;                                     // RETURN
  }
  // ==========================================================================
  /// call python Method  and convert resutl to double 
  double call_method ( PyObject* self , char* method )  
  {
    // check arguments
    Ostap::Assert ( self                          ,
                    "CallPython:invalid `object'" , 
                    "call_method"                 ,
                    Ostap::StatusCode(400)        ) ;
    // call Python
    PyObject* result = PyObject_CallMethod ( self , method , nullptr );    
    // error/exception ?
    return result_to_double ( result , method ) ;
  }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // CALLPYTHON_H
// ============================================================================
